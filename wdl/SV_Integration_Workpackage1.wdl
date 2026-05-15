version 1.0


# From raw single-sample calls to an intra-sample merged and re-genotyped
# callset, plus a subset of records marked for XGBoost training.
#
# Remark: this workflow is designed to process multiple samples in the same VM
# and to be robust to preemption.
#
workflow SV_Integration_Workpackage1 {
    input {
        File sv_integration_chunk_tsv
        Boolean has_pav = true
        String region = "all"
        String remote_outdir
        
        Int min_sv_length = 20
        Int max_sv_length = 10000
        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        Int ultralong_collapse_mode = 0
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File reference_agp
        File ploidy_bed_female
        File ploidy_bed_male
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        has_pav: "If false, skip PAV localization/canonization and merge only pbsv and Sniffles calls. The output schema still includes SUPP_PAV, set to zero."
        region: "Only consider VCF records in this genomic region. Set to 'all' to disable."
        remote_outdir: "Without final slash. Where the output of intra-sample truvari and kanpig is stored for each sample."
        max_sv_length: "Calls above this length are deemed 'ultralong', are not given to kanpig re-genotyping, and are processed separately."
        ultralong_collapse_mode: "0=do not use sequence similarity in truvari collapse; 1=use sequence similarity in truvari collapse."
        training_resource_vcf_gz: "We assume that the training resource VCF has already been subset to the correct length range upstream."
        training_resource_bed: "Training resource calls can belong only to these regions. Typically a high-confidence dipcall BED, or a BED derived from intersecting multiple dipcall BEDs."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            has_pav = has_pav,
            region = region,
            remote_outdir = remote_outdir,
            
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length,
            kanpig_params_singlesample = kanpig_params_singlesample,
            ultralong_collapse_mode = ultralong_collapse_mode,
            
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_tbi = training_resource_tbi,
            training_resource_bed = training_resource_bed,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            standard_chromosomes_bed = standard_chromosomes_bed,
            autosomes_bed = autosomes_bed,
            reference_agp = reference_agp,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            
            docker_image = docker_image
    }
    
    output {
    }
}


## Memory bottlenecks (measured on a 32GB VM):
#
# FixSymbolicRecords           5 GB
# CleanRefAltQual            250 MB
# RemoveRefAlt               200 MB
# truvari collapse           100 MB
# kanpig                     400 MB
#
## Multicore bottlenecks (measured on a 16-CPU VM):
# 
# bcftools merge             300 %
# bgzip                      500 %
# kanpig                     900 %
#
## Kanpig runtimes:
#
# genome, 16 CPUs, 32GB        2 m
# genome, 6 CPUs, 8GB          4 m
# chr6, 6 CPUs, 8GB           20 s
#
## BAM downloading (measured on a 6 CPUs, 8GB VM):
#
# gsutil -m cp                 8 m
# gcloud storage cp            4 m
#
## Truvari collapse ultralong:
#
# --pctseq 0                   2 s
# --pctseq 0.90                2 s to >=1 h
#
# Truvari bench (measured on a 6 CPUs, 8GB VM):
#
# Sequential                  10 m
# Parallel                     5 m
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        Boolean has_pav
        String region
        String remote_outdir
        
        Int min_sv_length
        Int max_sv_length
        String kanpig_params_singlesample
        Int ultralong_collapse_mode
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File reference_agp
        File ploidy_bed_female
        File ploidy_bed_male
        
        String docker_image
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 256
        Int preemptible_number = 4
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        export RUST_BACKTRACE="full"
        HAS_PAV="~{has_pav}"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # Remark: localizing the BAM could be avoided by running kanpig on the
        # remote BAM (likely only on the windows it considers). This does not
        # work on Google Cloud yet.
        #
        # @param 
        # $2 1=Localizes everything except the BAM. 2=Localizes just the BAM.
        # $3 A row of `sv_integration_chunk_tsv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local MODE=$2
            local LINE=$3
            
            local ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            local ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            local PAV_BED=$(echo ${LINE} | cut -d , -f 5)
            local PAV_TBI=$(echo ${LINE} | cut -d , -f 6)
            local PAV_VCF_GZ=$(echo ${LINE} | cut -d , -f 7)
            local PBSV_TBI=$(echo ${LINE} | cut -d , -f 8)
            local PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 9)
            local SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 10)
            local SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 11)
            if [ ${HAS_PAV} = "false" ]; then
                local N_FIELDS=$(echo ${LINE} | awk -F ',' '{ print NF }')
                if [ ${N_FIELDS} -lt 11 ]; then
                    PBSV_TBI=$(echo ${LINE} | cut -d , -f 5)
                    PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 6)
                    SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 7)
                    SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 8)
                fi
            fi
            
            if [ ${MODE} -eq 2 ]; then
                date 1>&2
                gcloud storage cp ${ALIGNED_BAM} ./${SAMPLE_ID}_aligned.bam
                date 1>&2
                gcloud storage cp ${ALIGNED_BAI} ./${SAMPLE_ID}_aligned.bam.bai
            else
                if [ ${HAS_PAV} = "true" ]; then
                    gcloud storage cp ${PAV_VCF_GZ} ./${SAMPLE_ID}_pav.vcf.gz
                    gcloud storage cp ${PAV_TBI} ./${SAMPLE_ID}_pav.vcf.gz.tbi
                fi
                gcloud storage cp ${PBSV_VCF_GZ} ./${SAMPLE_ID}_pbsv.vcf.gz
                gcloud storage cp ${PBSV_TBI} ./${SAMPLE_ID}_pbsv.vcf.gz.tbi
                gcloud storage cp ${SNIFFLES_VCF_GZ} ./${SAMPLE_ID}_sniffles.vcf.gz
                gcloud storage cp ${SNIFFLES_TBI} ./${SAMPLE_ID}_sniffles.vcf.gz.tbi
            fi
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_aligned.bam* ${SAMPLE_ID}_pav.vcf.gz* ${SAMPLE_ID}_pbsv.vcf.gz* ${SAMPLE_ID}_sniffles.vcf.gz*
        }
        
        
        # Builds a BED file that excludes every gap from the AGP file of
        # the reference.
        #
        function GetReferenceGaps() {
            # Computing non-gap regions
            awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if ( ( $1=="chr1" || $1=="chr2" || $1=="chr3" || $1=="chr4" || $1=="chr5" || $1=="chr6" || $1=="chr7" || $1=="chr8" || $1=="chr9" || $1=="chr10" || \
                           $1=="chr11" || $1=="chr12" || $1=="chr13" || $1=="chr14" || $1=="chr15" || $1=="chr16" || $1=="chr17" || $1=="chr18" || $1=="chr19" || $1=="chr20" || \
                           $1=="chr21" || $1=="chr22" || $1=="chrX" || $1=="chrY" || $1=="chrM" \
                         ) && $5=="N" \
                       ) print $0 \
                 }' ~{reference_agp} > gaps_unsorted.bed
            bedtools sort -i gaps_unsorted.bed -faidx ~{reference_fai} > gaps.bed
            bedtools complement -L -i gaps.bed -g ~{reference_fai} > not_gaps.bed
            
            # Intersecting non-gap regions with the training BED
            bedtools sort -i ~{training_resource_bed} -faidx ~{reference_fai} > training_resource_sorted.bed
            rm -f training_not_gaps_beds.wsv
            local ID="0"
            local ROW
            while read -u 4 ROW; do
                ID=$(( ${ID} + 1 ))
                echo "${ROW}" > ${ID}.bed
                bedtools intersect -a ${ID}.bed -b training_resource_sorted.bed -sorted -g ~{reference_fai} > training_not_gaps_${ID}.bed
                if [ -s training_not_gaps_${ID}.bed ]; then
                    echo "${ID} training_not_gaps_${ID}.bed" >> training_not_gaps_beds.wsv
                else
                    rm -f training_not_gaps_${ID}.bed
                fi
                rm -f ${ID}.bed
            done 4< not_gaps.bed
            ls -lht *.bed 1>&2
            
            # Removing temporary files
            rm -f gaps_unsorted.bed training_resource_sorted.bed
        }
        
        
        # Puts in canonical form a raw VCF from an SV caller. The procedure
        # creates sorted output files `SAMPLEID_CALLERID_X.vcf.gz`, where X is:
        #
        # sv: non-BND records with length in [MIN_SV_LENGTH..MAX_SV_LENGTH], in
        #     canonical form;
        # sv_ultralong: non-BND records with length >MAX_SV_LENGTH, devoid of
        #               sequence where possible to save space;
        # bnd: BND records, in their original form.
        #
        # Remark: the funtion outputs indexed `.vcf.gz` files, since they are
        # needed by `bcftools merge`.
        #
        function CanonizeVcf() {
            local INPUT_VCF_GZ=$1
            local INPUT_TBI=$2
            local SAMPLE_ID=$3
            local CALLER_ID=$4
            local MIN_SV_LENGTH=$5
            local MAX_SV_LENGTH=$6
            local STANDARD_CHROMOSOMES_BED=$7
            local NOT_GAPS_BED=$8
            
            # QUAL is used by truvari collapse to select a representation. We
            # assign values based on which representations we observed to be
            # more accurate in a few test examples.
            if [ ${CALLER_ID} = 'pav' ]; then
                local QUAL="4"
            elif [ ${CALLER_ID} = 'pbsv' ]; then
                local QUAL="3"
            elif [ ${CALLER_ID} = 'sniffles' ]; then
                local QUAL="2"
            fi
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi
            
            # Subsetting to standard chromosomes, or to a given region if any.
            if [ ~{region} != "all" ]; then
                ${TIME_COMMAND} bcftools view --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ~{region} --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            else
                ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Ensuring that SVLEN has the correct type for bcftools norm
            bcftools view --header-only ${SAMPLE_ID}_${CALLER_ID}_in.vcf | sed 's/ID=SVLEN,Number=.,/ID=SVLEN,Number=A,/g' > ${SAMPLE_ID}_${CALLER_ID}_header.txt
            ${TIME_COMMAND} bcftools reheader --header ${SAMPLE_ID}_${CALLER_ID}_header.txt --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Removing SNVs, if any.
            if [ ${CALLER_ID} = 'pav' ]; then
                ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV"' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            fi
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ${SAMPLE_ID}_${CALLER_ID}_in.vcf > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Isolating BNDs
            ${TIME_COMMAND} bcftools filter --include 'SVTYPE="BND"' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf
            ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="BND"' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Isolating ultra-long records and discarding short records
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>'${MAX_SV_LENGTH} --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1. Main VCF ------------------------------------------------------
            
            # 1.1 Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.2 Fixing symbolic records
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G FixSymbolicRecords ${SAMPLE_ID}_${CALLER_ID}_in.vcf ~{reference_fa} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.3 Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.4 Cleaning REF, ALT, QUAL, FILTER. 
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - QUAL is used by truvari collapse to select a representation. 
            #   Symbolic records are NOT given low quality (it was 1 in Phase 1)
            #   since e.g. all DEL records made by Sniffles are symbolic.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.5 Removing END, since its values may be inconsistent and make
            # GATK crash downstream.
            ${TIME_COMMAND} bcftools annotate --remove INFO/END --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.6 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz.tbi
            
            # 2. BND VCF -------------------------------------------------------
            
            # 2.1 Sorting 
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Remark: we do not run the following command, since it seems to
            # destroy BNDs ALTs (example: N]chr5:181473415] ->
            # GNcNNNNNNNNNNNNNN ):
            #
            # bcftools norm --check-ref s --fasta-ref ~{reference_fa}
            # --do-not-normalize
            
            # 2.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --rm-dup exact --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 2.3 Forcing every record to PASS and adding QUAL, since it is
            # used by `truvari collapse` to select a representation.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 2.4 Setting to 0/1 every non-ALT record, otherwise the
            # corresponding truvari collapse SUPP field becomes zero.
            ${TIME_COMMAND} bcftools +setGT --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf -- --target-gt q --include 'GT="ref" || GT="mis"' --new-gt c:0/1 --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz.tbi
            
            # 3. Ultralong VCF -------------------------------------------------
            
            # 3.1 Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 3.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 3.3 Removing sequence (lossless), forcing every record to PASS,
            # and setting QUAL, since it is used by `truvari collapse` to select
            # a representation.
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} RemoveRefAlt ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} ~{reference_fai} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            elif [ ~{ultralong_collapse_mode} -eq 1 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            fi
            
            # 3.4 Setting to 0/1 every non-ALT record, otherwise the
            # corresponding truvari collapse SUPP field becomes zero.
            ${TIME_COMMAND} bcftools +setGT --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf -- --target-gt q --include 'GT="ref" || GT="mis"' --new-gt c:0/1 --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_sv.vcf.gz`,
        # creating an output file `SAMPLEID_sv.vcf.gz`.
        #
        # Remark: the funtion's inputs are indexed `.vcf.gz`, since they are
        # needed by `bcftools merge`. It outputs a `.vcf.gz` since it's needed
        # downstream.
        #
        function IntrasampleMerge_sv() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ${HAS_PAV} = "true" ]; then
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pav_sv.vcf.gz ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz --output ${SAMPLE_ID}_out.vcf
            else
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz --output ${SAMPLE_ID}_out.vcf
            fi
            rm -f ${SAMPLE_ID}_*_sv.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Ensuring that every record has a unique ID, to enable joining by
            # CHROM,POS,ID in downstream calls to `bcftools annotate`. Using
            # CHROM,POS,REF,ALT can make `bcftools annotate` segfault, and the
            # speed of joining by CHROM,POS,ID is independent of SVLEN.
            #
            # Remark: we preserve the original ID just for debugging reasons.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf ; bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_sv.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_ultralong.vcf.gz`,
        # creating an output file `SAMPLEID_ultralong.vcf.gz`.
        #
        # Remark: if `truvari collapse` is run without taking sequence
        # similarity into account, different INS/DUP/CNV sequences of similar
        # length at similar POS may be wrongly collapsed. We tolerate this for
        # speed reasons.
        #
        # Remark: the function's inputs are indexed `.vcf.gz`, since they are
        # needed by `bcftools merge`. It outputs a `.vcf.gz` since it's needed
        # downstream.
        #
        function IntrasampleMerge_ultralong() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ${HAS_PAV} = "true" ]; then
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pav_ultralong.vcf.gz ${SAMPLE_ID}_pbsv_ultralong.vcf.gz ${SAMPLE_ID}_sniffles_ultralong.vcf.gz --output ${SAMPLE_ID}_out.vcf
            else
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pbsv_ultralong.vcf.gz ${SAMPLE_ID}_sniffles_ultralong.vcf.gz --output ${SAMPLE_ID}_out.vcf
            fi
            rm -f ${SAMPLE_ID}_*_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics -any --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Removing SVLEN from symbolic ALTs, in order not to interfere with
            # `truvari collapse`.
            local PCTSEQ_VALUE
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                bcftools view --header-only ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
                ${TIME_COMMAND} bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                    printf("%s",$1); \
                    for (i=2; i<=NF; i++) printf("\t%s",$i); \
                    printf("\n"); \
                }' >> ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
                PCTSEQ_VALUE="0"
            elif [ ~{ultralong_collapse_mode} -eq 1 ]; then
                PCTSEQ_VALUE="0.90"
            fi
            
            bgzip --compress-level 1 ${SAMPLE_ID}_in.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq ${PCTSEQ_VALUE} --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Adding SVLEN back into symbolic ALTs, to avoid overcollapse in the
            # cohort-level bcftools merge downstream.
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            fi
            
            # Ensuring that every record has a unique ID.
            # Remark: we preserve the original ID just for debugging reasons.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf ; bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_ultralong.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_bnd.vcf.gz`,
        # creating an output file `SAMPLEID_bnd.vcf.gz`.
        #
        # Remark: the funtion's inputs are indexed `.vcf.gz`, since they are
        # needed by `bcftools merge`. It outputs a `.vcf.gz` since it's needed
        # downstream.
        #
        function IntrasampleMerge_bnd() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ${HAS_PAV} = "true" ]; then
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pav_bnd.vcf.gz ${SAMPLE_ID}_pbsv_bnd.vcf.gz ${SAMPLE_ID}_sniffles_bnd.vcf.gz --output ${SAMPLE_ID}_out.vcf
            else
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pbsv_bnd.vcf.gz ${SAMPLE_ID}_sniffles_bnd.vcf.gz --output ${SAMPLE_ID}_out.vcf
            fi
            rm -f ${SAMPLE_ID}_*_bnd.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Ensuring that every record has a unique ID
            # Remark: we preserve the original ID just for debugging reasons.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf ; bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_bnd.vcf.gz.tbi
        }
        
        
        # Copies truvari's SUPP field from SAMPLE to three tags in INFO. This is
        # necessary, since kanpig overwrites the SAMPLE column.
        #
        # Remark: the funtion requires an indexed `.vcf.gz` in input, and it
        # outputs an indexed `.vcf.gz` of `bcf`, depending on `OUTPUT_FORMAT`
        # (`z` or `b`).
        #
        function CopySuppToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local OUTPUT_FORMAT=$3
            local OUTPUT_VCF_GZ=$4
            
            if [ ${HAS_PAV} = "true" ]; then
                bcftools query --format '%CHROM\t%POS\t%ID\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    printf("%s",$1); \
                    for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                    if ($4=="0") printf("\t0\t0\t0");
                    else if ($4=="1") printf("\t0\t0\t1");
                    else if ($4=="2") printf("\t0\t1\t0");
                    else if ($4=="3") printf("\t0\t1\t1");
                    else if ($4=="4") printf("\t1\t0\t0");
                    else if ($4=="5") printf("\t1\t0\t1");
                    else if ($4=="6") printf("\t1\t1\t0");
                    else if ($4=="7") printf("\t1\t1\t1");
                    printf("\n"); \
                }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            else
                bcftools query --format '%CHROM\t%POS\t%ID\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    printf("%s",$1); \
                    for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                    if ($4=="0") printf("\t0\t0\t0");
                    else if ($4=="1") printf("\t1\t0\t0");
                    else if ($4=="2") printf("\t0\t1\t0");
                    else if ($4=="3") printf("\t1\t1\t0");
                    printf("\n"); \
                }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            fi
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            # Remark: the order of the callers is now the reverse of the one in
            # which they were bcftools-merged.
            ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,~ID,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV,INFO/SUPP_PAV --output-type ${OUTPUT_FORMAT} ${INPUT_VCF_GZ} --output ${OUTPUT_VCF_GZ}
            if [ ${OUTPUT_FORMAT} = z ]; then
                bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_VCF_GZ}
            elif [ ${OUTPUT_FORMAT} = b ]; then
                bcftools index --threads ${N_THREADS} -f -c ${OUTPUT_VCF_GZ}
            fi
            
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt ${INPUT_VCF_GZ}*
        }
           
        
        # Remark: the function outputs a `.vcf.gz`, since it's needed by the 
        # following steps.
        #
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2
            local INPUT_VCF=$3
            local ALIGNMENTS_BAM=$4

            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            
            # Remark: kanpig needs --sizemin >= --kmer
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --sizemin 10 --sizemax ${INFINITY} --reference ~{reference_fa} --input ${INPUT_VCF} --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF} ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Discarding records that are not marked as present by kanpig
            local N_RECORDS_BEFORE_KANPIG=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            local N_RECORDS_AFTER_KANPIG=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi
            
            # Printing debug information
            local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_KANPIG} / ${N_RECORDS_BEFORE_KANPIG}" | bc )
            echo "${N_RECORDS_AFTER_KANPIG},${N_RECORDS_BEFORE_KANPIG},${PERCENT},Number of records that are marked as ALT by kanpig" > ${SAMPLE_ID}_kanpig.csv
            local N_HETS_IN_AUTOSOMES=$( bcftools query --format '%ID' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            local N_RECORDS_IN_AUTOSOMES=$( bcftools query --format '%ID' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            local PERCENT=$( echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_RECORDS_IN_AUTOSOMES}" | bc )
            echo "${N_HETS_IN_AUTOSOMES},${N_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${SAMPLE_ID}_kanpig.bed.gz
        }
        
        
        # Copies the following kanpig fields from SAMPLE to INFO:
        #
        # KS_1, KS_2, SQ, GQ, DP, AD_NON_ALT, AD_ALL
        #
        # This is necessary, since XGBoost downstream uses only INFO fields.
        #
        # Remark: the funtion requires an indexed `.vcf.gz` in input. It
        # outputs a `.vcf.gz` since it's needed by the following steps.
        #
        function CopyKanpigFieldsToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Creating new header lines
            touch ${SAMPLE_ID}_header.txt
            for FIELD in SQ GQ DP
            do
                bcftools view --header-only ${INPUT_VCF_GZ} | grep ID="${FIELD}," | sed -e 's/FORMAT/INFO/g' >> ${SAMPLE_ID}_header.txt
            done
            echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted to an integer in {0,1,2}.">' >> ${SAMPLE_ID}_header.txt
            
            # Copying fields from FORMAT to INFO. Every record is assumed to
            # have a distinct ID, which is enforced by the steps of the
            # pipeline upstream.
            bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\t%INFO/SUPP_PBSV\t%INFO/SUPP_SNIFFLES\t%INFO/SUPP_PAV\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                KS_1=-1; KS_2=-1; \
                p=0; \
                for (i=1; i<=length($4); i++) { \
                    if (substr($4,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { KS_1=$4; KS_2=$4; } \
                else { KS_1=substr($4,1,p-1); KS_2=substr($4,p+1); } \
                if (KS_1==".") KS_1=-1; \
                if (KS_2==".") KS_2=-1; \
                \
                SQ=$5; \
                if (SQ==".") SQ=-1; \
                \
                GQ=$6; \
                if (GQ==".") GQ=-1; \
                \
                DP=$7; \
                if (DP==".") DP=-1; \
                \
                AD_NON_ALT=-1; AD_ALL=1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { AD_NON_ALT=$8; AD_ALL=$8; } \
                else { AD_NON_ALT=substr($8,1,p-1); AD_ALL=substr($8,p+1); } \
                if (AD_NON_ALT==".") AD_NON_ALT=-1; \
                if (AD_ALL==".") AD_ALL=-1; \
                \
                GT_COUNT=-1; \
                if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|." || $9=="0" || $9==".") GT_COUNT=0; \
                else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|." || $9=="1") GT_COUNT=1; \
                else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
                \
                printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,$10,$11,$12); \
            }' | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_format.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,~ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV --output-type z ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2
            
            rm -f ${INPUT_VCF_GZ}*
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz 
        }
        
        
        cat << 'END' > truvari_bench.sh
#!/bin/bash
SAMPLE_ID=$1
INPUT_VCF_GZ=$2
TRAINING_RESOURCE_VCF_GZ=$3
INFINITY=$4
CHUNK_ID=$5
INCLUDE_BED=$6
${TIME_COMMAND} truvari bench -b ${TRAINING_RESOURCE_VCF_GZ} -c ${INPUT_VCF_GZ} --includebed ${INCLUDE_BED} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ${SAMPLE_ID}_truvari_${CHUNK_ID}/
END
        chmod +x truvari_bench.sh
        
        
        # Extracts every record that has a stringent `truvari bench` match with
        # some records in the resource.
        #
        # Remark: we use `--pick single` to force every resource record to be
        # matched with at most one sample record, which is hopefully the
        # most similar to it. This is because we assume that using a
        # contaminated training set in XGBoost downstream is worse than using a
        # slightly smaller training set. With `--pick multi` e.g. two records in
        # the sample VCF might be matched to the same record in the resource 
        # VCF (probably not good) and vice versa (good).
        #
        # Remark: multiple instances of `truvari bench` are run in parallel
        # using `not_gaps.bed`.
        #
        # Remark: in few anecdotal tests, `--pick multi` seems a bit faster than
        # `--pick single` (4m vs 5m with 6 hyperthreading cores).
        #
        # Remark: both the inputs and the output of the function are indexed
        # `.vcf.gz`, since they are needed by `truvari bench`.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Running in parallel
            ${TIME_COMMAND} xargs --arg-file=training_not_gaps_beds.wsv --max-lines=1 --max-procs=${N_THREADS} ./truvari_bench.sh ${SAMPLE_ID} ${INPUT_VCF_GZ} ~{training_resource_vcf_gz} ${INFINITY}
            
            # Concatenating outputs
            rm -f ${SAMPLE_ID}_outputs.txt
            while read -u 4 ROW; do
                ID=$(echo ${ROW} | cut -d ' ' -f 1)
                echo ${SAMPLE_ID}_truvari_${ID}/tp-comp.vcf.gz >> ${SAMPLE_ID}_outputs.txt
            done 4< training_not_gaps_beds.wsv
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list ${SAMPLE_ID}_outputs.txt --output-type z --output ${SAMPLE_ID}_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_training.vcf.gz
            
            # Removing temporary files
            rm -rf ${SAMPLE_ID}_script.sh ${SAMPLE_ID}_outputs.txt ./${SAMPLE_ID}_truvari_*/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        ~{docker_dir}/kanpig --version 1>&2
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Merging
            LocalizeSample ${SAMPLE_ID} 1 ${LINE}
            if [ ${HAS_PAV} = "true" ]; then
                CanonizeVcf ${SAMPLE_ID}_pav.vcf.gz ${SAMPLE_ID}_pav.vcf.gz.tbi ${SAMPLE_ID} pav ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            fi
            CanonizeVcf ${SAMPLE_ID}_pbsv.vcf.gz ${SAMPLE_ID}_pbsv.vcf.gz.tbi ${SAMPLE_ID} pbsv ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            CanonizeVcf ${SAMPLE_ID}_sniffles.vcf.gz ${SAMPLE_ID}_sniffles.vcf.gz.tbi ${SAMPLE_ID} sniffles ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            IntrasampleMerge_sv ${SAMPLE_ID}
            IntrasampleMerge_ultralong ${SAMPLE_ID}
            IntrasampleMerge_bnd ${SAMPLE_ID}
            
            # Genotyping and marking training records
            LocalizeSample ${SAMPLE_ID} 2 ${LINE}
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_sv.vcf.gz z ${SAMPLE_ID}_sv_supp.vcf.gz
            Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_sv_supp.vcf.gz ${SAMPLE_ID}_aligned.bam
            CopyKanpigFieldsToInfo ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            
            # Copying SUPP fields to INFO in the BND and ultralong VCFs as well,
            # just for uniformity. The original SUPP in FORMAT remains there,
            # and since these VCFs won't be re-genotyped with kanpig, it will be
            # correctly preserved by cohort-level truvari collapse.
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_bnd.vcf.gz b ${SAMPLE_ID}_bnd_supp.bcf
            mv ${SAMPLE_ID}_bnd_supp.bcf ${SAMPLE_ID}_bnd.bcf
            mv ${SAMPLE_ID}_bnd_supp.bcf.csi ${SAMPLE_ID}_bnd.bcf.csi
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_ultralong.vcf.gz b ${SAMPLE_ID}_ultralong_supp.bcf
            mv ${SAMPLE_ID}_ultralong_supp.bcf ${SAMPLE_ID}_ultralong.bcf
            mv ${SAMPLE_ID}_ultralong_supp.bcf.csi ${SAMPLE_ID}_ultralong.bcf.csi
            
            # Uploading
            gcloud storage mv ${SAMPLE_ID}_kanpig.vcf.'gz*' ${SAMPLE_ID}_kanpig.bed.gz ${SAMPLE_ID}_kanpig.csv ${SAMPLE_ID}_training.vcf.'gz*' ${SAMPLE_ID}_ultralong.'bcf*' ${SAMPLE_ID}_bnd.'bcf*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
