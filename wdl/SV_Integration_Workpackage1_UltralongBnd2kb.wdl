version 1.0


# From raw single-sample calls to intra-sample merged BND and 2kb+ ultralong
# callsets.
#
# Remark: this workflow is designed to process multiple samples in the same VM
# and to be robust to preemption.
#
workflow SV_Integration_Workpackage1_UltralongBnd2kb {
    input {
        File sv_integration_chunk_tsv
        Boolean has_pav = true
        String region = "all"
        String remote_outdir
        String requester_pays_project = ""
        
        Int min_ultralong_sv_length = 2000
        Int ultralong_collapse_mode = 0
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File reference_agp
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        has_pav: "If false, skip PAV localization/canonization and merge only pbsv and Sniffles calls. The output schema still includes SUPP_PAV, set to zero."
        region: "Only consider VCF records in this genomic region. Set to 'all' to disable."
        remote_outdir: "Without final slash. Where the merged ultralong and BND outputs are stored for each sample."
        requester_pays_project: "Google Cloud project to bill for requester-pays buckets. Leave blank for non-requester-pays buckets."
        min_ultralong_sv_length: "Non-BND calls with ABS(SVLEN) at or above this value are saved in the ultralong stream."
        ultralong_collapse_mode: "0=do not use sequence similarity in truvari collapse; 1=use sequence similarity in truvari collapse."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            has_pav = has_pav,
            region = region,
            remote_outdir = remote_outdir,
            requester_pays_project = requester_pays_project,
            
            min_ultralong_sv_length = min_ultralong_sv_length,
            ultralong_collapse_mode = ultralong_collapse_mode,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            standard_chromosomes_bed = standard_chromosomes_bed,
            reference_agp = reference_agp,
            
            docker_image = docker_image
    }
    
    output {
    }
}


## Memory bottlenecks (measured on a 32GB VM):
#
# CleanRefAltQual            250 MB
# RemoveRefAlt               200 MB
# truvari collapse           100 MB
#
## Multicore bottlenecks (measured on a 16-CPU VM):
# 
# bcftools merge             300 %
# bgzip                      500 %
#
## Truvari collapse ultralong:
#
# --pctseq 0                   2 s
# --pctseq 0.90                2 s to >=1 h
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        Boolean has_pav
        String region
        String remote_outdir
        String requester_pays_project
        
        Int min_ultralong_sv_length
        Int ultralong_collapse_mode
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File reference_agp
        
        String docker_image
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 128
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
        GCLOUD_STORAGE_BILLING_FLAGS=""
        if [ -n "~{requester_pays_project}" ]; then
            GCLOUD_STORAGE_BILLING_FLAGS="--billing-project=~{requester_pays_project}"
        fi
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param
        # $2 A row of `sv_integration_chunk_tsv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
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
            
            if [ ${HAS_PAV} = "true" ]; then
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PAV_VCF_GZ} ./${SAMPLE_ID}_pav.vcf.gz
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PAV_TBI} ./${SAMPLE_ID}_pav.vcf.gz.tbi
            fi
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PBSV_VCF_GZ} ./${SAMPLE_ID}_pbsv.vcf.gz
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PBSV_TBI} ./${SAMPLE_ID}_pbsv.vcf.gz.tbi
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${SNIFFLES_VCF_GZ} ./${SAMPLE_ID}_sniffles.vcf.gz
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${SNIFFLES_TBI} ./${SAMPLE_ID}_sniffles.vcf.gz.tbi
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_pav.vcf.gz* ${SAMPLE_ID}_pbsv.vcf.gz* ${SAMPLE_ID}_sniffles.vcf.gz*
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
            
            ls -lht *.bed 1>&2
            
            # Removing temporary files
            rm -f gaps_unsorted.bed
        }
        
        
        # Puts in canonical form a raw VCF from an SV caller. The procedure
        # creates sorted output files `SAMPLEID_CALLERID_X.vcf.gz`, where X is:
        #
        # sv_ultralong: non-BND records with length >=MIN_ULTRALONG_SV_LENGTH, devoid of
        #               sequence where possible to save space. With the default
        #               MIN_ULTRALONG_SV_LENGTH=2000, this includes records >=2kb.
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
            local MIN_ULTRALONG_SV_LENGTH=$5
            local STANDARD_CHROMOSOMES_BED=$6
            local NOT_GAPS_BED=$7
            
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
            
            # Removing uncalled ALT alleles before SVLEN is treated as Number=A.
            # This prevents malformed multiallelic records with too few SVLEN
            # values from crashing bcftools norm.
            ${TIME_COMMAND} bcftools view --output-type u --min-ac 1 --trim-alt-alleles ${SAMPLE_ID}_${CALLER_ID}_in.vcf | bcftools +fill-tags --output-type v --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf -- -t AC,AN
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
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
            
            # Isolating ultra-long records and discarding shorter records. The
            # default MIN_ULTRALONG_SV_LENGTH=2000 includes ABS(SVLEN) >= 2000.
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_ULTRALONG_SV_LENGTH} --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            # 1. BND VCF -------------------------------------------------------
            
            # 1.1 Sorting 
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # Remark: we do not run the following command, since it seems to
            # destroy BNDs ALTs (example: N]chr5:181473415] ->
            # GNcNNNNNNNNNNNNNN ):
            #
            # bcftools norm --check-ref s --fasta-ref ~{reference_fa}
            # --do-not-normalize
            
            # 1.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --rm-dup exact --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.3 Forcing every record to PASS and adding QUAL, since it is
            # used by `truvari collapse` to select a representation.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 1.4 Removing records not called ALT by this caller.
            ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz.tbi
            
            # 2. Ultralong VCF -------------------------------------------------
            
            # 2.1 Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 2.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            
            # 2.3 Removing sequence (lossless), forcing every record to PASS,
            # and setting QUAL, since it is used by `truvari collapse` to select
            # a representation.
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} RemoveRefAlt ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} ~{reference_fai} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            elif [ ~{ultralong_collapse_mode} -eq 1 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            fi
            
            # 2.4 Removing records not called ALT by this caller.
            ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz.tbi
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
        # This preserves support by caller after truvari collapse.
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
           
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ${GCLOUD_STORAGE_BILLING_FLAGS} ~{remote_outdir}/${SAMPLE_ID}.ultralong_bnd_2kb.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Merging only BND and 2kb+ ultralong records. This intentionally
            # stops before processing shorter SVs, BAM localization, Kanpig,
            # or training-record extraction.
            LocalizeSample ${SAMPLE_ID} ${LINE}
            if [ ${HAS_PAV} = "true" ]; then
                CanonizeVcf ${SAMPLE_ID}_pav.vcf.gz ${SAMPLE_ID}_pav.vcf.gz.tbi ${SAMPLE_ID} pav ~{min_ultralong_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            fi
            CanonizeVcf ${SAMPLE_ID}_pbsv.vcf.gz ${SAMPLE_ID}_pbsv.vcf.gz.tbi ${SAMPLE_ID} pbsv ~{min_ultralong_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            CanonizeVcf ${SAMPLE_ID}_sniffles.vcf.gz ${SAMPLE_ID}_sniffles.vcf.gz.tbi ${SAMPLE_ID} sniffles ~{min_ultralong_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            IntrasampleMerge_ultralong ${SAMPLE_ID}
            IntrasampleMerge_bnd ${SAMPLE_ID}
            
            # Copying SUPP fields to INFO in the BND and ultralong VCFs as well,
            # matching the original WP1 BCF output shape.
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_bnd.vcf.gz b ${SAMPLE_ID}_bnd_supp.bcf
            mv ${SAMPLE_ID}_bnd_supp.bcf ${SAMPLE_ID}_bnd.bcf
            mv ${SAMPLE_ID}_bnd_supp.bcf.csi ${SAMPLE_ID}_bnd.bcf.csi
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_ultralong.vcf.gz b ${SAMPLE_ID}_ultralong_supp.bcf
            mv ${SAMPLE_ID}_ultralong_supp.bcf ${SAMPLE_ID}_ultralong.bcf
            mv ${SAMPLE_ID}_ultralong_supp.bcf.csi ${SAMPLE_ID}_ultralong.bcf.csi
            
            # Uploading
            gcloud storage mv ${GCLOUD_STORAGE_BILLING_FLAGS} ${SAMPLE_ID}_ultralong.'bcf*' ${SAMPLE_ID}_bnd.'bcf*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.ultralong_bnd_2kb.done
            gcloud storage mv ${GCLOUD_STORAGE_BILLING_FLAGS} ${SAMPLE_ID}.ultralong_bnd_2kb.done ~{remote_outdir}/
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
