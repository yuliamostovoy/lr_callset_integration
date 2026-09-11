version 1.0


# Annotates BND VCFs with multiple BAM- and genotyper-derived features.
#
# Remark: for BNDs, `feature_extraction.py` outputs empty values for features:
#
# FEX_DEPTH_MAD, FEX_READ_LEN_MED, FEX_STRAND_BIAS
#
# These (along with SVLEN and SUPP_PAV) should not be used for scoring (PAV does
# not emit BND records).
#
# Remark: BNDs are not canonized at this stage. Every BND record from each 
# sample is kept intact and simply annotated.
#
workflow SV_Integration_BndAnnotate {
    input {
        File chunk_tsv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai

        String min_mapq = "0,60"
        
        Int custom_breakpoint_window_bp = 500
        Int custom_breakpoint_window_slack_bp = 150
        Int custom_min_clip_length = 200
        String custom_min_indel_length = "500,1000,5000"
        Int custom_adjacency_slack_bp = 300

        File feature_extraction_py
        Int use_cutefc = 0

        File tr_bed
        File segdup_bed
        File gc_content_bed
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 3
    }
    parameter_meta {
        chunk_tsv: "Format: `ID,mean_coverage,bai,bam`."
        tr_bed: "From: https://github.com/PacificBiosciences/pbsv/tree/master/annotations"
        segdup_bed: "From: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/"
        gc_content_bed: "From: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/"
        min_mapq: "Comma-separated list of MAPQs. Custom annotations only use alignments with MAPQ>=this (except for avg. MAPQ)."
        custom_min_indel_length: "Comma-separated list of min CIGAR INDEL lengths. Custom annotations only use CIGAR INDELs >=this."
    }
    
    call Impl {
        input:
            chunk_tsv = chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,

            min_mapq = min_mapq,

            custom_breakpoint_window_bp = custom_breakpoint_window_bp,
            custom_breakpoint_window_slack_bp = custom_breakpoint_window_slack_bp,
            custom_min_clip_length = custom_min_clip_length,
            custom_min_indel_length = custom_min_indel_length,
            custom_adjacency_slack_bp = custom_adjacency_slack_bp,
    
            feature_extraction_py = feature_extraction_py,
            use_cutefc = use_cutefc,

            tr_bed = tr_bed,
            segdup_bed = segdup_bed,
            gc_content_bed = gc_content_bed,

            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    output {
    }
}


# Performance on 15x HPRC+HGSVC samples using a 4-core, 8GB VM.
#
# TOOL                                                CPU     RAM     TIME
# BAM download                                                          5m
#
# samtools bedcov (4 threads)                        350%    900M       5m
# annotate_mapq_secondary.sh                         400%     15M      30s
# bcftools annotate                                  100%     15M      50s
# java BndGetBins                                    200%     50M       1s
# java BndCreateBedcovAnnotations                    200%     50M       1s
# annotate_clipped_alignments_1.sh                   400%    600M       2m
# annotate_clipped_alignments_2.sh                   400%     50M       2m
#
# feature_extraction.py (1 thread)                   100%      2G       1m
# cutefc (2 threads)                                  50%    1.5G      30m
#
task Impl {
    input {
        File chunk_tsv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai

        String min_mapq
        
        Int custom_breakpoint_window_bp
        Int custom_breakpoint_window_slack_bp
        Int custom_min_clip_length
        String custom_min_indel_length
        Int custom_adjacency_slack_bp

        File feature_extraction_py
        Int use_cutefc

        File tr_bed
        File segdup_bed
        File gc_content_bed
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 50
        Int preemptible_number
    }
    parameter_meta {
        n_cpu: "4 is good enough, since only custom annotations use all available threads, and they are not the bottleneck."
        ram_size_gb: "16GB is needed by `feature_extraction.py`."
        disk_size_gb: "50GB is needed by the BAMs."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param 
        # $2 A row of `chunk_tsv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_bnd.bcf ./${SAMPLE_ID}.bcf
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_bnd.bcf.csi ./${SAMPLE_ID}.bcf.csi

            # Converting to .vcf.gz for downstream tools
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}.vcf.gz
            rm -f ${SAMPLE_ID}.bcf*

            local N_RECORDS=$( bcftools index --nrecords ${SAMPLE_ID}.vcf.gz.tbi )
            if [ "${N_RECORDS}" != "0" ]; then
                local ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
                local ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
                ${TIME_COMMAND} gcloud storage cp ${ALIGNED_BAM} ./${SAMPLE_ID}.bam
                gcloud storage cp ${ALIGNED_BAI} ./${SAMPLE_ID}.bam.bai    
            fi
        }
        
        
        # Deletes all files related to the given sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}*.bam* ${SAMPLE_ID}*.bcf* ${SAMPLE_ID}*.vcf* ${SAMPLE_ID}*.tsv* ${SAMPLE_ID}*.csv* ${SAMPLE_ID}*.txt
        }
        
        
        # Ensures that the VCF is correctly formatted.
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            gunzip -c ${INPUT_VCF_GZ} > ${SAMPLE_ID}_in.vcf
            rm -f ${INPUT_VCF_GZ}*
            
            # 1. Cleaning REF, ALT, QUAL, FILTER.
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQualBnd ${SAMPLE_ID}_in.vcf ${DEFAULT_QUAL} > ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 2. Making sure IDs will be distinct at the inter-sample level 
            # (they are already distinct at the intra-sample level, thanks to
            # the steps upstream).
            ${TIME_COMMAND} bcftools annotate --set-id ${SAMPLE_ID}'_%ID' --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            mv ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}_canonized.vcf
        }
        
        
        
        
        # ------------------------- Custom annotations -------------------------        
        
        # Given an input VCF containing only BND records, the procedure
        # annotates each record with the coverage of windows around POS and ALT
        # extracted from a BAM.
        #
        # Remark: `samtools bedcov` skips reads with any of the following flags
        # set: UNMAP, SECONDARY, QCFAIL, DUP
        #
        function AnnotateCoverageBins() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local MEAN_COVERAGE=$5
            local MIN_MAPQ=$6

            ${TIME_COMMAND} java -cp ~{docker_dir} BndGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} 0 > ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} samtools bedcov -j --min-MQ ${MIN_MAPQ} ${SAMPLE_ID}_bins.bed ${INPUT_BAM} > ${SAMPLE_ID}_counts.bed
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} java -cp ~{docker_dir} BndCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed ${BREAKPOINT_WINDOW_BP} ${MEAN_COVERAGE} | sort -k 1,1 -k 2,2n > ${SAMPLE_ID}_tags.tsv
            rm -f ${SAMPLE_ID}_counts.bed
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ${INPUT_VCF} | sort -k 5,5 > ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 5 -2 1 ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_tags.tsv | sort -k 1,1 -k 6,6 | paste - - - - | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $7,$14,$21,$28); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_tags.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POS_0_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin before the breakpoint at POS">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POS_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin after the breakpoint at POS">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POS_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin before the breakpoint at ALT">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POS_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin after the breakpoint at ALT">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/BIN_POS_0_'${MIN_MAPQ}',INFO/BIN_POS_1_'${MIN_MAPQ}',INFO/BIN_POS_2_'${MIN_MAPQ}',INFO/BIN_POS_3_'${MIN_MAPQ}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }


        cat << 'END' > annotate_mapq_secondary.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
CHROM=$2
START=$3
END=$4
RECORD_ID=$5
BIN_ID=$6

samtools view --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} | awk '{ sum+=$5; count++ } END { printf("%.4f\n", count>0?sum/count:0) }' > ${RECORD_ID}_${BIN_ID}_mapq.txt
N_SECONDARY=$(samtools view --count ${INPUT_BAM} --require-flags 256 ${CHROM}:${START}-${END})
N_TOTAL=$(samtools view --count ${INPUT_BAM} ${CHROM}:${START}-${END})
if [ ${N_TOTAL} -eq 0 ]; then
    RATIO="0"
else
    RATIO=$(printf "%.4f\n" $(echo "scale=4; ${N_SECONDARY} / ${N_TOTAL}" | bc))
fi
echo "${RATIO}" > ${RECORD_ID}_${BIN_ID}_secondary.txt
END
        chmod +x annotate_mapq_secondary.sh
        
        
        # Given an input VCF containing only BND records, the procedure
        # annotates each record with the average MAPQ and the number of
        # secondary alignments (=repeat-induced multi-mappings) over each of its
        # breakpoints, extracted from a BAM.
        #
        # Remark: we do not collect the number of supplementary alignments,
        # since more specific counts are already captured by
        # `AnnotateClippedAlignments()`.
        #
        function AnnotateMapqSecondary() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
    
            ${TIME_COMMAND} java -cp ~{docker_dir} BndGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} 0 | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            rm -f *_mapq.txt *_secondary.txt
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_counts.tsv
            while read -u 6 ID || [ -n "${ID}" ]; do
                local MAPQ_POINT_0=$(cat ${ID}_0_mapq.txt)
                local MAPQ_POINT_1=$(cat ${ID}_1_mapq.txt)
                local MAPQ_POINT_2=$(cat ${ID}_2_mapq.txt)
                local MAPQ_POINT_3=$(cat ${ID}_3_mapq.txt)
                local SECONDARY_POINT_0=$(cat ${ID}_0_secondary.txt)
                local SECONDARY_POINT_1=$(cat ${ID}_1_secondary.txt)
                local SECONDARY_POINT_2=$(cat ${ID}_2_secondary.txt)
                local SECONDARY_POINT_3=$(cat ${ID}_3_secondary.txt)
                echo -e "${ID}\t${MAPQ_POINT_0}\t${SECONDARY_POINT_0}\t${MAPQ_POINT_1}\t${SECONDARY_POINT_1}\t${MAPQ_POINT_2}\t${SECONDARY_POINT_2}\t${MAPQ_POINT_3}\t${SECONDARY_POINT_3}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 6< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ${INPUT_VCF} | sort -k 5,5 > ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 5 -2 1 ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $6,$7,$8,$9,$10,$11,$12,$13); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POINT_MAPQ_0,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_MAPQ_1,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_MAPQ_2,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_MAPQ_3,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY_0,Number=1,Type=Float,Description="Breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY_1,Number=1,Type=Float,Description="Breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY_2,Number=1,Type=Float,Description="Breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY_3,Number=1,Type=Float,Description="Breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/BIN_POINT_MAPQ_0,INFO/BIN_POINT_SECONDARY_0,INFO/BIN_POINT_MAPQ_1,INFO/BIN_POINT_SECONDARY_1,INFO/BIN_POINT_MAPQ_2,INFO/BIN_POINT_SECONDARY_2,INFO/BIN_POINT_MAPQ_3,INFO/BIN_POINT_SECONDARY_3'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
     

        cat << 'END' > annotate_clipped_alignments_1.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
CLASSPATH=$2
MIN_CLIP_LENGTH=$3
MIN_INDEL_LENGTH=$4
MIN_MAPQ=$5
CHROM=$6
START=$7
END=$8
RECORD_ID=$9
BIN_ID=${10}

samtools view --min-MQ ${MIN_MAPQ} --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} > ${RECORD_ID}_${BIN_ID}.sam
echo $(wc -l < ${RECORD_ID}_${BIN_ID}.sam) > ${RECORD_ID}_${BIN_ID}_n_alignments.txt
java -cp ${CLASSPATH} UltralongIntervalGetClips ${RECORD_ID}_${BIN_ID}.sam ${START} ${END} ${MIN_CLIP_LENGTH} ${MIN_INDEL_LENGTH} ${RECORD_ID}_${BIN_ID}
rm -f ${RECORD_ID}_${BIN_ID}.sam
sort -k 1,1 ${RECORD_ID}_${BIN_ID}_leftmaximal.txt > ${RECORD_ID}_${BIN_ID}_leftmaximal_sorted.txt
sort -k 1,1 ${RECORD_ID}_${BIN_ID}_rightmaximal.txt > ${RECORD_ID}_${BIN_ID}_rightmaximal_sorted.txt
rm -f ${RECORD_ID}_${BIN_ID}_leftmaximal.txt ${RECORD_ID}_${BIN_ID}_rightmaximal.txt
END
        chmod +x annotate_clipped_alignments_1.sh

        
        cat << 'END' > annotate_clipped_alignments_2.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
ADJACENCY_SLACK_BP=$2
MEAN_COVERAGE=$3
NORMALIZATION_MODE=$4
MIN_INDEL_LENGTH=$5
RECORD_ID=$6

# Counting
if [ "${NORMALIZATION_MODE}" -eq 0 ]; then
    N_ALIGNMENTS_0=${MEAN_COVERAGE}
    N_ALIGNMENTS_1=${MEAN_COVERAGE}
    N_ALIGNMENTS_2=${MEAN_COVERAGE}
    N_ALIGNMENTS_3=${MEAN_COVERAGE}
    NORMALIZATION_FACTOR_02=${MEAN_COVERAGE}
    NORMALIZATION_FACTOR_03=${MEAN_COVERAGE}
    NORMALIZATION_FACTOR_12=${MEAN_COVERAGE}
    NORMALIZATION_FACTOR_13=${MEAN_COVERAGE}
else
    N_ALIGNMENTS_0=$(cat ${RECORD_ID}_0_n_alignments.txt)
    N_ALIGNMENTS_1=$(cat ${RECORD_ID}_1_n_alignments.txt)
    N_ALIGNMENTS_2=$(cat ${RECORD_ID}_2_n_alignments.txt)
    N_ALIGNMENTS_3=$(cat ${RECORD_ID}_3_n_alignments.txt)
    NORMALIZATION_FACTOR_02=$(echo "scale=4; ( ${N_ALIGNMENTS_0} + ${N_ALIGNMENTS_2} ) / 2" | bc)
    NORMALIZATION_FACTOR_03=$(echo "scale=4; ( ${N_ALIGNMENTS_0} + ${N_ALIGNMENTS_3} ) / 2" | bc)
    NORMALIZATION_FACTOR_12=$(echo "scale=4; ( ${N_ALIGNMENTS_1} + ${N_ALIGNMENTS_2} ) / 2" | bc)
    NORMALIZATION_FACTOR_13=$(echo "scale=4; ( ${N_ALIGNMENTS_1} + ${N_ALIGNMENTS_3} ) / 2" | bc)
fi
rm -f ${RECORD_ID}_*_n_alignments.txt

_0_R=$(wc -l < ${RECORD_ID}_0_rightmaximal_sorted.txt)
_1_L=$(wc -l < ${RECORD_ID}_1_leftmaximal_sorted.txt)
_2_R=$(wc -l < ${RECORD_ID}_2_rightmaximal_sorted.txt)
_3_L=$(wc -l < ${RECORD_ID}_3_leftmaximal_sorted.txt)

_0_R_2_R=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_0_rightmaximal_sorted.txt ${_0_R} 0 ${RECORD_ID}_2_rightmaximal_sorted.txt ${_2_R} 0 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR_02} | tr ',' '\t')
_0_R_3_L=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_0_rightmaximal_sorted.txt ${_0_R} 0 ${RECORD_ID}_3_leftmaximal_sorted.txt ${_3_L} 1 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR_03} | tr ',' '\t')
_1_L_2_R=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_1_leftmaximal_sorted.txt ${_1_L} 1 ${RECORD_ID}_2_rightmaximal_sorted.txt ${_2_R} 0 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR_12} | tr ',' '\t')
_1_L_3_L=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_1_leftmaximal_sorted.txt ${_1_L} 1 ${RECORD_ID}_3_leftmaximal_sorted.txt ${_3_L} 1 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR_13} | tr ',' '\t')

rm -f ${RECORD_ID}_*maximal_sorted.txt

# Normalizing
if [ $(echo "${N_ALIGNMENTS_0} == 0" | bc) -eq 1 ]; then
    _0_R="0"
else
    _0_R=$(printf "%.4f\n" $(echo "scale=4; ${_0_R} / ${N_ALIGNMENTS_0}" | bc))
fi
if [ $(echo "${N_ALIGNMENTS_1} == 0" | bc) -eq 1 ]; then
    _1_L="0"
else
    _1_L=$(printf "%.4f\n" $(echo "scale=4; ${_1_L} / ${N_ALIGNMENTS_1}" | bc))
fi
if [ $(echo "${N_ALIGNMENTS_2} == 0" | bc) -eq 1 ]; then
    _2_R="0"
else
    _2_R=$(printf "%.4f\n" $(echo "scale=4; ${_2_R} / ${N_ALIGNMENTS_2}" | bc))
fi
if [ $(echo "${N_ALIGNMENTS_3} == 0" | bc) -eq 1 ]; then
    _3_L="0"
else
    _3_L=$(printf "%.4f\n" $(echo "scale=4; ${_3_L} / ${N_ALIGNMENTS_3}" | bc))
fi
OUTPUT="${RECORD_ID}\t${_0_R}\t${_1_L}\t${_2_R}\t${_3_L}\t${_0_R_2_R}\t${_0_R_3_L}\t${_1_L_2_R}\t${_1_L_3_L}"
i=0
for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
    i=$(( ${i} + 1 )); _0_INS=$(cut -f ${i} ${RECORD_ID}_0_indel.txt); _1_INS=$(cut -f ${i} ${RECORD_ID}_1_indel.txt); _2_INS=$(cut -f ${i} ${RECORD_ID}_2_indel.txt); _3_INS=$(cut -f ${i} ${RECORD_ID}_3_indel.txt)
    i=$(( ${i} + 1 )); _0_DEL_START=$(cut -f ${i} ${RECORD_ID}_0_indel.txt); _1_DEL_START=$(cut -f ${i} ${RECORD_ID}_1_indel.txt); _2_DEL_START=$(cut -f ${i} ${RECORD_ID}_2_indel.txt); _3_DEL_START=$(cut -f ${i} ${RECORD_ID}_3_indel.txt)
    i=$(( ${i} + 1 )); _0_DEL_END=$(cut -f ${i} ${RECORD_ID}_0_indel.txt); _1_DEL_END=$(cut -f ${i} ${RECORD_ID}_1_indel.txt); _2_DEL_END=$(cut -f ${i} ${RECORD_ID}_2_indel.txt); _3_DEL_END=$(cut -f ${i} ${RECORD_ID}_3_indel.txt)
    if [ $(echo "${N_ALIGNMENTS_0} == 0" | bc) -eq 1 ]; then
        _0_INS="0"; _0_DEL_START="0"; _0_DEL_END="0" 
    else
        _0_INS=$(printf "%.4f\n" $(echo "scale=4; ${_0_INS} / ${N_ALIGNMENTS_0}" | bc))
        _0_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${_0_DEL_START} / ${N_ALIGNMENTS_0}" | bc))
        _0_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${_0_DEL_END} / ${N_ALIGNMENTS_0}" | bc))
    fi
    if [ $(echo "${N_ALIGNMENTS_1} == 0" | bc) -eq 1 ]; then
        _1_INS="0"; _1_DEL_START="0"; _1_DEL_END="0"
    else
        _1_INS=$(printf "%.4f\n" $(echo "scale=4; ${_1_INS} / ${N_ALIGNMENTS_1}" | bc))
        _1_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${_1_DEL_START} / ${N_ALIGNMENTS_1}" | bc))
        _1_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${_1_DEL_END} / ${N_ALIGNMENTS_1}" | bc))
    fi
    if [ $(echo "${N_ALIGNMENTS_2} == 0" | bc) -eq 1 ]; then
        _2_INS="0"; _2_DEL_START="0"; _2_DEL_END="0"
    else
        _2_INS=$(printf "%.4f\n" $(echo "scale=4; ${_2_INS} / ${N_ALIGNMENTS_2}" | bc))
        _2_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${_2_DEL_START} / ${N_ALIGNMENTS_2}" | bc))
        _2_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${_2_DEL_END} / ${N_ALIGNMENTS_2}" | bc))
    fi
    if [ $(echo "${N_ALIGNMENTS_3} == 0" | bc) -eq 1 ]; then
        _3_INS="0"; _3_DEL_START="0"; _3_DEL_END="0"
    else
        _3_INS=$(printf "%.4f\n" $(echo "scale=4; ${_3_INS} / ${N_ALIGNMENTS_3}" | bc))
        _3_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${_3_DEL_START} / ${N_ALIGNMENTS_3}" | bc))
        _3_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${_3_DEL_END} / ${N_ALIGNMENTS_3}" | bc))
    fi
    OUTPUT="${OUTPUT}\t${_0_INS}\t${_0_DEL_START}\t${_0_DEL_END}\t${_1_INS}\t${_1_DEL_START}\t${_1_DEL_END}\t${_2_INS}\t${_2_DEL_START}\t${_2_DEL_END}\t${_3_INS}\t${_3_DEL_START}\t${_3_DEL_END}"
done
rm -f ${RECORD_ID}_*indel.txt
echo -e "${OUTPUT}" > ${RECORD_ID}_counts.txt
END
        chmod +x annotate_clipped_alignments_2.sh


        # Given an input VCF containing only BND records, the procedure
        # annotates it with clipped alignment measures extracted from a BAM.
        #
        function AnnotateClippedAlignments() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local ADJACENCY_SLACK_BP=$5
            local MIN_CLIP_LENGTH=$6
            local MIN_INDEL_LENGTH=$7
            local CLIP_SLACK_BP=$8
            local MEAN_COVERAGE=$9
            local MIN_MAPQ=${10}
            local NORMALIZATION_MODE=${11}
    
            ${TIME_COMMAND} java -cp ~{docker_dir} BndGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} ${CLIP_SLACK_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ${MIN_CLIP_LENGTH} ${MIN_INDEL_LENGTH} ${MIN_MAPQ}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
            rm -f *_counts.txt
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2.sh ~{docker_dir} ${ADJACENCY_SLACK_BP} ${MEAN_COVERAGE} ${NORMALIZATION_MODE} ${MIN_INDEL_LENGTH}
            rm -f ${SAMPLE_ID}_variantID.txt
            cat *_counts.txt | sort -k 1,1 > ${SAMPLE_ID}_counts.tsv
            rm -f *_counts.txt
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ${INPUT_VCF} | sort -k 5,5 > ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 5 -2 1 ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s\t%d\t%s\t%s\t%s",$2,$3,$4,$5,$1); \
                for (i=6; i<=NF; i++) printf("\t%.4f",$i); \
                printf("\n"); \
            }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz

            echo '##INFO=<ID=C0R_'${MIN_MAPQ}',Number=1,Type=Float,Description="First breakpoint, left bin: number of right-clipped alignments.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_'${MIN_MAPQ}',Number=1,Type=Float,Description="First breakpoint, right bin: number of left-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C2R_'${MIN_MAPQ}',Number=1,Type=Float,Description="Second breakpoint, left bin: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C3L_'${MIN_MAPQ}',Number=1,Type=Float,Description="Second breakpoint, right bin: number of left-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            
            echo '##INFO=<ID=C0R_C2R_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C0R_C2R_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C0R_C2R_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C0R_C2R_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt

            echo '##INFO=<ID=C0R_C3L_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C0R_C3L_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C0R_C3L_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C0R_C3L_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt

            echo '##INFO=<ID=C1L_C2R_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_C2R_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_C2R_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_C2R_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt

            echo '##INFO=<ID=C1L_C3L_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_C3L_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_C3L_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=C1L_C3L_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of clipped reads across bins.">' >> ${SAMPLE_ID}_header.txt

            for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
                echo '##INFO=<ID=C0_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS of length >='${LENGTH}' in the left bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C0_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start of length >='${LENGTH}' in the left bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C0_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end of length >='${LENGTH}' in the left bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C1_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS of length >='${LENGTH}' in the right bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C1_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start of length >='${LENGTH}' in the right bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C1_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end of length >='${LENGTH}' in the right bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C2_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS of length >='${LENGTH}' in the left bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C2_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start of length >='${LENGTH}' in the left bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C2_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end of length >='${LENGTH}' in the left bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C3_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS of length >='${LENGTH}' in the right bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C3_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start of length >='${LENGTH}' in the right bin.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=C3_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end of length >='${LENGTH}' in the right bin.">' >> ${SAMPLE_ID}_header.txt
            done
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/C0R_'${MIN_MAPQ}',INFO/C1L_'${MIN_MAPQ}',INFO/C2R_'${MIN_MAPQ}',INFO/C3L_'${MIN_MAPQ}',INFO/C0R_C2R_1_'${MIN_MAPQ}',INFO/C0R_C2R_2_'${MIN_MAPQ}',INFO/C0R_C2R_3_'${MIN_MAPQ}',INFO/C0R_C2R_4_'${MIN_MAPQ}',INFO/C0R_C3L_1_'${MIN_MAPQ}',INFO/C0R_C3L_2_'${MIN_MAPQ}',INFO/C0R_C3L_3_'${MIN_MAPQ}',INFO/C0R_C3L_4_'${MIN_MAPQ}',INFO/C1L_C2R_1_'${MIN_MAPQ}',INFO/C1L_C2R_2_'${MIN_MAPQ}',INFO/C1L_C2R_3_'${MIN_MAPQ}',INFO/C1L_C2R_4_'${MIN_MAPQ}',INFO/C1L_C3L_1_'${MIN_MAPQ}',INFO/C1L_C3L_2_'${MIN_MAPQ}',INFO/C1L_C3L_3_'${MIN_MAPQ}',INFO/C1L_C3L_4_'${MIN_MAPQ}
            for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
                COLUMNS=${COLUMNS}',INFO/C0_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/C0_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/C0_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',INFO/C1_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/C1_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/C1_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',INFO/C2_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/C2_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/C2_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',INFO/C3_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/C3_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/C3_DEL_END_'${MIN_MAPQ}'_'${LENGTH}
            done
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }

        
        function AnnotateCustom() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local MEAN_COVERAGE=$3
            local MIN_MAPQ=$4
            local NORMALIZATION_MODE=$5

            local MAPQS=$( echo ${MIN_MAPQ} | tr ',' ' ' )
            
            mv ${INPUT_VCF} ${SAMPLE_ID}_in.vcf
            
            for MAPQ in ${MAPQS}; do
                AnnotateCoverageBins ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ${MEAN_COVERAGE} ${MAPQ}
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            done
            AnnotateMapqSecondary ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            for MAPQ in ${MAPQS}; do
                AnnotateClippedAlignments ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ~{custom_adjacency_slack_bp} ~{custom_min_clip_length} ~{custom_min_indel_length} ~{custom_breakpoint_window_slack_bp} ${MEAN_COVERAGE} ${MAPQ} ${NORMALIZATION_MODE}
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            done

            mv ${SAMPLE_ID}_in.vcf ${INPUT_VCF}
        }


        # Adds field `INFO/GT_COUNT` to the input VCF, which is overwritten.
        #
        function AnnotateGtCount() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2

            bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT]\n' ${INPUT_VCF} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($6=="0/0" || $6=="0|0" || $6=="./."  || $6==".|." || $6=="./0" || $6==".|0" || $6=="0/." || $6=="0|." || $6=="0" || $6==".") GT_COUNT=0; \
                else if ($6=="0/1" || $6=="0|1" || $6=="1/0" || $6=="1|0" || $6=="./1" || $6==".|1" || $6=="1/." || $6=="1|." || $6=="1") GT_COUNT=1; \
                else if ($6=="1/1" || $6=="1|1") GT_COUNT=2; \
                printf("%s\t%d\t%s\t%s\t%s\t%d\n",$1,$2,$3,$4,$5,  GT_COUNT); \
            }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="Original GT from the callers, converted to an integer in {0,1,2}.">' > ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/GT_COUNT'
            bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt ${INPUT_VCF}
            mv ${SAMPLE_ID}_annotated.vcf ${INPUT_VCF}
        }
        

        
        
        # ------------------- Annotations from Kalra et al. --------------------
        
        cat << 'END' > feature_extraction_thread.sh
#!/bin/bash
set -euxo pipefail

ALIGNMENTS_BAM=$1
INPUT_VCF=$2

python ~{feature_extraction_py} ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${INPUT_VCF}_features.csv 1>&2
END
        chmod +x feature_extraction_thread.sh


        # Runs a slightly modified version of the code from:
        #
        # Kalra, Paulin, Sedlazeck. "A systematic assessment of machine
        # learning for structural variant filtering." bioRxiv (2026): 2026-01.
        #
        # Remark: running the script in parallel, either over VCF chunks or per-
        # variant using xargs, does not seem to reduce runtime significantly, 
        # suggesting that some calls are the bottleneck.
        #
        # Remark: we could try to re-implement the Python script e.g. with
        # samtools. For simplicity we leave this to the future.
        #
        function FeatureExtraction() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            local PARALLELIZE=$4
            
            if [ ${PARALLELIZE} -eq 1 ]; then
                # # Chunk-based parallelization
                # bcftools view --header-only ${INPUT_VCF} > ${SAMPLE_ID}_header.txt
                # local N_RECORDS=$(bcftools view --no-header ${INPUT_VCF} | wc -l)
                # local N_RECORDS_PER_THREAD=$(( (${N_RECORDS} + ${N_THREADS} - 1) / ${N_THREADS} ))
                # bcftools view --no-header ${INPUT_VCF} | split -l ${N_RECORDS_PER_THREAD} - ${SAMPLE_ID}_chunk_
                # for FILE in ${SAMPLE_ID}_chunk_*; do
                #     cat ${SAMPLE_ID}_header.txt ${FILE} > ${FILE}.vcf
                #     rm -f ${FILE}
                #     ${TIME_COMMAND} python ~{feature_extraction_py} ${FILE}.vcf ${ALIGNMENTS_BAM} ~{reference_fa} ${FILE}_features.csv 1>&2 &
                # done
                # wait
                # rm -f ${SAMPLE_ID}_features.csv
                # for FILE in $( ls ${SAMPLE_ID}_chunk_*_features.csv | sort -V ); do
                #     cat ${FILE} >> ${SAMPLE_ID}_features.csv
                #     rm -f ${FILE}
                # done

                # Call-based parallelization
                bcftools view --header-only ${INPUT_VCF} > ${SAMPLE_ID}_header.txt
                bcftools view --no-header ${INPUT_VCF} | split -l 1 - ${SAMPLE_ID}_chunk_
                for FILE in ${SAMPLE_ID}_chunk_* ; do
                    cat ${SAMPLE_ID}_header.txt ${FILE} > ${FILE}_vcf
                    rm -f ${FILE}
                done
                ls -laht 1>&2
                rm -f ${SAMPLE_ID}_header.txt
                ls ${SAMPLE_ID}_chunk_*_vcf > list.txt
                ${TIME_COMMAND} xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} ./feature_extraction_thread.sh ${ALIGNMENTS_BAM}
                ls -laht 1>&2
                rm -f ${SAMPLE_ID}_feature_extraction_features.csv
                for FILE in $( ls ${SAMPLE_ID}_chunk_*_features.csv | sort -V ); do
                    cat ${FILE} >> ${SAMPLE_ID}_feature_extraction_features.csv
                    rm -f ${FILE}
                done
            else
                ${TIME_COMMAND} python ~{feature_extraction_py} ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${SAMPLE_ID}_feature_extraction_features.csv 1>&2
            fi
            head -n 10 ${SAMPLE_ID}_feature_extraction_features.csv 1>&2 || echo "0"
            tail -n +2 ${SAMPLE_ID}_feature_extraction_features.csv | tr ',' '\t' | cut -f 6,8-19 | sort -k 1,1 > ${SAMPLE_ID}_feature_extraction_features.tsv
            rm -f ${SAMPLE_ID}_feature_extraction_features.csv
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_feature_extraction_features.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_feature_extraction_features.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz
        }


        # A separate function just to enable running `feature_extraction_py` in 
        # parallel with other tools, if needed.
        #
        function FeatureExtraction_Annotate() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2

            echo '##INFO=<ID=FEX_DEPTH_RATIO,Number=1,Type=Float,Description="depth_ratio from feature_extraction">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_DEPTH_MAD,Number=1,Type=Float,Description="depth_mad from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_AB,Number=1,Type=Float,Description="ab from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_CN_SLOP,Number=1,Type=Float,Description="cn_slop from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_MQ_DROP,Number=1,Type=Float,Description="mq_drop from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_CLIP_FRAC,Number=1,Type=Float,Description="clip_frac from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_SPLIT_READS,Number=1,Type=Float,Description="split_reads from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_READ_LEN_MED,Number=1,Type=Float,Description="read_len_med from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_STRAND_BIAS,Number=1,Type=Float,Description="strand_bias from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_GC_FRAC,Number=1,Type=Float,Description="gc_frac from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_HOMOPOLYMER_MAX,Number=1,Type=Float,Description="homopolymer_max from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_LCR_MASK,Number=1,Type=Float,Description="lcr_mask from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/FEX_DEPTH_RATIO,INFO/FEX_DEPTH_MAD,INFO/FEX_AB,INFO/FEX_CN_SLOP,INFO/FEX_MQ_DROP,INFO/FEX_CLIP_FRAC,INFO/FEX_SPLIT_READS,INFO/FEX_READ_LEN_MED,INFO/FEX_STRAND_BIAS,INFO/FEX_GC_FRAC,INFO/FEX_HOMOPOLYMER_MAX,INFO/FEX_LCR_MASK'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        
        
        # -------------------- Annotations from genotypers ---------------------
        
        # Remark: extracting the BAM with a 1kbp slack around each call does not
        # make cuteFC output exactly the same annotations as with the original
        # BAM. Maybe such annotations are still useful for filtering, but we
        # skip this analysis for simplicity.
        #
        # Remark: cuteFC does not seem to run faster when more threads than CPU
        # cores are given to it (hyperthreading). In general, cuteFC seems to
        # make poor use of multiple threads.
        #
        function Cutefc() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            local CUTEFC_N_THREADS=$4
            local EXTRACT_BAM=$5

            local SLACK_BP=1000  # Arbitrary
            
            mkdir ./cutefc_dir/
            ${TIME_COMMAND} cuteFC --threads ${CUTEFC_N_THREADS} --genotype --max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -Ivcf ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${SAMPLE_ID}_cutefc.vcf ./cutefc_dir
            rm -rf ./cutefc_dir
            bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t[%PL]\t%INFO/CIPOS\t%INFO/CILEN\t%INFO/RE\t%INFO/STRAND\n' ${SAMPLE_ID}_cutefc.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($6=="0/0" || $6=="0|0" || $6=="./."  || $6==".|." || $6=="./0" || $6==".|0" || $6=="0/." || $6=="0|." || $6=="0" || $6==".") GT_COUNT=0; \
                else if ($6=="0/1" || $6=="0|1" || $6=="1/0" || $6=="1|0" || $6=="./1" || $6==".|1" || $6=="1/." || $6=="1|." || $6=="1") GT_COUNT=1; \
                else if ($6=="1/1" || $6=="1|1") GT_COUNT=2; \
                \
                if ($7==".") GQ=-1; \
                else GQ=$7; \
                \
                if ($8==".") DR=-1; \
                else DR=$8; \
                \
                if ($9==".") DV=-1; \
                else DV=$9; \
                \
                PL_1=-1; PL_2=-1; PL_3=-1; \
                p=0; \
                for (i=1; i<=length($10); i++) { \
                    if (substr($10,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    PL_1=substr($10,1,p-1); \
                    q=0; \
                    for (i=p+1; i<=length($10); i++) { \
                        if (substr($10,i,1)==",") { q=i; break; } \
                    } \
                    if (q>0) { \
                        PL_2=substr($10,p+1,q-1-p); \
                        PL_3=substr($10,q+1); \
                    } \
                    else { PL_2=substr($10,p+1); } \
                } \
                else { PL_1=$10; }
                \
                CIPOS_1=-1; CIPOS_2=-1; \
                p=0; \
                for (i=1; i<=length($11); i++) { \
                    if (substr($11,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CIPOS_1=substr($11,1,p-1); \
                    CIPOS_2=substr($11,p+1); \
                } \
                else { CIPOS_1=$11; } \
                \
                CILEN_1=-1; CILEN_2=-1; \
                p=0; \
                for (i=1; i<=length($12); i++) { \
                    if (substr($12,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CILEN_1=substr($12,1,p-1); \
                    CILEN_2=substr($12,p+1); \
                } \
                else { CILEN_1=$12; } \
                \
                if ($13==".") RE=-1; \
                else RE=$13; \
                \
                STRAND=-1; \
                if ($14=="--") STRAND=0; \
                else if ($14=="-+") STRAND=1; \
                else if ($14=="+-") STRAND=2; \
                else if ($14=="++") STRAND=3; \
                \
                printf("%s\t%d\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,$4,$5,  GT_COUNT,GQ,DR/RE,DV/RE,PL_1,PL_2,PL_3,CIPOS_1,CIPOS_2,CILEN_1,CILEN_2,STRAND); \
            }' | bgzip -c > ${SAMPLE_ID}_cutefc_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_cutefc.vcf
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_cutefc_annotations.tsv.gz
            if [ ${EXTRACT_BAM} -eq 1 ]; then
                rm -f ${SAMPLE_ID}_extracted.bam*
            fi
        }


        # A separate function just to enable running cuteFC in parallel with 
        # other tools, if needed.
        #
        function Cutefc_Annotate() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
        
            echo '##INFO=<ID=CUTEFC_GT_COUNT,Number=1,Type=Integer,Description="Cutefc GT converted to an integer in {0,1,2}.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_GQ,Number=1,Type=Integer,Description="Genotype quality according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_DR,Number=1,Type=Float,Description="High-quality reference reads according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_DV,Number=1,Type=Float,Description="High-quality variant reads according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_1,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_2,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_3,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_1,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_2,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_1,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_2,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_STRAND,Number=1,Type=Integer,Description="Cutefc strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--) converted to an integer in {0,1,2,3}.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/CUTEFC_GT_COUNT,INFO/CUTEFC_GQ,INFO/CUTEFC_DR,INFO/CUTEFC_DV,INFO/CUTEFC_PL_1,INFO/CUTEFC_PL_2,INFO/CUTEFC_PL_3,INFO/CUTEFC_CIPOS_1,INFO/CUTEFC_CIPOS_2,INFO/CUTEFC_CILEN_1,INFO/CUTEFC_CILEN_2,INFO/CUTEFC_STRAND'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_cutefc_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_cutefc_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        



        # --------------------- Repeat track annotations -----------------------

        function VcfToBed_StartEnd() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2

            ${TIME_COMMAND} java -cp ~{docker_dir} BndGetBins ${INPUT_VCF} ~{reference_fai} 0 0 | sort -k1,1 -k2,2n > ${SAMPLE_ID}_bins.bed
            grep -P '\t0$' ${SAMPLE_ID}_bins.bed > ${SAMPLE_ID}_start.bed
            grep -P '\t2$' ${SAMPLE_ID}_bins.bed > ${SAMPLE_ID}_end.bed
            rm -f ${SAMPLE_ID}_bins.bed
        }


        function AnnotateTrack_Bnd() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local START_BED=$3
            local END_BED=$4
            local TRACK_BED=$5
            local TRACK_ID=$6

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${START_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t1\n",$4); }' > ${SAMPLE_ID}_start_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -a ${START_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t0\n",$4); }' >> ${SAMPLE_ID}_start_track.tsv
            sort -k 1,1 ${SAMPLE_ID}_start_track.tsv > ${SAMPLE_ID}_start_track_sorted.tsv ; rm -f ${SAMPLE_ID}_start_track.tsv

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${END_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t1\n",$4); }' > ${SAMPLE_ID}_end_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -a ${END_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t0\n",$4); }' >> ${SAMPLE_ID}_end_track.tsv
            sort -k 1,1 ${SAMPLE_ID}_end_track.tsv > ${SAMPLE_ID}_end_track_sorted.tsv ; rm -f ${SAMPLE_ID}_end_track.tsv

            ${TIME_COMMAND} join -t $'\t' -1 1 -2 1 ${SAMPLE_ID}_start_track_sorted.tsv ${SAMPLE_ID}_end_track_sorted.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\n",$1,  $2,$3); }' | sort -k 1,1 > ${SAMPLE_ID}_track.tsv
            rm -f ${SAMPLE_ID}_start_track_sorted.tsv ${SAMPLE_ID}_end_track_sorted.tsv
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_track.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%d\t%d\n",$2,$3,$4,$5,$1,  $6,$7); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_track.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_track.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_track.tsv.gz
            echo '##INFO=<ID=START_'${TRACK_ID}',Number=1,Type=Integer,Description="START breakpoint is contained in a '${TRACK_ID}'">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=END_'${TRACK_ID}',Number=1,Type=Integer,Description="END breakpoint is contained in a '${TRACK_ID}'">' >> ${SAMPLE_ID}_header.txt
            COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/START_'${TRACK_ID}',INFO/END_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_track.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_track* ${SAMPLE_ID}_header.txt
        }

        
        

        # ---------------------------- Main program ----------------------------
        
        DEFAULT_QUAL="60"   # Arbitrary
        samtools --version 1>&2
        bcftools --version 1>&2
        cuteFC --version 1>&2
        df -h 1>&2
        
        cat ~{chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            # Skipping the sample if it has already been processed
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ "${TEST}" != "0" ]; then
                continue
            fi
            MEAN_COVERAGE=$(echo ${LINE} | cut -d , -f 2)

            # 1. Canonizing the VCF
            LocalizeSample ${SAMPLE_ID} ${LINE}
            df -h 1>&2
            N_RECORDS=$( bcftools index --nrecords ${SAMPLE_ID}.vcf.gz.tbi )
            if [ ${N_RECORDS} -eq 0 ]; then
                mv ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_bnd.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_bnd.vcf.gz
                gcloud storage mv ${SAMPLE_ID}_bnd.vcf.'gz*' ~{remote_outdir}/
                touch ${SAMPLE_ID}.done
                gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
                DelocalizeSample ${SAMPLE_ID}
                continue
            fi
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            
            # 2. Custom annotations
            AnnotateCustom ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${MEAN_COVERAGE} ~{min_mapq} 0

            # 3. Repeat tracks
            VcfToBed_StartEnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf
            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_start.bed ${SAMPLE_ID}_end.bed ~{tr_bed} "TR"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_start.bed ${SAMPLE_ID}_end.bed ~{segdup_bed} "SD"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_start.bed ${SAMPLE_ID}_end.bed ~{gc_content_bed} "GC"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            rm -f ${SAMPLE_ID}_start.bed ${SAMPLE_ID}_end.bed

            # 4. Original GT from the callers
            AnnotateGtCount ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf

            # 5. Adding annotations from cuteFC and Kalra et al.
            # Remark: we could run these two annotations in parallel, since they
            # are both slow and use threads inefficiently. In practice this does
            # not decrease total runtime.
            FeatureExtraction ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}.bam 0
            FeatureExtraction_Annotate ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            if [ ~{use_cutefc} -eq 1 ]; then
                Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}.bam $(( ${N_THREADS} / 2 )) 0
                Cutefc_Annotate ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf
                rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            fi

            # Uploading
            bcftools view --output-type z ${SAMPLE_ID}_canonized.vcf --output ${SAMPLE_ID}_bnd.vcf.gz
            bcftools index -f -t ${SAMPLE_ID}_bnd.vcf.gz
            gcloud storage mv ${SAMPLE_ID}_bnd.vcf.'gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht 1>&2
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
    }
}
