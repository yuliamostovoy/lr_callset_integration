version 1.0


# Annotates BND VCFs with multiple BAM- and genotyper-derived features.
#
# Remark: for BNDs, `feature_extraction.py` outputs empty values for features:
#
# FEX_DEPTH_MAD, FEX_READ_LEN_MED, FEX_STRAND_BIAS
#
# These (along with SVLEN) should not be used for scoring.
#
workflow SV_Integration_BndAnnotate {
    input {
        File chunk_tsv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai
        
        Int custom_breakpoint_window_bp = 500
        Int custom_min_clip_length = 200
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
        chunk_tsv: "Format: `ID,?,bai,bam,?,...,?` where `?` means a single string."
        tr_bed: "From: https://github.com/PacificBiosciences/pbsv/tree/master/annotations"
        segdup_bed: "From: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/"
        gc_content_bed: "From: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/"
    }
    
    call Impl {
        input:
            chunk_tsv = chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,

            custom_breakpoint_window_bp = custom_breakpoint_window_bp,
            custom_min_clip_length = custom_min_clip_length,
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
# java UltralongBndGetBins                           200%     50M       1s
# java UltralongBndCreateBedcovAnnotations           200%     50M       1s
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
        
        Int custom_breakpoint_window_bp
        Int custom_min_clip_length
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
            
            local ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            local ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            
            ${TIME_COMMAND} gcloud storage cp ${ALIGNED_BAM} ./${SAMPLE_ID}.bam
            gcloud storage cp ${ALIGNED_BAI} ./${SAMPLE_ID}.bam.bai
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_bnd.bcf ./${SAMPLE_ID}.bcf
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_bnd.bcf.csi ./${SAMPLE_ID}.bcf.csi
            
            # Converting to .vcf.gz for downstream tools
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}.vcf.gz
            rm -f ${SAMPLE_ID}.bcf*
        }
        
        
        # Deletes all files related to the given sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}*.bam* ${SAMPLE_ID}*.bcf* ${SAMPLE_ID}*.vcf* ${SAMPLE_ID}*.tsv* ${SAMPLE_ID}*.csv*
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
        # annotates each record with the coverage of a window around POS and ALT
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

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBndGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} samtools bedcov ${SAMPLE_ID}_bins.bed ${INPUT_BAM} > ${SAMPLE_ID}_counts.bed
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBndCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed ${BREAKPOINT_WINDOW_BP} | sort -k 1,1 -k 2,2n > ${SAMPLE_ID}_tags.tsv
            rm -f ${SAMPLE_ID}_counts.bed
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\n' ${INPUT_VCF} | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv | sort -k 1,1 -k 4,4 | paste - - | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$5,$10); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POS_0,Number=1,Type=Float,Description="Coverage of the bin around the breakpoint at POS">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POS_1,Number=1,Type=Float,Description="Coverage of the bin around the breakpoint at ALT">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/BIN_POS_0,INFO/BIN_POS_1'
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

samtools view --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} | awk '{ sum+=$5; count++ } END { print (count>0?sum/count:0) }' > ${RECORD_ID}_${BIN_ID}_mapq.txt
samtools view --count ${INPUT_BAM} --require-flags 256 ${CHROM}:${START}-${END} > ${RECORD_ID}_${BIN_ID}_secondary.txt
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
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBndGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            rm -f *_mapq.txt *_secondary.txt
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_counts.tsv
            while read -u 6 ID; do
                local MAPQ_POINT_0=$(cat ${ID}_0_mapq.txt)
                local MAPQ_POINT_1=$(cat ${ID}_1_mapq.txt)
                local SECONDARY_POINT_0=$(cat ${ID}_0_secondary.txt)
                local SECONDARY_POINT_1=$(cat ${ID}_1_secondary.txt)
                echo -e "${ID}\t${MAPQ_POINT_0}\t${SECONDARY_POINT_0}\t${MAPQ_POINT_1}\t${SECONDARY_POINT_1}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 6< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$4,$5,$6,$7); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POINT_MAPQ_0,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_MAPQ_1,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY_0,Number=1,Type=Integer,Description="Breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY_1,Number=1,Type=Integer,Description="Breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/BIN_POINT_MAPQ_0,INFO/BIN_POINT_SECONDARY_0,INFO/BIN_POINT_MAPQ_1,INFO/BIN_POINT_SECONDARY_1'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
     

        cat << 'END' > annotate_clipped_alignments_1.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
CLASSPATH=$2
MIN_CLIP_LENGTH=$3
CHROM=$4
START=$5
END=$6
RECORD_ID=$7
BIN_ID=$8

samtools view --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} > ${RECORD_ID}_${BIN_ID}.sam
java -cp ${CLASSPATH} UltralongIntervalGetClips ${RECORD_ID}_${BIN_ID}.sam ${START} ${END} ${MIN_CLIP_LENGTH} ${RECORD_ID}_${BIN_ID}
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
RECORD_ID=$3

LL=$(wc -l < ${RECORD_ID}_0_leftmaximal_sorted.txt)
LR=$(wc -l < ${RECORD_ID}_0_rightmaximal_sorted.txt)
RL=$(wc -l < ${RECORD_ID}_1_leftmaximal_sorted.txt)
RR=$(wc -l < ${RECORD_ID}_1_rightmaximal_sorted.txt)
LL_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_0_leftmaximal_sorted.txt ${LL} 1 ${RECORD_ID}_1_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
LL_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_0_leftmaximal_sorted.txt ${LL} 1 ${RECORD_ID}_1_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
LR_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_0_rightmaximal_sorted.txt ${LR} 0 ${RECORD_ID}_1_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
LR_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${RECORD_ID}_0_rightmaximal_sorted.txt ${LR} 0 ${RECORD_ID}_1_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
echo -e "${RECORD_ID}\t${LL}\t${LR}\t${RL}\t${RR}\t${LL_RL}\t${LL_RR}\t${LR_RL}\t${LR_RR}" > ${RECORD_ID}_counts.txt
rm -f ${RECORD_ID}_*maximal_sorted.txt
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
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBndGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ${MIN_CLIP_LENGTH}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
            rm -f *_counts.txt
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2.sh ~{docker_dir} ${ADJACENCY_SLACK_BP}
            rm -f ${SAMPLE_ID}_variantID.txt
            cat *_counts.txt | sort -k 1,1 > ${SAMPLE_ID}_counts.tsv
            rm -f *_counts.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=LL,Number=1,Type=Integer,Description="Left breakpoint: number of left-clipped alignments.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR,Number=1,Type=Integer,Description="Left breakpoint: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=RL,Number=1,Type=Integer,Description="Right breakpoint: number of left-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=RR,Number=1,Type=Integer,Description="Right breakpoint: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_1,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_2,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_3,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_4,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_1,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_2,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_3,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_4,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_1,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_2,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_3,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_4,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a left-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_1,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_2,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_3,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_4,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left breakpoint and a right-clipped alignment in the right breakpoint.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/LL,INFO/LR,INFO/RL,INFO/RR,INFO/LL_RL_1,INFO/LL_RL_2,INFO/LL_RL_3,INFO/LL_RL_4,INFO/LL_RR_1,INFO/LL_RR_2,INFO/LL_RR_3,INFO/LL_RR_4,INFO/LR_RL_1,INFO/LR_RL_2,INFO/LR_RL_3,INFO/LR_RL_4,INFO/LR_RR_1,INFO/LR_RR_2,INFO/LR_RR_3,INFO/LR_RR_4'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }

        
        function AnnotateCustom() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            
            mv ${INPUT_VCF} ${SAMPLE_ID}_in.vcf
            
            AnnotateCoverageBins ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateMapqSecondary ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateClippedAlignments ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ~{custom_adjacency_slack_bp} ~{custom_min_clip_length}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf

            mv ${SAMPLE_ID}_in.vcf ${INPUT_VCF}
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
            tail -n +2 ${SAMPLE_ID}_feature_extraction_features.csv | tr ',' '\t' | cut -f 1,2,6,8-19 | bgzip -c > ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_feature_extraction_features.csv
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
            echo '##INFO=<ID=FEX_SPLIT_READS,Number=1,Type=Integer,Description="split_reads from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_READ_LEN_MED,Number=1,Type=Float,Description="read_len_med from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_STRAND_BIAS,Number=1,Type=Float,Description="strand_bias from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_GC_FRAC,Number=1,Type=Float,Description="gc_frac from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_HOMOPOLYMER_MAX,Number=1,Type=Integer,Description="homopolymer_max from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_LCR_MASK,Number=1,Type=Integer,Description="lcr_mask from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/FEX_DEPTH_RATIO,INFO/FEX_DEPTH_MAD,INFO/FEX_AB,INFO/FEX_CN_SLOP,INFO/FEX_MQ_DROP,INFO/FEX_CLIP_FRAC,INFO/FEX_SPLIT_READS,INFO/FEX_READ_LEN_MED,INFO/FEX_STRAND_BIAS,INFO/FEX_GC_FRAC,INFO/FEX_HOMOPOLYMER_MAX,INFO/FEX_LCR_MASK'
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
            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t[%PL]\t%INFO/CIPOS\t%INFO/CILEN\t%INFO/RE\t%INFO/STRAND\n' ${SAMPLE_ID}_cutefc.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                \
                if ($5==".") GQ=-1; \
                else GQ=$5; \
                \
                if ($6==".") DR=-1; \
                else DR=$6; \
                \
                if ($7==".") DV=-1; \
                else DV=$7; \
                \
                PL_1=-1; PL_2=-1; PL_3=-1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    PL_1=substr($8,1,p-1); \
                    q=0; \
                    for (i=p+1; i<=length($8); i++) { \
                        if (substr($8,i,1)==",") { q=i; break; } \
                    } \
                    if (q>0) { \
                        PL_2=substr($8,p+1,q-1-p); \
                        PL_3=substr($8,q+1); \
                    } \
                    else { PL_2=substr($8,p+1); } \
                } \
                else { PL_1=$8; }
                \
                CIPOS_1=-1; CIPOS_2=-1; \
                p=0; \
                for (i=1; i<=length($9); i++) { \
                    if (substr($9,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CIPOS_1=substr($9,1,p-1); \
                    CIPOS_2=substr($9,p+1); \
                } \
                else { CIPOS_1=$9; } \
                \
                CILEN_1=-1; CILEN_2=-1; \
                p=0; \
                for (i=1; i<=length($10); i++) { \
                    if (substr($10,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CILEN_1=substr($10,1,p-1); \
                    CILEN_2=substr($10,p+1); \
                } \
                else { CILEN_1=$10; } \
                \
                if ($11==".") RE=-1; \
                else RE=$11; \
                \
                STRAND=-1; \
                if ($12=="--") STRAND=0; \
                else if ($12=="-+") STRAND=1; \
                else if ($12=="+-") STRAND=2; \
                else if ($12=="++") STRAND=3; \
                \
                printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,GT_COUNT,GQ,DR,DV,PL_1,PL_2,PL_3,CIPOS_1,CIPOS_2,CILEN_1,CILEN_2,RE,STRAND); \
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
            echo '##INFO=<ID=CUTEFC_DR,Number=1,Type=Integer,Description="High-quality reference reads according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_DV,Number=1,Type=Integer,Description="High-quality variant reads according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_1,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_2,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_3,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_1,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_2,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_1,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_2,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_RE,Number=1,Type=Integer,Description="Number of read support this record according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_STRAND,Number=1,Type=Integer,Description="Cutefc strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--) converted to an integer in {0,1,2,3}.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/CUTEFC_GT_COUNT,INFO/CUTEFC_GQ,INFO/CUTEFC_DR,INFO/CUTEFC_DV,INFO/CUTEFC_PL_1,INFO/CUTEFC_PL_2,INFO/CUTEFC_PL_3,INFO/CUTEFC_CIPOS_1,INFO/CUTEFC_CIPOS_2,INFO/CUTEFC_CILEN_1,INFO/CUTEFC_CILEN_2,INFO/CUTEFC_RE,INFO/CUTEFC_STRAND'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_cutefc_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_cutefc_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        



        # --------------------- Repeat track annotations -----------------------

        function VcfToBed_StartEnd() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBndGetBins ${INPUT_VCF} ~{reference_fai} 0 | sort -k1,1 -k2,2n > ${SAMPLE_ID}_bins.bed
            grep -P '\t0$' ${SAMPLE_ID}_bins.bed > ${SAMPLE_ID}_start.bed
            grep -P '\t1$' ${SAMPLE_ID}_bins.bed > ${SAMPLE_ID}_end.bed
            rm -f ${SAMPLE_ID}_bins.bed
        }


        function AnnotateTrack_Bnd() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local POINT_BED=$3
            local POINT_ID=$4
            local TRACK_BED=$5
            local TRACK_ID=$6

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${POINT_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t1\n",$1,$2,$4); }' > ${SAMPLE_ID}_${POINT_ID}_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -a ${POINT_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t0\n",$1,$2,$4); }' >> ${SAMPLE_ID}_${POINT_ID}_track.tsv
            sort -k 1,1 -k 2,2n ${SAMPLE_ID}_${POINT_ID}_track.tsv | bgzip > ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz
            tabix -@ ${N_THREADS} -0 -f -s1 -b2 -e2 ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz
            echo '##INFO=<ID='${POINT_ID}'_'${TRACK_ID}',Number=1,Type=Integer,Description="'${POINT_ID}' breakpoint is contained in a '${TRACK_ID}'">' > ${SAMPLE_ID}_header.txt
            COLUMNS='CHROM,POS,~ID,INFO/'${POINT_ID}'_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_${POINT_ID}_track* ${SAMPLE_ID}_header.txt
        }

        
        
                
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"  # Arbitrary
        DEFAULT_QUAL="60"   # Arbitrary
        samtools --version 1>&2
        bcftools --version 1>&2
        cuteFC --version 1>&2
        df -h 1>&2
        
        cat ~{chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            # Skipping the sample if it has already been processed
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi

            # 1. Canonizing the VCF
            LocalizeSample ${SAMPLE_ID} ${LINE}
            df -h 1>&2
            N_RECORDS=$(bcftools view --no-header ${SAMPLE_ID}.vcf.gz | wc -l)
            if [ ${N_RECORDS} -eq 0 ]; then
                continue
            fi
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            
            # 2. Custom annotations
            AnnotateCustom ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf

            # 3. Repeat tracks
            VcfToBed_StartEnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf

            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_start.bed "START" ~{tr_bed} "TR"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_end.bed "END" ~{tr_bed} "TR"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf

            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_start.bed "START" ~{segdup_bed} "SD"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_end.bed "END" ~{segdup_bed} "SD"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf

            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_start.bed "START" ~{gc_content_bed} "GC"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf
            AnnotateTrack_Bnd ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf ${SAMPLE_ID}_end.bed "END" ~{gc_content_bed} "GC"
            rm -f ${SAMPLE_ID}_canonized.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_canonized.vcf

            # 4. Adding annotations from cuteFC and Kalra et al.
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
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
