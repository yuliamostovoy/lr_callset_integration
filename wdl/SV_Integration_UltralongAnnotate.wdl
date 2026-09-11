version 1.0


# Annotates ultralong VCFs with multiple BAM- and genotyper-derived features,
# and splits them by SVTYPE.
#
# Remark: the entire ultralong pipeline supports symbolic <INS>, which are never
# removed.
#
workflow SV_Integration_UltralongAnnotate {
    input {
        File chunk_tsv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai
        
        String min_mapq = "0,60"

        Int ins2dup_bin_length = 100
        Float ins2dup_bin_coverage_ratio = 1.5
        Int convert_ins_to_dup = 1
        
        Int custom_breakpoint_window_bp = 500
        Int custom_min_clip_length = 200
        String custom_min_indel_length = "500,1000,5000"
        Int custom_adjacency_slack_bp = 300

        Int use_cutefc = 0

        File tr_bed
        File segdup_bed
        File gc_content_bed
        Float repeat_overlap_fraction = 0.8
        File mei_bed_gz
        File mei_bed_tbi
        Int mei_bed_tolerance = 300
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 3
    }
    parameter_meta {
        chunk_tsv: "Format: `ID,mean_coverage,bai,bam`."
        convert_ins_to_dup: "1=INS records that correspond to duplications are rewritten as DUP records"
        tr_bed: "From: https://github.com/PacificBiosciences/pbsv/tree/master/annotations"
        segdup_bed: "From: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/"
        gc_content_bed: "From: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/"
        mei_bed_gz: "A sorted, indexed, 3-column BED. Can be built with: gunzip -c repeatmasker.tsv.gz | awk 'BEGIN { FS = \"\t\"; OFS = \"\t\"; } $12==\"LINE\" || $12==\"SINE\" || $12==\"LTR\" || $12==\"DNA\" { print $6, $7, $8 }' | sort -k1,1 -k2,2n mei.bed | bzgip > mei.bed.gz ; tabix -p bed mei.bed.gz"
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

            ins2dup_bin_length = ins2dup_bin_length,
            ins2dup_bin_coverage_ratio = ins2dup_bin_coverage_ratio,
            convert_ins_to_dup = convert_ins_to_dup,

            custom_breakpoint_window_bp = custom_breakpoint_window_bp,
            custom_min_clip_length = custom_min_clip_length,
            custom_min_indel_length = custom_min_indel_length,
            custom_adjacency_slack_bp = custom_adjacency_slack_bp,
    
            use_cutefc = use_cutefc,

            tr_bed = tr_bed,
            segdup_bed = segdup_bed,
            gc_content_bed = gc_content_bed,
            repeat_overlap_fraction = repeat_overlap_fraction,
            mei_bed_gz = mei_bed_gz,
            mei_bed_tbi = mei_bed_tbi,
            mei_bed_tolerance = mei_bed_tolerance,

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
# samtools bedcov (1 thread)                          80%    900M      15m
# samtools bedcov (4 threads)                        350%    900M       5m
# annotate_mapq_secondary.sh                         400%     15M      10s
# bcftools annotate                                  100%     15M      50s
# java UltralongIntervalGetBins                      200%     50M       1s
# java UltralongIntervalCreateBedcovAnnotations      200%     50M       1s
# annotate_clipped_alignments_1.sh                   300%    500M       2m
# annotate_clipped_alignments_2.sh                   300%     50M       2m
#
# cutefc (1 thread)                                   30%    1.5G      50m
# cutefc (2 threads)                                  50%    1.5G      25m
# cutefc (4 threads)                                  30%    1.5G      50m
# cutefc (2 threads, extracted BAM)                  200%    900M       2m
# samtools view (for extracted BAM, 2 threads)       100%     20M      10m
#
# feature_extraction.py (1 thread)                   100%     12G       7m
# feature_extraction.py (4 threads, 4 chunks)         50%   2-12G     3-9m 
# feature_extraction.py (4 threads, per-call paral.) 300%     12G      10m
#
# UltralongInsGetIntervals                           200%     50M       1m
# xargs interval_2_breakpoints.sh                    200%    1.5G       1m
# UltralongInsExtractDups                            200%     60M       1s
# truvari collapse                                   100%    100M       1m
#
task Impl {
    input {
        File chunk_tsv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai

        String min_mapq

        Int ins2dup_bin_length
        Float ins2dup_bin_coverage_ratio
        Int convert_ins_to_dup
        
        Int custom_breakpoint_window_bp
        Int custom_min_clip_length
        String custom_min_indel_length
        Int custom_adjacency_slack_bp

        Int use_cutefc

        File tr_bed
        File segdup_bed
        File gc_content_bed
        Float repeat_overlap_fraction
        File mei_bed_gz
        File mei_bed_tbi
        Int mei_bed_tolerance
        
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
            
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_ultralong.bcf ./${SAMPLE_ID}.bcf
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_ultralong.bcf.csi ./${SAMPLE_ID}.bcf.csi

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
            
            # 1. Removing SVLEN from symbolic ALTs and fixing END.
            # This is necessary for `truvari bench` to work on symbolic records:
            # https://github.com/ACEnglish/truvari/issues/290
            java -cp ~{docker_dir} FixUltralongRecords ${INPUT_VCF_GZ} ~{reference_fai} > ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF_GZ}* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 2. Cleaning REF, ALT, QUAL, FILTER.
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_in.vcf ${DEFAULT_QUAL} > ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 3. Making sure IDs will be distinct at the inter-sample level 
            # (they are already distinct at the intra-sample level, thanks to
            # the steps upstream).
            ${TIME_COMMAND} bcftools annotate --set-id ${SAMPLE_ID}'_%ID' --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_canonized.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_canonized.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
        }
        
        
        
        
        # ------------------------- Custom annotations -------------------------

cat << 'END' > samtools_bedcov_thread.sh
#!/bin/bash
set -euxo pipefail

MIN_MAPQ=$1
INPUT_BAM=$2
INPUT_BED=$3

samtools bedcov -j --min-MQ ${MIN_MAPQ} ${INPUT_BED} ${INPUT_BAM} > ${INPUT_BED}_counts.bed
END
        chmod +x samtools_bedcov_thread.sh


        # Given an input VCF containing only interval records, the procedure
        # annotates each record with coverage measures extracted from a BAM.
        #
        # Remark: `samtools bedcov` skips reads with any of the following flags
        # set: UNMAP, SECONDARY, QCFAIL, DUP
        #
        function AnnotateCoverageBins_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local N_BINS=$4
            local BREAKPOINT_WINDOW_BP=$5
            local BEDCOV_N_THREADS=$6
            local MEAN_COVERAGE=$7
            local MIN_MAPQ=$8

            # Running `samtools bedcov` in parallel, since it can be very slow
            # for large intervals.
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} ${N_BINS} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
            split -d -a 5 -l 1 ${SAMPLE_ID}_bins.bed ${SAMPLE_ID}_chunk_
            ls ${SAMPLE_ID}_chunk_* | sort -V > ${SAMPLE_ID}_list.txt
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_list.txt --max-lines=1 --max-procs=${BEDCOV_N_THREADS} ./samtools_bedcov_thread.sh ${MIN_MAPQ} ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_counts.bed
            while read -u 4 CHUNK || [ -n "${CHUNK}" ]; do
                cat ${CHUNK}_counts.bed >> ${SAMPLE_ID}_counts.bed
            done 4< ${SAMPLE_ID}_list.txt
            rm -f ${SAMPLE_ID}_chunk_* ${SAMPLE_ID}_list.txt 
            
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed $(( ${N_BINS} + 4 )) ${MEAN_COVERAGE} | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
            rm -f ${SAMPLE_ID}_counts.bed
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ${INPUT_VCF} | sort -k 5,5 > ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 5 -2 1 ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_tags.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_tags.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_BEFORE_COVERAGE_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin before the call">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_LEFT_COVERAGE_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the left-breakpoint bin">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/BIN_BEFORE_COVERAGE_'${MIN_MAPQ}',INFO/BIN_LEFT_COVERAGE_'${MIN_MAPQ}
            for i in $( seq 1 ${N_BINS} ); do
                echo '##INFO=<ID=BIN_'${i}'_COVERAGE_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the i-th bin">' >> ${SAMPLE_ID}_header.txt
                COLUMNS=${COLUMNS}',INFO/BIN_'${i}'_COVERAGE_'${MIN_MAPQ}
            done
            echo '##INFO=<ID=BIN_RIGHT_COVERAGE_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the right-breakpoint bin">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_AFTER_COVERAGE_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin after the call">' >> ${SAMPLE_ID}_header.txt
            COLUMNS=${COLUMNS}',INFO/BIN_RIGHT_COVERAGE_'${MIN_MAPQ}',INFO/BIN_AFTER_COVERAGE_'${MIN_MAPQ}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        # Given an input VCF containing only point records, the procedure
        # annotates each record with the coverage of a window around POS 
        # extracted from a BAM.
        #
        # Remark: `samtools bedcov` skips reads with any of the following flags
        # set: UNMAP, SECONDARY, QCFAIL, DUP
        #
        function AnnotateCoverageBins_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local MEAN_COVERAGE=$5
            local MIN_MAPQ=$6

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} samtools bedcov -j --min-MQ ${MIN_MAPQ} ${SAMPLE_ID}_bins.bed ${INPUT_BAM} > ${SAMPLE_ID}_counts.bed
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed ${BREAKPOINT_WINDOW_BP} ${MEAN_COVERAGE} | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
            rm -f ${SAMPLE_ID}_counts.bed
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ${INPUT_VCF} | sort -k 5,5 > ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 5 -2 1 ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_tags.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\n",$2,$3,$4,$5,$1,  $6); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_ref_alt_id.tsv ${SAMPLE_ID}_tags.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POS_'${MIN_MAPQ}',Number=1,Type=Float,Description="Coverage of the bin around the POS breakpoint">' > ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/BIN_POS_'${MIN_MAPQ}
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
ID=$5

samtools view --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} | awk '{ sum+=$5; count++ } END { printf("%.4f\n", count>0?sum/count:0) }' > ${ID}_mapq.txt
N_SECONDARY=$(samtools view --count ${INPUT_BAM} --require-flags 256 ${CHROM}:${START}-${END})
N_TOTAL=$(samtools view --count ${INPUT_BAM} ${CHROM}:${START}-${END})
if [ ${N_TOTAL} -eq 0 ]; then
    RATIO="0"
else
    RATIO=$(printf "%.4f\n" $(echo "scale=4; ${N_SECONDARY} / ${N_TOTAL}" | bc))
fi
echo "${RATIO}" > ${ID}_secondary.txt

END
        chmod +x annotate_mapq_secondary.sh


        # Given an input VCF containing only interval records, the procedure
        # annotates each record with the average MAPQ and the fraction of
        # secondary alignments (=repeat-induced multi-mappings) over its left
        # and right breakpoint, extracted from a BAM.
        #
        # Remark: we do not collect the number of supplementary alignments,
        # since more specific counts are already captured by
        # `AnnotateClippedAlignments()`.
        #
        function AnnotateMapqSecondary_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_counts.tsv
            while read -u 5 ID || [ -n "${ID}" ]; do
                local MAPQ_LEFT=$(cat ${ID}_left_mapq.txt)
                local MAPQ_RIGHT=$(cat ${ID}_right_mapq.txt)
                local SECONDARY_LEFT=$(cat ${ID}_left_secondary.txt)
                local SECONDARY_RIGHT=$(cat ${ID}_right_secondary.txt)
                echo -e "${ID}\t${MAPQ_LEFT}\t${MAPQ_RIGHT}\t${SECONDARY_LEFT}\t${SECONDARY_RIGHT}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 5< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $6,$7,$8,$9); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_LEFT_MAPQ,Number=1,Type=Float,Description="Left breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_RIGHT_MAPQ,Number=1,Type=Float,Description="Right breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_LEFT_SECONDARY,Number=1,Type=Float,Description="Left breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_RIGHT_SECONDARY,Number=1,Type=Float,Description="Right breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/BIN_LEFT_MAPQ,INFO/BIN_RIGHT_MAPQ,INFO/BIN_LEFT_SECONDARY,INFO/BIN_RIGHT_SECONDARY'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }

        
        # Given an input VCF containing only point records, the procedure
        # annotates each record with the average MAPQ and the fraction of
        # secondary alignments (=repeat-induced multi-mappings) over its POS
        # breakpoint, extracted from a BAM.
        #
        # Remark: we do not collect the number of supplementary alignments,
        # since more specific counts are already captured by
        # `AnnotateClippedAlignments()`.
        #
        function AnnotateMapqSecondary_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_counts.tsv
            while read -u 6 ID || [ -n "${ID}" ]; do
                local MAPQ_POINT=$(cat ${ID}_point_mapq.txt)
                local SECONDARY_POINT=$(cat ${ID}_point_secondary.txt)
                echo -e "${ID}\t${MAPQ_POINT}\t${SECONDARY_POINT}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 6< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $6,$7); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POINT_MAPQ,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY,Number=1,Type=Float,Description="Breakpoint window: fraction of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/BIN_POINT_MAPQ,INFO/BIN_POINT_SECONDARY'
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
ID=$9

samtools view --min-MQ ${MIN_MAPQ} --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} > ${ID}.sam
echo $(wc -l < ${ID}.sam) > ${ID}_n_alignments.txt
java -cp ${CLASSPATH} UltralongIntervalGetClips ${ID}.sam ${START} ${END} ${MIN_CLIP_LENGTH} ${MIN_INDEL_LENGTH} ${ID}
rm -f ${ID}.sam
sort -k 1,1 ${ID}_leftmaximal.txt > ${ID}_leftmaximal_sorted.txt
sort -k 1,1 ${ID}_rightmaximal.txt > ${ID}_rightmaximal_sorted.txt
rm -f ${ID}_leftmaximal.txt ${ID}_rightmaximal.txt
END
        chmod +x annotate_clipped_alignments_1.sh
    
    
        cat << 'END' > annotate_clipped_alignments_2_interval.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
ADJACENCY_SLACK_BP=$2
MEAN_COVERAGE=$3
NORMALIZATION_MODE=$4
MIN_INDEL_LENGTH=$5
ID=$6

# Counting
if [ "${NORMALIZATION_MODE}" -eq 0 ]; then
    N_ALIGNMENTS_L=${MEAN_COVERAGE}
    N_ALIGNMENTS_R=${MEAN_COVERAGE}
    NORMALIZATION_FACTOR=${MEAN_COVERAGE}
else
    N_ALIGNMENTS_L=$(cat ${ID}_left_n_alignments.txt)
    N_ALIGNMENTS_R=$(cat ${ID}_right_n_alignments.txt)
    NORMALIZATION_FACTOR=$(echo "scale=4; ( ${N_ALIGNMENTS_L} + ${N_ALIGNMENTS_R} ) / 2" | bc)
fi
rm -f ${ID}_*_n_alignments.txt
LL=$(wc -l < ${ID}_left_leftmaximal_sorted.txt)
LR=$(wc -l < ${ID}_left_rightmaximal_sorted.txt)
RL=$(wc -l < ${ID}_right_leftmaximal_sorted.txt)
RR=$(wc -l < ${ID}_right_rightmaximal_sorted.txt)
LL_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_leftmaximal_sorted.txt ${LL} 1 ${ID}_right_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR} | tr ',' '\t')
LL_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_leftmaximal_sorted.txt ${LL} 1 ${ID}_right_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR} | tr ',' '\t')
LR_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_rightmaximal_sorted.txt ${LR} 0 ${ID}_right_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR} | tr ',' '\t')
LR_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_rightmaximal_sorted.txt ${LR} 0 ${ID}_right_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} 0 ${NORMALIZATION_FACTOR} | tr ',' '\t')
rm -f ${ID}_*maximal_sorted.txt

# Normalizing
if [ $(echo "${N_ALIGNMENTS_L} == 0" | bc) -eq 1 ]; then
    LL="0"; LR="0";
else
    LL=$(printf "%.4f\n" $(echo "scale=4; ${LL} / ${N_ALIGNMENTS_L}" | bc))
    LR=$(printf "%.4f\n" $(echo "scale=4; ${LR} / ${N_ALIGNMENTS_L}" | bc))
fi
if [ $(echo "${N_ALIGNMENTS_R} == 0" | bc) -eq 1 ]; then
    RL="0"; RR="0";
else
    RL=$(printf "%.4f\n" $(echo "scale=4; ${RL} / ${N_ALIGNMENTS_R}" | bc))
    RR=$(printf "%.4f\n" $(echo "scale=4; ${RR} / ${N_ALIGNMENTS_R}" | bc))
fi
OUTPUT="${ID}\t${LL}\t${LR}\t${RL}\t${RR}\t${LL_RL}\t${LL_RR}\t${LR_RL}\t${LR_RR}"
i=0
for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
    i=$(( ${i} + 1 )); L_INS=$(cut -f ${i} ${ID}_left_indel.txt); R_INS=$(cut -f ${i} ${ID}_right_indel.txt)
    i=$(( ${i} + 1 )); L_DEL_START=$(cut -f ${i} ${ID}_left_indel.txt); R_DEL_START=$(cut -f ${i} ${ID}_right_indel.txt)
    i=$(( ${i} + 1 )); L_DEL_END=$(cut -f ${i} ${ID}_left_indel.txt); R_DEL_END=$(cut -f ${i} ${ID}_right_indel.txt)
    if [ $(echo "${N_ALIGNMENTS_L} == 0" | bc) -eq 1 ]; then
        L_INS="0"; L_DEL_START="0"; L_DEL_END="0"
    else
        L_INS=$(printf "%.4f\n" $(echo "scale=4; ${L_INS} / ${N_ALIGNMENTS_L}" | bc))
        L_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${L_DEL_START} / ${N_ALIGNMENTS_L}" | bc))
        L_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${L_DEL_END} / ${N_ALIGNMENTS_L}" | bc))
    fi
    if [ $(echo "${N_ALIGNMENTS_R} == 0" | bc) -eq 1 ]; then
        R_INS="0"; R_DEL_START="0"; R_DEL_END="0"
    else
        R_INS=$(printf "%.4f\n" $(echo "scale=4; ${R_INS} / ${N_ALIGNMENTS_R}" | bc))
        R_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${R_DEL_START} / ${N_ALIGNMENTS_R}" | bc))
        R_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${R_DEL_END} / ${N_ALIGNMENTS_R}" | bc))
    fi
    OUTPUT="${OUTPUT}\t${L_INS}\t${L_DEL_START}\t${L_DEL_END}\t${R_INS}\t${R_DEL_START}\t${R_DEL_END}"
done
rm -f ${ID}_left_indel.txt ${ID}_right_indel.txt
echo -e "${OUTPUT}" > ${ID}_counts.txt
END
        chmod +x annotate_clipped_alignments_2_interval.sh

        
        cat << 'END' > annotate_clipped_alignments_2_point.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
ADJACENCY_SLACK_BP=$2
MEAN_COVERAGE=$3
NORMALIZATION_MODE=$4
MIN_INDEL_LENGTH=$5
ID=$6

# Counting
if [ "${NORMALIZATION_MODE}" -eq 0 ]; then
    N_ALIGNMENTS=${MEAN_COVERAGE}
else
    N_ALIGNMENTS=$(cat ${ID}_point_n_alignments.txt)
fi
rm -f ${ID}_point_n_alignments.txt
PL=$(wc -l < ${ID}_point_leftmaximal_sorted.txt)
PR=$(wc -l < ${ID}_point_rightmaximal_sorted.txt)
PL_PL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_point_leftmaximal_sorted.txt ${PL} 1 ${ID}_point_leftmaximal_sorted.txt ${PL} 1 ${ADJACENCY_SLACK_BP} 1 ${N_ALIGNMENTS} | tr ',' '\t')
PL_PR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_point_leftmaximal_sorted.txt ${PL} 1 ${ID}_point_rightmaximal_sorted.txt ${PR} 0 ${ADJACENCY_SLACK_BP} 1 ${N_ALIGNMENTS} | tr ',' '\t')
PR_PR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_point_rightmaximal_sorted.txt ${PR} 0 ${ID}_point_rightmaximal_sorted.txt ${PR} 0 ${ADJACENCY_SLACK_BP} 1 ${N_ALIGNMENTS} | tr ',' '\t')
rm -f ${ID}_*maximal_sorted.txt

# Normalizing
if [ $(echo "${N_ALIGNMENTS} == 0" | bc) -eq 1 ]; then
    PL="0"; PR="0";
else
    PL=$(printf "%.4f\n" $(echo "scale=4; ${PL} / ${N_ALIGNMENTS}" | bc))
    PR=$(printf "%.4f\n" $(echo "scale=4; ${PR} / ${N_ALIGNMENTS}" | bc))
fi
OUTPUT="${ID}\t${PL}\t${PR}\t${PL_PL}\t${PL_PR}\t${PR_PR}"
i=0
for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
    i=$(( ${i} + 1 )); P_INS=$(cut -f ${i} ${ID}_point_indel.txt)
    i=$(( ${i} + 1 )); P_DEL_START=$(cut -f ${i} ${ID}_point_indel.txt)
    i=$(( ${i} + 1 )); P_DEL_END=$(cut -f ${i} ${ID}_point_indel.txt)
    if [ $(echo "${N_ALIGNMENTS} == 0" | bc) -eq 1 ]; then
        P_INS="0"; P_DEL_START="0"; P_DEL_END="0"
    else
        P_INS=$(printf "%.4f\n" $(echo "scale=4; ${P_INS} / ${N_ALIGNMENTS}" | bc))
        P_DEL_START=$(printf "%.4f\n" $(echo "scale=4; ${P_DEL_START} / ${N_ALIGNMENTS}" | bc))
        P_DEL_END=$(printf "%.4f\n" $(echo "scale=4; ${P_DEL_END} / ${N_ALIGNMENTS}" | bc))
    fi
    OUTPUT="${OUTPUT}\t${P_INS}\t${P_DEL_START}\t${P_DEL_END}"
done
rm -f ${ID}_point_indel.txt
echo -e "${OUTPUT}" > ${ID}_counts.txt
END
        chmod +x annotate_clipped_alignments_2_point.sh


        # Given an input VCF containing only interval calls, the procedure
        # annotates it with clipped alignment measures extracted from a BAM.
        #
        function AnnotateClippedAlignments_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local ADJACENCY_SLACK_BP=$5
            local MIN_CLIP_LENGTH=$6
            local MIN_INDEL_LENGTH=$7
            local MEAN_COVERAGE=$8
            local MIN_MAPQ=$9
            local NORMALIZATION_MODE=${10}

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ${MIN_CLIP_LENGTH} ${MIN_INDEL_LENGTH} ${MIN_MAPQ}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
            rm -f *_counts.txt
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2_interval.sh ~{docker_dir} ${ADJACENCY_SLACK_BP} ${MEAN_COVERAGE} ${NORMALIZATION_MODE} ${MIN_INDEL_LENGTH}
            rm -f ${SAMPLE_ID}_variantID.txt
            cat *_counts.txt | sort -k 1,1 > ${SAMPLE_ID}_counts.tsv
            rm -f *_counts.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s\t%d\t%s\t%s\t%s",$2,$3,$4,$5,$1); \
                for (i=6; i<=NF; i++) printf("\t%.4f",$i); \
                printf("\n"); \
            }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=LL_'${MIN_MAPQ}',Number=1,Type=Float,Description="Left window: number of left-clipped alignments.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_'${MIN_MAPQ}',Number=1,Type=Float,Description="Left window: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=RL_'${MIN_MAPQ}',Number=1,Type=Float,Description="Right window: number of left-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=RR_'${MIN_MAPQ}',Number=1,Type=Float,Description="Right window: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RL_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LL_RR_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RL_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=LR_RR_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
            for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
                echo '##INFO=<ID=L_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS in the left window of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=L_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start in the left window of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=L_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end in the left window of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=R_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS in the right window of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=R_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start in the right window of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=R_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end in the right window of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
            done
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/LL_'${MIN_MAPQ}',INFO/LR_'${MIN_MAPQ}',INFO/RL_'${MIN_MAPQ}',INFO/RR_'${MIN_MAPQ}',INFO/LL_RL_1_'${MIN_MAPQ}',INFO/LL_RL_2_'${MIN_MAPQ}',INFO/LL_RL_3_'${MIN_MAPQ}',INFO/LL_RL_4_'${MIN_MAPQ}',INFO/LL_RR_1_'${MIN_MAPQ}',INFO/LL_RR_2_'${MIN_MAPQ}',INFO/LL_RR_3_'${MIN_MAPQ}',INFO/LL_RR_4_'${MIN_MAPQ}',INFO/LR_RL_1_'${MIN_MAPQ}',INFO/LR_RL_2_'${MIN_MAPQ}',INFO/LR_RL_3_'${MIN_MAPQ}',INFO/LR_RL_4_'${MIN_MAPQ}',INFO/LR_RR_1_'${MIN_MAPQ}',INFO/LR_RR_2_'${MIN_MAPQ}',INFO/LR_RR_3_'${MIN_MAPQ}',INFO/LR_RR_4_'${MIN_MAPQ}
            for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
                COLUMNS=${COLUMNS}',INFO/L_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/L_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/L_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',INFO/R_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/R_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/R_DEL_END_'${MIN_MAPQ}'_'${LENGTH}
            done
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }


        # Given an input VCF containing only point records, the procedure
        # annotates it with clipped alignment measures extracted from a BAM.
        #
        function AnnotateClippedAlignments_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local ADJACENCY_SLACK_BP=$5
            local MIN_CLIP_LENGTH=$6
            local MIN_INDEL_LENGTH=$7
            local MEAN_COVERAGE=$8
            local MIN_MAPQ=$9
            local NORMALIZATION_MODE=${10}
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ${MIN_CLIP_LENGTH} ${MIN_INDEL_LENGTH} ${MIN_MAPQ}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
            rm -f *_counts.txt
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2_point.sh ~{docker_dir} ${ADJACENCY_SLACK_BP} ${MEAN_COVERAGE} ${NORMALIZATION_MODE} ${MIN_INDEL_LENGTH}
            rm -f ${SAMPLE_ID}_variantID.txt
            cat *_counts.txt | sort -k 1,1 > ${SAMPLE_ID}_counts.tsv
            rm -f *_counts.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s\t%d\t%s\t%s\t%s",$2,$3,$4,$5,$1); \
                for (i=6; i<=NF; i++) printf("\t%.4f",$i); \
                printf("\n"); \
            }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=PL_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of left-clipped alignments.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_1_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_2_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_3_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_4_'${MIN_MAPQ}',Number=1,Type=Float,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
                echo '##INFO=<ID=P_INS_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR INS of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=P_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL start of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
                echo '##INFO=<ID=P_DEL_END_'${MIN_MAPQ}'_'${LENGTH}',Number=1,Type=Float,Description="Number of reads with a CIGAR DEL end of length >='${LENGTH}'.">' >> ${SAMPLE_ID}_header.txt
            done
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/PL_'${MIN_MAPQ}',INFO/PR_'${MIN_MAPQ}',INFO/PL_PL_1_'${MIN_MAPQ}',INFO/PL_PL_2_'${MIN_MAPQ}',INFO/PL_PL_3_'${MIN_MAPQ}',INFO/PL_PL_4_'${MIN_MAPQ}',INFO/PL_PR_1_'${MIN_MAPQ}',INFO/PL_PR_2_'${MIN_MAPQ}',INFO/PL_PR_3_'${MIN_MAPQ}',INFO/PL_PR_4_'${MIN_MAPQ}',INFO/PR_PR_1_'${MIN_MAPQ}',INFO/PR_PR_2_'${MIN_MAPQ}',INFO/PR_PR_3_'${MIN_MAPQ}',INFO/PR_PR_4_'${MIN_MAPQ}
            for LENGTH in $(echo ${MIN_INDEL_LENGTH} | tr ',' ' '); do
                COLUMNS=${COLUMNS}',INFO/P_INS_'${MIN_MAPQ}'_'${LENGTH}',INFO/P_DEL_START_'${MIN_MAPQ}'_'${LENGTH}',INFO/P_DEL_END_'${MIN_MAPQ}'_'${LENGTH}
            done
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }


        function AnnotateCustom_NotIns() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local MEAN_COVERAGE=$3
            local MIN_MAPQ=$4
            local NORMALIZATION_MODE=$5

            local N_COVERAGE_BINS="10"
            local MAPQS=$( echo ${MIN_MAPQ} | tr ',' ' ' )
            
            mv ${INPUT_VCF} ${SAMPLE_ID}_in.vcf

            for MAPQ in ${MAPQS}; do
                AnnotateCoverageBins_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ${N_COVERAGE_BINS} ~{custom_breakpoint_window_bp} $(( ${N_THREADS} / 2 )) ${MEAN_COVERAGE} ${MAPQ}
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            done
            AnnotateMapqSecondary_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            for MAPQ in ${MAPQS}; do
                AnnotateClippedAlignments_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ~{custom_adjacency_slack_bp} ~{custom_min_clip_length} ~{custom_min_indel_length} ${MEAN_COVERAGE} ${MAPQ} ${NORMALIZATION_MODE}
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            done

            mv ${SAMPLE_ID}_in.vcf ${INPUT_VCF}
        }

        
        function AnnotateCustom_Ins() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local MEAN_COVERAGE=$3
            local MIN_MAPQ=$4
            local NORMALIZATION_MODE=$5

            local MAPQS=$( echo ${MIN_MAPQ} | tr ',' ' ' )
            
            mv ${INPUT_VCF} ${SAMPLE_ID}_in.vcf
            
            for MAPQ in ${MAPQS}; do
                AnnotateCoverageBins_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ${MEAN_COVERAGE} ${MAPQ}
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            done
            AnnotateMapqSecondary_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            for MAPQ in ${MAPQS}; do
                AnnotateClippedAlignments_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ~{custom_adjacency_slack_bp} ~{custom_min_clip_length} ~{custom_min_indel_length} ${MEAN_COVERAGE} ${MAPQ} ${NORMALIZATION_MODE}
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
        
        # Output fields:
        #
        # chrom,pos,end,svlen,svtype,id,log_svlen,depth_ratio,depth_mad,ab,cn_slop,mq_drop,clip_frac,split_reads,read_len_med,strand_bias,gc_frac,homopolymer_max,lcr_mask,support_read,svtype_DEL,label,dataset
        #
        cat << 'END' > feature_extraction_thread.sh
#!/bin/bash
set -euxo pipefail

ALIGNMENTS_BAM=$1
INPUT_VCF=$2

python ~{docker_dir}/feature_extraction.py ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${INPUT_VCF}_features.csv 1>&2
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
                #     ${TIME_COMMAND} python ~{docker_dir}/feature_extraction.py ${FILE}.vcf ${ALIGNMENTS_BAM} ~{reference_fa} ${FILE}_features.csv 1>&2 &
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
                ${TIME_COMMAND} python ~{docker_dir}/feature_extraction.py ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${SAMPLE_ID}_feature_extraction_features.csv 1>&2
            fi
            head -n 10 ${SAMPLE_ID}_feature_extraction_features.csv 1>&2 || echo "0"
            tail -n +2 ${SAMPLE_ID}_feature_extraction_features.csv | tr ',' '\t' | cut -f 6,8-19 | sort -k 1,1 > ${SAMPLE_ID}_feature_extraction_features.tsv
            rm -f ${SAMPLE_ID}_feature_extraction_features.csv
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_feature_extraction_features.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",$2,$3,$4,$5,$1,  $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_feature_extraction_features.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_feature_extraction_annotations.tsv.gz
        }


        # A separate function just to enable running `feature_extraction.py` in 
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

            if [ ${EXTRACT_BAM} -eq 1 ]; then
                java -cp ~{docker_dir} UltralongGetIntervals ${INPUT_VCF} ~{reference_fai} ${SLACK_BP} > ${SAMPLE_ID}_regions.bed
                ${TIME_COMMAND} samtools view --threads ${N_THREADS} --bam --regions-file ${SAMPLE_ID}_regions.bed ${ALIGNMENTS_BAM} --output ${SAMPLE_ID}_extracted.bam
                ${TIME_COMMAND} samtools index --threads ${N_THREADS} ${SAMPLE_ID}_extracted.bam
                ALIGNMENTS_BAM=${SAMPLE_ID}_extracted.bam
            fi
            
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

        function VcfToBed_Start() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
      
            bcftools query --format '%CHROM\t%POS\t%ID\n' ${INPUT_VCF} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1,$2,$2,$3); }' | sort -k1,1 -k2,2n > ${SAMPLE_ID}_start.bed
        }


        function VcfToBed_StartEndInterval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
      
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%INFO/SVLEN\t%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_matrix.tsv
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1,$2,$2,$4); }' ${SAMPLE_ID}_matrix.tsv | sort -k1,1 -k2,2n > ${SAMPLE_ID}_start.bed
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1, $2 + $3 - 1, $2 + $3 - 1, $4); }' ${SAMPLE_ID}_matrix.tsv | sort -k1,1 -k2,2n > ${SAMPLE_ID}_end.bed
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1,$2, $2 + $3 - 1, $4); }' ${SAMPLE_ID}_matrix.tsv | sort -k1,1 -k2,2n > ${SAMPLE_ID}_interval.bed
            rm -f ${SAMPLE_ID}_matrix.tsv
        }


        function AnnotateTrack_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local POINT_BED=$3
            local POINT_ID=$4
            local TRACK_BED=$5
            local TRACK_ID=$6

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${POINT_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t1\n",$4); }' > ${SAMPLE_ID}_${POINT_ID}_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -a ${POINT_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t0\n",$4); }' >> ${SAMPLE_ID}_${POINT_ID}_track.tsv
            sort -k 1,1 ${SAMPLE_ID}_${POINT_ID}_track.tsv > ${SAMPLE_ID}_${POINT_ID}_track_sorted.tsv
            rm -f ${SAMPLE_ID}_${POINT_ID}_track.tsv
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_${POINT_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_${POINT_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_${POINT_ID}_track_sorted.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%d\n",$2,$3,$4,$5,$1,  $6); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz
            rm -f ${SAMPLE_ID}_${POINT_ID}_track_sorted.tsv ${SAMPLE_ID}_${POINT_ID}_chrom_pos_id_ref_alt.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz
            echo '##INFO=<ID='${POINT_ID}'_'${TRACK_ID}',Number=1,Type=Integer,Description="'${POINT_ID}' breakpoint is contained in a '${TRACK_ID}'">' > ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/'${POINT_ID}'_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_${POINT_ID}_track* ${SAMPLE_ID}_header.txt
        }


        function AnnotateTrack_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INTERVAL_BED=$3
            local TRACK_BED=$4
            local TRACK_ID=$5
            local OVERLAP_FRACTION=$6

            ${TIME_COMMAND} bedtools intersect -wa -u -f ${OVERLAP_FRACTION} -a ${INTERVAL_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t1\n",$4); }' > ${SAMPLE_ID}_interval_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -f ${OVERLAP_FRACTION} -a ${INTERVAL_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t0\n",$4); }' >> ${SAMPLE_ID}_interval_track.tsv
            sort -k 1,1 -k 2,2n ${SAMPLE_ID}_interval_track.tsv | bgzip > ${SAMPLE_ID}_interval_track.tsv.gz
            sort -k 1,1 ${SAMPLE_ID}_interval_track.tsv > ${SAMPLE_ID}_interval_track_sorted.tsv
            rm -f ${SAMPLE_ID}_interval_track.tsv
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-5 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv ${SAMPLE_ID}_interval_track_sorted.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%s\t%s\t%d\n",$2,$3,$4,$5,$1,  $6); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_interval_track.tsv.gz
            rm -f ${SAMPLE_ID}_interval_track_sorted.tsv ${SAMPLE_ID}_chrom_pos_id_ref_alt.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_interval_track.tsv.gz
            echo '##INFO=<ID=INTERVAL_'${TRACK_ID}',Number=1,Type=Integer,Description="Interval overlaps the '${TRACK_ID}' track by at least '${OVERLAP_FRACTION}'">' > ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,REF,ALT,~ID,INFO/INTERVAL_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_interval_track.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_interval_track* ${SAMPLE_ID}_header.txt
        }


        
        
        # ---------------- Handling alternative representations ----------------
        #
        # 1. Some INS records may correspond to a DUP. Such records could be 
        # located at the beginning/end of their DUP, or in the middle of the 
        # DUP. INS records located at DUP breakpoints would see peculiar BAM 
        # features at their POS, while those located in the middle of their DUP 
        # would likely see no significant BAM features.
        #
        # 2. Assume that we keep such records in the INS VCF, without any 
        # special treatment.
        # 2.1 Assume that we then use a stringent match against both the INS and
        # the DUP truth VCF (with --dup-to-ins) to mark the true INS records. 
        # This would mark as true only INS records at the start (and perhaps at 
        # the end) of their DUP. The model would then learn to mark as true the 
        # INS records that have the BAM features of the true INS records, as 
        # well as those that have the BAM features of DUP breakpoints. This is 
        # fine, but INS records located in the middle of their DUP would not 
        # have such features and would likely be marked as false at test time:
        # this may result in losing some DUPs completely.
        # 2.2 Assume instead that we use a lenient match proportional to SVLEN
        # to mark the true INS records. This would mark as true also INS in the
        # middle of their DUP. This is fine, but since such INS have no BAM 
        # features, the model's precision may decrease.
        #
        # 3. We try to rewrite INS records as DUP using a simple heuristic that 
        # detects a longest interval with increased BAM coverage and position
        # and length compatible with the INS.
        # 3.1 Note that an INS could encode a complex DUP that corresponds to a
        # permutation of the source interval, and coverage inside the source 
        # interval may be uneven depending on how many times each source 
        # fragment occurs in the INS sequence. Our method puts in the DUP class
        # both simple and complex DUPs.
        # 3.2 We could have used `truvari anno remap`, rather than depth: this 
        # might have given more accurate breakpoints, and it might have allowed 
        # us to separate simple DUPs from complex DUPs. However, we want to 
        # perform the same classification at testing time, and we observed that
        # mapping some INS can take a large amount of RAM and thus make the task
        # expensive.
        #
        # 4. None of these problems affect the main (non-ultralong) VCF, since
        # features come from Kanpig which simply takes paths in the variation
        # graph, regardless of how they are represented in the VCF, and since 
        # the records that are marked as true by `truvari bench` are likely to 
        # capture the Kanpig features of true variants.
        #
        # 5. Alternative representations should have been handled in
        # `SV_Integration_Workpackage1.wdl` before truvari collapse. We do it
        # here just to reuse the intra-sample VCFs we already have from that
        # workflow.

        cat << 'END' > interval_2_breakpoints.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
INPUT_BAM=$2
BIN_LENGTH=$3
BIN_COVERAGE_RATIO=$4
RAM_MB=$5
MEI_BED_GZ=$6
MEI_BED_TOLERANCE=$7
MIN_DEPTH_BREAKPOINT=$8
MIN_DEPTH_INTERVAL=$9
REGION=${10}
ID=${11}

# Remark: we do not filter by MAPQ here, to get a clearer breakpoint signal.
samtools depth -aa -r ${REGION} ${INPUT_BAM} -o ${ID}_depth.tsv
tabix ${MEI_BED_GZ} ${REGION} > ${ID}_mei.bed
java -cp ${CLASSPATH} -Xmx${RAM_MB}m UltralongDepthGetBreakpoints ${ID}_depth.tsv $(wc -l < ${ID}_depth.tsv) ${BIN_LENGTH} ${BIN_COVERAGE_RATIO} ${ID}_mei.bed $(wc -l < ${ID}_mei.bed) ${MEI_BED_TOLERANCE} ${MIN_DEPTH_BREAKPOINT} ${MIN_DEPTH_INTERVAL} > ${ID}_breakpoints.tsv
rm -f ${ID}_depth.tsv ${ID}_mei.bed
END
        chmod +x interval_2_breakpoints.sh


        # Remark: this procedure creates the following files, which may be
        # empty:  `_ins.vcf`, `_ins_dup.vcf`, `_not_ins.vcf`.
        #
        # Remark: ~40% of all the ultralong INS of a sample are reclassified as
        # DUP. This is not too far from the results of `truvari anno remap`.
        #
        function Ins2Dup() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_BAM=$3
            local BIN_LENGTH=$4
            local BIN_COVERAGE_RATIO=$5
            local DEPTH_N_THREADS=$6
            local MEAN_COVERAGE=$7

            local RAM_PER_THREAD_MB=$(( ( (~{ram_size_gb} - 1) * 1024) / ${DEPTH_N_THREADS} ))
            local MIN_DEPTH_BREAKPOINT=$( echo "scale=4; ${MEAN_COVERAGE} / 4" | bc )  # Arbitrary lower bound
            local MIN_DEPTH_INTERVAL=${MEAN_COVERAGE}  # Arbitrary lower bound

            if [ ~{convert_ins_to_dup} -eq 1 ]; then
                local INSDUP_MODE=0
            else
                local INSDUP_MODE=1
            fi
            local INSDUP_QUAL=1

            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE!=\"INS\"" --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_not_ins.vcf
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE=\"INS\"" --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_ins.vcf
            local N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins.vcf | wc -l)
            if [ ${N_INS} -eq 0 ]; then
                bcftools view --header-only --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_ins_dup.vcf
                return
            fi
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongInsGetIntervals ${SAMPLE_ID}_ins.vcf ~{reference_fai} > ${SAMPLE_ID}_ins_intervals.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_ins_intervals.wsv --max-lines=1 --max-procs=${DEPTH_N_THREADS} ./interval_2_breakpoints.sh ~{docker_dir} ${INPUT_BAM} ${BIN_LENGTH} ${BIN_COVERAGE_RATIO} ${RAM_PER_THREAD_MB} ~{mei_bed_gz} ~{mei_bed_tolerance} ${MIN_DEPTH_BREAKPOINT} ${MIN_DEPTH_INTERVAL}
            rm -f ${SAMPLE_ID}_ins_intervals.wsv
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongInsExtractDups ${SAMPLE_ID}_ins.vcf . ${SAMPLE_ID}_ins_ins.vcf ${SAMPLE_ID}_ins_dup.vcf ${INSDUP_QUAL} ${INSDUP_MODE}
            rm -f *_breakpoints.tsv
            local N_INS_DUP=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins_dup.vcf | wc -l)
            if [ ${N_INS_DUP} -eq 0 ]; then
                rm -f ${SAMPLE_ID}_ins_ins.vcf
            else
                rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_ins_ins.vcf ${SAMPLE_ID}_ins.vcf
            fi
        }


        # Merges `_ins_dup.vcf` with `_not_ins.vcf` or `_ins.vcf`, removing
        # `_ins_dup.vcf` afterwards.
        #
        # Remark: `truvari collapse` is run with the same parameters as in
        # `SV_Integration_Workpackage1.wdl`. INS->DUP records are assigned
        # lower QUAL to make truvari choose an original DUP record as a 
        # cluster representative. If the DUP record has a GT it will be 
        # selected, since its column after bcftools merge is the first one.
        # This collapses ~3-5 variants per sample.
        #
        # Remark: truvari needs `bcftools merge` and it does not work with 
        # `bcftools concat`.
        function MergeInsdups() {
            local SAMPLE_ID=$1
            local N_NOT_INS=$2
            local N_INS_DUP=$3
            local N_INS=$4

            if [ ${N_INS_DUP} -eq 0 ]; then
                rm -f ${SAMPLE_ID}_ins_dup.vcf
                return
            fi
            if [ ~{convert_ins_to_dup} -eq 1 ]; then
                if [ ${N_NOT_INS} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_not_ins.vcf
                    mv ${SAMPLE_ID}_ins_dup.vcf ${SAMPLE_ID}_not_ins.vcf
                    return
                fi
                # Adding SVLEN to symbolic ALTs, to avoid overcollapse in
                # `bcftools merge`.
                ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${SAMPLE_ID}_ins_dup.vcf > ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_ins_dup.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_ins_dup.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${SAMPLE_ID}_not_ins.vcf | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_not_ins.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_ins.vcf.gz
                ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_ins_dup.vcf --output ${SAMPLE_ID}_ins_dup.vcf.gz
                rm -f ${SAMPLE_ID}_ins_dup.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins_dup.vcf.gz
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --force-samples --output-type v ${SAMPLE_ID}_not_ins.vcf.gz ${SAMPLE_ID}_ins_dup.vcf.gz --output ${SAMPLE_ID}_out.vcf
                # Removing SVLEN from symbolic ALTs
                bcftools view --header-only ${SAMPLE_ID}_out.vcf --output ${SAMPLE_ID}_not_ins.vcf
                ${TIME_COMMAND} bcftools view --no-header ${SAMPLE_ID}_out.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                    printf("%s",$1); \
                    for (i=2; i<=NF; i++) printf("\t%s",$i); \
                    printf("\n"); \
                }' >> ${SAMPLE_ID}_not_ins.vcf
                rm -f ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_not_ins.vcf.gz* ${SAMPLE_ID}_ins_dup.vcf.gz* ; bgzip --threads ${N_THREADS} --compress-level 1 ${SAMPLE_ID}_not_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_ins.vcf.gz
                ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_not_ins.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf --removed-output /dev/null
                rm -f ${SAMPLE_ID}_not_ins.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_not_ins.vcf
            else
                if [ ${N_INS} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_ins.vcf
                    mv ${SAMPLE_ID}_ins_dup.vcf ${SAMPLE_ID}_ins.vcf
                    return
                fi
                # In this case both VCFs have non-symbolic ALTs
                ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level 1 ${SAMPLE_ID}_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins.vcf.gz
                ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level 1 ${SAMPLE_ID}_ins_dup.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins_dup.vcf.gz
                ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type v ${SAMPLE_ID}_ins.vcf.gz ${SAMPLE_ID}_ins_dup.vcf.gz --output ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_ins.vcf* ${SAMPLE_ID}_ins_dup.vcf* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_ins.vcf
            fi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"  # Arbitrary
        DEFAULT_QUAL="60"   # Arbitrary
        samtools --version 1>&2
        bcftools --version 1>&2
        cuteFC --version 1>&2
        truvari --help 1>&2
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

            # 1. Canonizing the VCF and handling alternative representations
            LocalizeSample ${SAMPLE_ID} ${LINE}
            df -h 1>&2
            N_RECORDS=$( bcftools index --nrecords ${SAMPLE_ID}.vcf.gz.tbi )
            if [ ${N_RECORDS} -eq 0 ]; then
                cp ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_del.vcf.gz
                cp ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_dup.vcf.gz
                cp ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_insdup.vcf.gz
                cp ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_inv.vcf.gz
                cp ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_ins.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_del.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_dup.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_insdup.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_inv.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_ins.vcf.gz
                gcloud storage mv ${SAMPLE_ID}_'del.vcf.gz*' ${SAMPLE_ID}_'inv.vcf.gz*' ${SAMPLE_ID}_'dup.vcf.gz*' ${SAMPLE_ID}_'insdup.vcf.gz*' ${SAMPLE_ID}_'ins.vcf.gz*' ~{remote_outdir}/
                touch ${SAMPLE_ID}.done
                gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
                DelocalizeSample ${SAMPLE_ID}
                continue
            fi
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            Ins2Dup ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam ~{ins2dup_bin_length} ~{ins2dup_bin_coverage_ratio} $(( ${N_THREADS} / 2 )) ${MEAN_COVERAGE}
            rm -f ${SAMPLE_ID}_canonized.vcf.gz

            # 2. Adding custom annotations
            N_NOT_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_not_ins.vcf | wc -l)
            if [ ${N_NOT_INS} -gt 0 ]; then
                AnnotateCustom_NotIns ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${MEAN_COVERAGE} ~{min_mapq} 0
            fi
            N_INS_DUP=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins_dup.vcf | wc -l)
            if [ ${N_INS_DUP} -gt 0 ]; then
                if [ ~{convert_ins_to_dup} -eq 1 ]; then
                    AnnotateCustom_NotIns ${SAMPLE_ID} ${SAMPLE_ID}_ins_dup.vcf ${MEAN_COVERAGE} ~{min_mapq} 1
                else
                    AnnotateCustom_Ins ${SAMPLE_ID} ${SAMPLE_ID}_ins_dup.vcf ${MEAN_COVERAGE} ~{min_mapq} 1
                fi
            fi
            N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins.vcf | wc -l)
            if [ ${N_INS} -gt 0 ]; then
                AnnotateCustom_Ins ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf ${MEAN_COVERAGE} ~{min_mapq} 0
            fi
            MergeInsdups ${SAMPLE_ID} ${N_NOT_INS} ${N_INS_DUP} ${N_INS}

            # 3. Adding repeat track annotations
            # 3.1 Not INS
            N_NOT_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_not_ins.vcf | wc -l)
            if [ ${N_NOT_INS} -gt 0 ]; then
                VcfToBed_StartEndInterval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf

                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{tr_bed} "TR"
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_end.bed "END" ~{tr_bed} "TR"
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
                AnnotateTrack_Interval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_interval.bed ~{tr_bed} "TR" ~{repeat_overlap_fraction}
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf

                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{segdup_bed} "SD"
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_end.bed "END" ~{segdup_bed} "SD"
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
                AnnotateTrack_Interval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_interval.bed ~{segdup_bed} "SD" ~{repeat_overlap_fraction}
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf

                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{gc_content_bed} "GC"
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_end.bed "END" ~{gc_content_bed} "GC"
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
                AnnotateTrack_Interval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_interval.bed ~{gc_content_bed} "GC" ~{repeat_overlap_fraction}
                rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf

                rm -f ${SAMPLE_ID}_start.bed ${SAMPLE_ID}_end.bed ${SAMPLE_ID}_interval.bed
            fi
            # 3.2 INS
            N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins.vcf | wc -l)
            if [ ${N_INS} -gt 0 ]; then
                VcfToBed_Start ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf

                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{tr_bed} "TR"
                rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_ins.vcf

                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{segdup_bed} "SD"
                rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_ins.vcf

                AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{gc_content_bed} "GC"
                rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_ins.vcf

                rm -f ${SAMPLE_ID}_start.bed
            fi

            # 4. Merging INS and non-INS VCFs, so that both of them get
            # annotated in a single pass over the BAM later.
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level 1 ${SAMPLE_ID}_not_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_ins.vcf.gz
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level 1 ${SAMPLE_ID}_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins.vcf.gz
            ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type v ${SAMPLE_ID}_not_ins.vcf.gz ${SAMPLE_ID}_ins.vcf.gz --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_not_ins.vcf* ${SAMPLE_ID}_ins.vcf* ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf

            # 5. Adding the original GT from the callers
            AnnotateGtCount ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf

            # 6. Adding annotations from cuteFC and Kalra et al.
            # Remark: we could run these two annotations in parallel, since they
            # are both slow and use threads inefficiently. In practice this does
            # not decrease total runtime.
            FeatureExtraction ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam 0
            FeatureExtraction_Annotate ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            if [ ~{use_cutefc} -eq 1 ]; then
                Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam $(( ${N_THREADS} / 2 )) 0
                Cutefc_Annotate ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            fi
            
            # Ensuring that the VCF has the correct header
            bcftools view --header-only ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_header.txt
            grep -q INSDUP ${SAMPLE_ID}_header.txt && FOUND="1" || FOUND="0"
            if [ ${FOUND} = "0" ]; then
                N_ROWS=$(wc -l < ${SAMPLE_ID}_header.txt)
                head -n $(( ${N_ROWS} - 1 )) ${SAMPLE_ID}_header.txt > ${SAMPLE_ID}_header_new.txt
                echo '##INFO=<ID=INSDUP,Number=0,Type=Flag,Description="The record is the result of an INS-DUP conversion">' >> ${SAMPLE_ID}_header_new.txt
                echo '##INFO=<ID=INS_POS,Number=1,Type=Integer,Description="The POS of the original INS record">' >> ${SAMPLE_ID}_header_new.txt
                echo '##INFO=<ID=INS_ALT,Number=1,Type=String,Description="The ALT allele of the original INS record">' >> ${SAMPLE_ID}_header_new.txt
                echo '##INFO=<ID=INS_SVLEN,Number=1,Type=Integer,Description="The length of the ALT allele of the original INS record">' >> ${SAMPLE_ID}_header_new.txt
                echo '##INFO=<ID=INS_QUAL,Number=1,Type=Integer,Description="The QUAL of the original INS record">' >> ${SAMPLE_ID}_header_new.txt
                tail -n 1 ${SAMPLE_ID}_header.txt >> ${SAMPLE_ID}_header_new.txt
                bcftools reheader --header ${SAMPLE_ID}_header_new.txt ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
                rm -f ${SAMPLE_ID}_header_new.txt
            fi
            rm -f ${SAMPLE_ID}_header.txt

            # Splitting by SVTYPE and uploading
            bcftools filter --include 'SVTYPE="DEL"' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_del.vcf.gz &
            bcftools filter --include 'SVTYPE="INV"' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_inv.vcf.gz &
            if [ ~{convert_ins_to_dup} -eq 1 ]; then
                bcftools filter --include 'SVTYPE="INS"' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_ins.vcf.gz &
                bcftools filter --include 'SVTYPE="DUP" && INFO/INSDUP=0' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_dup.vcf.gz &
                bcftools filter --include 'SVTYPE="DUP" && INFO/INSDUP=1' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_insdup.vcf.gz &
            else
                bcftools filter --include 'SVTYPE="DUP"' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_dup.vcf.gz &
                bcftools filter --include 'SVTYPE="INS" && INFO/INSDUP=0' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_ins.vcf.gz &
                bcftools filter --include 'SVTYPE="INS" && INFO/INSDUP=1' --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_insdup.vcf.gz &
            fi
            wait
            bcftools index -f -t ${SAMPLE_ID}_del.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_dup.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_insdup.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_inv.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_ins.vcf.gz &
            wait
            gcloud storage mv ${SAMPLE_ID}_'del.vcf.gz*' ${SAMPLE_ID}_'inv.vcf.gz*' ${SAMPLE_ID}_'dup.vcf.gz*' ${SAMPLE_ID}_'insdup.vcf.gz*' ${SAMPLE_ID}_'ins.vcf.gz*' ~{remote_outdir}/
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
