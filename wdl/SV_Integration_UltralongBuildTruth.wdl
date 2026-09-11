version 1.0


# Given a set of per-sample, canonized svim-asm VCFs and dipcall BEDs, the 
# program prepares per-sample, per-SVTYPE truth VCFs to be later used by 
# `SV_Integration_UltralongGetTrainingIntervals.wdl` .
#
workflow SV_Integration_UltralongBuildTruth {
    input {
        File samples_tsv
        
        String remote_indir_svimasm
        String remote_outdir
        
        Int svimasm_min_sv_length = 5000
        Int svimasm_convert_ins_to_dup = 1

        Int svimasm_ins_use_gaps = 0
        Int svimasm_ins_use_gaps_slack_bp = 200
        Float svimasm_ins_use_gaps_length_similarity = 0.9

        Int svimasm_ins_use_remap = 1
        Int svimasm_ins_remap_max_length = 2000000
        Float svimasm_ins_remap_cov_threshold = 0.8

        File reference_fa
        File reference_fai
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID, DIPCALL_BED"
        remote_indir_svimasm: "Without final slash. Contains per-sample canonized and filtered svim-asm VCFs."
        svimasm_min_sv_length: "Removes SVs shorter than this from the svim-asm VCF to speed up all steps."
        svimasm_convert_ins_to_dup: "1=SVIM-asm's INS records that correspond to duplications are rewritten as DUP records"
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,

            remote_indir_svimasm = remote_indir_svimasm,
            remote_outdir = remote_outdir,

            svimasm_min_sv_length = svimasm_min_sv_length,
            svimasm_convert_ins_to_dup = svimasm_convert_ins_to_dup,

            svimasm_ins_use_gaps = svimasm_ins_use_gaps,
            svimasm_ins_use_gaps_slack_bp = svimasm_ins_use_gaps_slack_bp,
            svimasm_ins_use_gaps_length_similarity = svimasm_ins_use_gaps_length_similarity,

            svimasm_ins_use_remap = svimasm_ins_use_remap,
            svimasm_ins_remap_max_length = svimasm_ins_remap_max_length,
            svimasm_ins_remap_cov_threshold = svimasm_ins_remap_cov_threshold,

            reference_fa = reference_fa,
            reference_fai = reference_fai,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 2-core, 8GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# UltralongSvimasmInsExtractDups            200%        70M           1s        
# bcftools sort                              70%        10M           1s
# bcftools concat                           100%        20M           1s
# truvari anno remap                        200%        20G          10m
# UltralongSvimasmInsExtractDupsPrime       200%        70M           1s
#
task Impl {
    input {
        File samples_tsv
        
        String remote_indir_svimasm
        String remote_outdir
        
        Int svimasm_min_sv_length
        Int svimasm_convert_ins_to_dup

        Int svimasm_ins_use_gaps
        Int svimasm_ins_use_gaps_slack_bp
        Float svimasm_ins_use_gaps_length_similarity

        Int svimasm_ins_use_remap
        Int svimasm_ins_remap_max_length
        Float svimasm_ins_remap_cov_threshold

        File reference_fa
        File reference_fai
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 32
        Int disk_size_gb = 20
        Int preemptible_number = 3
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
        


        
        # ---------------------- Steps of the pipeline -------------------------

        # Runs `bcftools concat`, making sure that symbolic ALTs are not over-
        # collapsed and that the output is sorted.
        #
        # Remark: the procedure deletes the input VCFs.
        #
        function Concat() {
            local SAMPLE_ID=$1
            local INPUT1_VCF_GZ=$2
            local INPUT2_VCF_GZ=$3
            local OUTPUT_SUFFIX=$4

            # Adding SVLEN to symbolic ALTs
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${INPUT1_VCF_GZ} | bgzip --compress-level 1 > ${SAMPLE_ID}_out1.vcf.gz
            rm -f ${INPUT1_VCF_GZ} ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_out1.vcf.gz
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${INPUT2_VCF_GZ} | bgzip --compress-level 1 > ${SAMPLE_ID}_out2.vcf.gz
            rm -f ${INPUT2_VCF_GZ} ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_out2.vcf.gz

            # Concatenating
            ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type v ${SAMPLE_ID}_out1.vcf.gz ${SAMPLE_ID}_out2.vcf.gz --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_out1.vcf.gz* ${SAMPLE_ID}_out2.vcf.gz*

            # Removing SVLEN from symbolic ALTs
            bcftools view --header-only ${SAMPLE_ID}_out.vcf --output ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf
            ${TIME_COMMAND} bcftools view --no-header ${SAMPLE_ID}_out.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf
            rm -f ${SAMPLE_ID}_out.vcf
            bcftools sort --output-type z ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf --output ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf.gz
        }




        # ---------------------------- Main program ----------------------------

        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2

        if [ ~{svimasm_convert_ins_to_dup} -eq 1 ]; then
            INSDUP_MODE="0"
        else
            INSDUP_MODE="1"
        fi
        cat ~{samples_tsv} | tr '\t' ',' > samples.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            DIPCALL_BED=$(echo ${LINE} | cut -d , -f 2)
            
            # Skipping the sample if it has already been processed or if there 
            # is no truth
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "1" )
            if [ ${TEST} != "1" ]; then
                continue
            fi
            TEST=$( gcloud storage ls ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.gz || echo "1" )
            if [ ${TEST} = "1" ]; then
                continue
            fi

            # Downloading the canonized svim-asm VCF
            gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.'gz*' .
            bcftools filter --include 'ABS(SVLEN)>='~{svimasm_min_sv_length} --output-type z ${SAMPLE_ID}_canonized.vcf.gz --output ${SAMPLE_ID}_svimasm.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm.vcf.gz
            rm -f ${SAMPLE_ID}_canonized.vcf.gz*

            # Extracting gaps from the dipcall BED
            # Remark: we do not handle gaps in the reference explicitly,
            # since we assume that calls in reference gaps have already been 
            # removed from the query VCFs upstream.
            N_GAPS="0"
            if [ ~{svimasm_ins_use_gaps} -eq 1 ]; then
                gcloud storage cp ${DIPCALL_BED} ./${SAMPLE_ID}_dipcall.bed
                ${TIME_COMMAND} bedtools sort -i ${SAMPLE_ID}_dipcall.bed -g ~{reference_fai} > ${SAMPLE_ID}_dipcall_sorted.bed
                ${TIME_COMMAND} bedtools complement -L -i ${SAMPLE_ID}_dipcall_sorted.bed -g ~{reference_fai} > ${SAMPLE_ID}_gaps.bed
                rm -f ${SAMPLE_ID}_dipcall.bed ${SAMPLE_ID}_dipcall_sorted.bed
                N_GAPS=$(wc -l < ${SAMPLE_ID}_gaps.bed)
                echo "The dipcall BED of ${SAMPLE_ID} has ${N_GAPS} gaps" 1>&2
            fi

            # 1. DEL
            bcftools filter --include 'SVTYPE=="DEL"' --output-type z ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_del.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_del.vcf.gz

            # 2. DUP
            # We include the intervals of both DUP:TANDEM and of DUP:INT, even
            # though the latter (interspersed duplication) behaves like DUP only
            # in terms of depth, not in terms of breakpoints.
            #
            # Remark: SVIM-asm outputs a `FORMAT:CN` field with the copy number 
            # of a DUP. We could use it to separate simple DUPs from 
            # triplications etc.: this is left to the future.
            bcftools filter --include 'SVTYPE=="DUP"' --output-type v ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_out.vcf
            java -cp ~{docker_dir} UltralongForceDup ${SAMPLE_ID}_out.vcf | bgzip > ${SAMPLE_ID}_svimasm_dup.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_dup.vcf.gz

            # 3. INS
            # Even though svim-asm explicitly emits DUPs, it might also 
            # emit INS records that are actually DUP. We reclassify such 
            # INS as DUP as follows.
            #
            # We reclassify as DUP every INS that is compatible with
            # a dipcall BED gap. Recall that a base is included in the 
            # dipcall BED iff: (1) it is covered by one alignment >=50kb 
            # mapQ>=5 in each haplotype; and (2) it is not covered by other 
            # >=10kb alignments in either haplotype. A simple DUP >10kb might 
            # satisfy (1) but it likely does not satisfy (2). This heuristic
            # is very fast.
            #
            # However, a complex DUP might satisfy (2), since its sub-
            # intervals might be covered by alignments <10kb, so the 
            # heuristic above would only mark simple DUPs as true. If 
            # svim-asm DUP records do not contain complex DUPs, then the 
            # latter would never be marked as true and the model would never
            # learn their BAM patterns.
            #
            # Thus, we use `truvari anno remap` on the remaining INS.
            # Working only on the remaining INS hopefully reduces the
            # runtime of this slow step. However, using two methods might 
            # give rise to two populations of DUPs with different breakpoint
            # decisions. Since these two populations are likely simple and
            # complex DUPs, which are likely to have distinct BAM patterns 
            # anyway, we accept this issue in exchange for the advantage in 
            # runtime.
            #
            # Remark: we assume that, for every DUP:INT, svim-asm outputs also
            # an INS, i.e. that it was run with `--interspersed_duplications_as_
            # insertions` too. I.e. the svim-asm VCF represents both the source
            # interval and the destination point of every DUP:INS.
            bcftools filter --include 'SVTYPE=="INS"' --output-type v ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_ins.vcf
            N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins.vcf | wc -l)
            if [ ${N_INS} -gt 0 ]; then
                # 3.1 Filtering all INS with the dipcall BED
                # Approx. 5% of all INS get marked as DUP in this step.
                if [ ~{svimasm_ins_use_gaps} -eq 1 -a ${N_GAPS} -gt 0 ]; then
                    ${TIME_COMMAND} java -cp ~{docker_dir} UltralongSvimasmInsExtractDups ${SAMPLE_ID}_svimasm_ins.vcf ${SAMPLE_ID}_gaps.bed $(wc -l < ${SAMPLE_ID}_gaps.bed) ~{svimasm_ins_use_gaps_slack_bp} ~{svimasm_ins_use_gaps_length_similarity} ${SAMPLE_ID}_svimasm_ins_dup.vcf ${SAMPLE_ID}_svimasm_ins_ins.vcf ${INSDUP_MODE}
                    rm -f ${SAMPLE_ID}_svimasm_ins.vcf ; mv ${SAMPLE_ID}_svimasm_ins_ins.vcf ${SAMPLE_ID}_svimasm_ins.vcf
                    N_INS_DUP=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins_dup.vcf | wc -l)
                    if [ ${N_INS_DUP} -gt 0 ]; then
                        ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_svimasm_ins_dup.vcf --output ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_ins_dup.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                    fi
                fi

                # 3.2 Filtering the remaining INS with `truvari anno remap`.
                # Approx. 55% of the remaining INS get marked as DUP.
                N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins.vcf | wc -l)
                if [ ~{svimasm_ins_use_remap} -eq 1 -a ${N_INS} -gt 0 ]; then
                    bgzip --compress-level 1 ${SAMPLE_ID}_svimasm_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins.vcf.gz
                    ${TIME_COMMAND} truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length ~{svimasm_ins_remap_max_length} --cov-threshold ~{svimasm_ins_remap_cov_threshold} -r ~{reference_fa} ${SAMPLE_ID}_svimasm_ins.vcf.gz -o ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz
                    rm -f ${SAMPLE_ID}_svimasm_ins.vcf.gz*
                    ${TIME_COMMAND} java -cp ~{docker_dir} UltralongSvimasmInsExtractDupsPrime ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf ${SAMPLE_ID}_svimasm_ins_ins.vcf ${INSDUP_MODE}
                    rm -f ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz* ; mv ${SAMPLE_ID}_svimasm_ins_ins.vcf ${SAMPLE_ID}_svimasm_ins.vcf
                    ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf --output ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz
                    rm -f ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz
                    if [ -e ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ]; then
                        Concat ${SAMPLE_ID} ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz svimasm_ins_dup
                    else
                        mv ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        mv ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz.tbi ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz.tbi
                    fi
                fi
            fi
            bgzip --compress-level 1 ${SAMPLE_ID}_svimasm_ins.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins.vcf.gz
            if [ ! -e ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ]; then
                bcftools view --header-only --output-type z ${SAMPLE_ID}_svimasm_ins.vcf.gz --output ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
            fi

            # 4. INV
            bcftools filter --include 'SVTYPE=="INV"' --output-type z ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_inv.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_inv.vcf.gz

            # Uploading and deallocating the sample
            gcloud storage mv ${SAMPLE_ID}_'*.vcf.gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -rf ${SAMPLE_ID}_*
            ls -laht 1>&2
        done 3< samples.csv
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
