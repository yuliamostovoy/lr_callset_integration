version 1.0


# Given a set of annotated, per-sample, per-SVTYPE query VCFs, and their 
# corresponding svim-asm truth VCFs and dipcall BEDs, the program selects, 
# for each sample, the query records that are likely true according to the 
# assemblies.
#
# Remark: sequence similarity is never used to compute matches.
#
workflow SV_Integration_UltralongGetTrainingIntervals {
    input {
        File samples_tsv
        
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir

        Int convert_ins_to_dup = 1

        Int match_to_gaps = 0
        Int match_insdups_to_dups = 1
        Int match_dups_to_insdups = 0
        Int match_ins_to_dup = 1
        Int match_ins_to_insdup = 0
        Int match_ins_to_dup_slack_bp = 200

        File reference_fai
        
        Int truvari_refdist = 500
        Float truvari_pctsize_loose = 0.4
        Float truvari_pctsize_strict = 0.9
        Float truvari_pctovl_loose = 0.4

        Int max_read_length = 25000
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID, DIPCALL_BED"
        remote_indir_query: "Without final slash. Contains per-sample annotated VCFs created by `SV_Integration_UltralongAnnotate.wdl`."
        remote_indir_svimasm: "Without final slash. Contains per-sample canonized and filtered svim-asm VCFs."
        convert_ins_to_dup: "1=SVIM-asm's INS records that correspond to duplications are rewritten as DUP records"
        truvari_pctsize_loose: "Should be set for loose matches to capture somewhat similar calls, since breakpoint distance is explicitly checked after `truvari bench` when appropriate."
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,

            remote_indir_query = remote_indir_query,
            remote_indir_svimasm = remote_indir_svimasm,
            remote_outdir = remote_outdir,

            convert_ins_to_dup = convert_ins_to_dup,

            match_to_gaps = match_to_gaps,
            match_insdups_to_dups = match_insdups_to_dups,
            match_dups_to_insdups = match_dups_to_insdups,
            match_ins_to_dup = match_ins_to_dup,
            match_ins_to_insdup = match_ins_to_insdup,
            match_ins_to_dup_slack_bp = match_ins_to_dup_slack_bp,

            reference_fai = reference_fai,

            truvari_refdist = truvari_refdist,
            truvari_pctsize_loose = truvari_pctsize_loose,
            truvari_pctsize_strict = truvari_pctsize_strict,
            truvari_pctovl_loose = truvari_pctovl_loose,

            max_read_length = max_read_length,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 2-core, 8GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# bcftools sort                              70%        10M           1s
# bcftools concat                           100%        20M           1s
# truvari bench                             100%       100M           1s
# UltralongBed2IntervalVcf                  200%        50M           1s
#
task Impl {
    input {
        File samples_tsv
        
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir

        Int convert_ins_to_dup

        Int match_to_gaps
        Int match_insdups_to_dups
        Int match_dups_to_insdups
        Int match_ins_to_dup
        Int match_ins_to_insdup
        Int match_ins_to_dup_slack_bp

        File reference_fai

        Int truvari_refdist
        Float truvari_pctsize_loose
        Float truvari_pctsize_strict
        Float truvari_pctovl_loose

        Int max_read_length
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 4
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
            rm -f ${INPUT1_VCF_GZ}* ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_out1.vcf.gz
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${INPUT2_VCF_GZ} | bgzip --compress-level 1 > ${SAMPLE_ID}_out2.vcf.gz
            rm -f ${INPUT2_VCF_GZ}* ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_out2.vcf.gz

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
            rm -f ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${OUTPUT_SUFFIX}.vcf.gz
        }


        # Runs `truvari bench`, additionally making sure that the output is 
        # sorted and optionally enforcing a max distance between each 
        # corresponding pair of breakpoints.
        #
        # Remark: `StartDistance=EndDistance` for INS.
        #
        function Bench() {
            local SAMPLE_ID=$1
            local PCTSIZE=$2
            local PCTOVL=$3
            local REFDIST=$4
            local ENFORCE_BREAKPOINT_DISTANCE=$5
            local QUERY_VCF_GZ=$6
            local TRUTH_VCF_GZ=$7
            local OUTPUT_VCF_GZ=$8
            local PRINT_REMAP_PERC_VALUES=$9

            local DEFAULT_TRUVARI_CHUNKSIZE="1000"
            if [ ${REFDIST} -gt ${DEFAULT_TRUVARI_CHUNKSIZE} ]; then
                CHUNKSIZE=${REFDIST}
            else
                CHUNKSIZE=${DEFAULT_TRUVARI_CHUNKSIZE}
            fi
            ${TIME_COMMAND} truvari bench -b ${TRUTH_VCF_GZ} -c ${QUERY_VCF_GZ} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ${REFDIST} --chunksize ${CHUNKSIZE} --pctseq 0 --pctsize ${PCTSIZE} --pctovl ${PCTOVL} --pick single -o ./${SAMPLE_ID}_truvari/
            if [ ${ENFORCE_BREAKPOINT_DISTANCE} -eq 1 ]; then
                bcftools filter --include 'ABS(StartDistance)<='${REFDIST}' && ABS(EndDistance)<='${REFDIST} --output-type z ${SAMPLE_ID}_truvari/tp-comp.vcf.gz --output ${OUTPUT_VCF_GZ}
            else
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${OUTPUT_VCF_GZ}
            fi
            ${TIME_COMMAND} bcftools sort --output-type z ${OUTPUT_VCF_GZ} --output ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${OUTPUT_VCF_GZ}
            bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_VCF_GZ}
            if [ ${PRINT_REMAP_PERC_VALUES} -eq 1 ]; then
                bcftools query --format '%INFO/remap_classification\t%INFO/remap_perc\n' ${SAMPLE_ID}_truvari/tp-base.vcf.gz 1>&2
            fi
            rm -rf ${SAMPLE_ID}_truvari/
        }




        # ---------------------------- Main program ----------------------------

        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2

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
            TEST=$( gcloud storage ls ~{remote_indir_svimasm}/${SAMPLE_ID}.done || echo "1" )
            if [ ${TEST} = "1" ]; then
                continue
            fi
            gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_'*.vcf.gz*' .
            gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_'*.vcf.gz*' .

            # Extracting gaps from the dipcall BED
            # Remark: we do not handle gaps in the reference explicitly,
            # since we assume that calls in reference gaps have already been 
            # removed from the query VCFs upstream.
            N_GAPS="0"
            if [ ~{match_to_gaps} -eq 1 ]; then
                gcloud storage cp ${DIPCALL_BED} ./${SAMPLE_ID}_dipcall.bed
                ${TIME_COMMAND} bedtools sort -i ${SAMPLE_ID}_dipcall.bed -g ~{reference_fai} > ${SAMPLE_ID}_dipcall_sorted.bed
                ${TIME_COMMAND} bedtools complement -L -i ${SAMPLE_ID}_dipcall_sorted.bed -g ~{reference_fai} > ${SAMPLE_ID}_gaps.bed
                rm -f ${SAMPLE_ID}_dipcall.bed ${SAMPLE_ID}_dipcall_sorted.bed
                N_GAPS=$(wc -l < ${SAMPLE_ID}_gaps.bed)
                echo "The dipcall BED of ${SAMPLE_ID} has ${N_GAPS} gaps" 1>&2
            fi

            # 1. DEL
            # Remark: svim-asm does not emit calls at some query DELs that 
            # correspond to gaps in the dipcall BED. DELs >10kb are likely
            # to create gaps in dipcall's BED by the definition of confident
            # BED, since such DELs induce split alignments (see above). However, 
            # marking as true every query DEL that matches a gap in dipcall's 
            # BED does not improve performance, and it decreases it both within
            # the confident BED and on the whole genome.
            # Approx. 8% of all DEL get marked as true by a dipcall gap.
            Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ~{truvari_refdist} 1 ${SAMPLE_ID}_del.vcf.gz ${SAMPLE_ID}_svimasm_del.vcf.gz ${SAMPLE_ID}_del1.vcf.gz 0        
            if [ ~{match_to_gaps} -eq 1 -a ${N_GAPS} -gt 0 ]; then
                bcftools view --header-only ${SAMPLE_ID}_del.vcf.gz > ${SAMPLE_ID}_gaps.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBed2IntervalVcf ${SAMPLE_ID}_gaps.bed DEL >> ${SAMPLE_ID}_gaps.vcf
                bgzip -@ ${N_THREADS} ${SAMPLE_ID}_gaps.vcf
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_gaps.vcf.gz
                Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ~{truvari_refdist} 1 ${SAMPLE_ID}_del.vcf.gz ${SAMPLE_ID}_gaps.vcf.gz ${SAMPLE_ID}_del2.vcf.gz 0
                N_DEL1=$(bcftools index --nrecords ${SAMPLE_ID}_del1.vcf.gz)
                N_DEL2=$(bcftools index --nrecords ${SAMPLE_ID}_del2.vcf.gz)
                if [ ${N_DEL1} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_del1.vcf.gz*
                    mv ${SAMPLE_ID}_del2.vcf.gz ${SAMPLE_ID}_del_training.vcf.gz
                    mv ${SAMPLE_ID}_del2.vcf.gz.tbi ${SAMPLE_ID}_del_training.vcf.gz.tbi
                elif [ ${N_DEL2} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_del2.vcf.gz*
                    mv ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del_training.vcf.gz
                    mv ${SAMPLE_ID}_del1.vcf.gz.tbi ${SAMPLE_ID}_del_training.vcf.gz.tbi
                else
                    Concat ${SAMPLE_ID} ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del2.vcf.gz del_training
                fi
                rm -f ${SAMPLE_ID}_gaps.vcf.gz*
            else
                mv ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del_training.vcf.gz
                mv ${SAMPLE_ID}_del1.vcf.gz.tbi ${SAMPLE_ID}_del_training.vcf.gz.tbi
            fi

            # 2. DUP
            # Remark: DUP records from svim-asm seem to be generally 
            # comprehensive, and they do not need to be augmented with e.g. gaps
            # in dipcall's BED.
            Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ~{truvari_refdist} 1 ${SAMPLE_ID}_dup.vcf.gz ${SAMPLE_ID}_svimasm_dup.vcf.gz ${SAMPLE_ID}_dup_training.vcf.gz 0
            if [ ~{convert_ins_to_dup} -eq 1 -a ~{match_dups_to_insdups} -eq 1 ]; then
                # Remark: records that are originally called as DUP on both 
                # query and truth are likely enriched in simple DUPs, whereas
                # INS->DUP records are likely enriched in complex DUPs. However,
                # it could still happen that some truth INSDUPs represent simple
                # DUPs.
                Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ~{truvari_refdist} 1 ${SAMPLE_ID}_dup.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ${SAMPLE_ID}_dup_training_prime.vcf.gz 0
                N_DUP_TRAINING=$(bcftools index --nrecords ${SAMPLE_ID}_dup_training.vcf.gz)
                N_DUP_TRAINING_PRIME=$(bcftools index --nrecords ${SAMPLE_ID}_dup_training_prime.vcf.gz)
                if [ ${N_DUP_TRAINING} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_dup_training.vcf.gz*
                    mv ${SAMPLE_ID}_dup_training_prime.vcf.gz ${SAMPLE_ID}_dup_training.vcf.gz
                    mv ${SAMPLE_ID}_dup_training_prime.vcf.gz.tbi ${SAMPLE_ID}_dup_training.vcf.gz.tbi
                elif [ ${N_DUP_TRAINING_PRIME} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_dup_training_prime.vcf.gz*
                else
                    Concat ${SAMPLE_ID} ${SAMPLE_ID}_dup_training.vcf.gz ${SAMPLE_ID}_dup_training_prime.vcf.gz dup_training
                fi
            fi

            # 3. INSDUP
            # Remark: we do not enforce breakpoint distance in INSDUP-INSDUP 
            # comparisons, since the query breakpoints computed with the depth 
            # heuristic might be very different from the truth ones, computed 
            # with `truvari anno remap`. Hopefully the query breakpoints can 
            # still create useful feature vectors.
            if [ ~{convert_ins_to_dup} -eq 1 ]; then
                TRUVARI_REFDIST=~{truvari_refdist}
                ENFORCE_BREAKPOINT_DISTANCE="0"
            else
                TRUVARI_REFDIST="2000000"
                ENFORCE_BREAKPOINT_DISTANCE="1"
            fi
            Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ${TRUVARI_REFDIST} ${ENFORCE_BREAKPOINT_DISTANCE} ${SAMPLE_ID}_insdup.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ${SAMPLE_ID}_insdup_training.vcf.gz 1
            if [ ~{convert_ins_to_dup} -eq 1 -a ~{match_insdups_to_dups} -eq 1 ]; then
                # We compare the query INS->DUPs to svim-asm's DUPs, since the
                # former might contain simple DUPs.
                # Remark: in practice there are few matches to svim-asm's DUPs
                # so this comparison might be irrelevant. This is probably due
                # to the fact that we truvari-collapsed query INSDUPs to query
                # DUPs, favoring the latter representation.
                Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ${TRUVARI_REFDIST} 1 ${SAMPLE_ID}_insdup.vcf.gz ${SAMPLE_ID}_svimasm_dup.vcf.gz ${SAMPLE_ID}_insdup_training_prime.vcf.gz 0
                N_INSDUP_TRAINING=$(bcftools index --nrecords ${SAMPLE_ID}_insdup_training.vcf.gz)
                N_INSDUP_TRAINING_PRIME=$(bcftools index --nrecords ${SAMPLE_ID}_insdup_training_prime.vcf.gz)
                if [ ${N_INSDUP_TRAINING} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_insdup_training.vcf.gz*
                    mv ${SAMPLE_ID}_insdup_training_prime.vcf.gz ${SAMPLE_ID}_insdup_training.vcf.gz
                    mv ${SAMPLE_ID}_insdup_training_prime.vcf.gz.tbi ${SAMPLE_ID}_insdup_training.vcf.gz.tbi
                elif [ ${N_INSDUP_TRAINING_PRIME} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_insdup_training_prime.vcf.gz*
                else
                    Concat ${SAMPLE_ID} ${SAMPLE_ID}_insdup_training.vcf.gz ${SAMPLE_ID}_insdup_training_prime.vcf.gz insdup_training
                fi
            fi
            
            # 4. INS
            # Some query INS might be short DUPs (simple or complex) that are 
            # represented as CIGAR INS operations in the BAM and do not change
            # the local depth (i.e. their BAM pattern is identical to the one of
            # short INS). They may correspond to a DUP or INSDUP in svim-asm,
            # and if we do not mark them as true, the model may not learn their
            # BAM pattern if there are no other short INS.
            #
            # Remark: in practice matching query INS to truth DUP and INSDUP
            # does not improve performance. There are few matches per sample, 
            # mostly with truth INSDUPs, which might introduce spurious BAM 
            # patterns in the training set.
            #
            # Remark: it is not useful to convert to intervals short DUPs that 
            # are represented as INS records in the query VCF, since there is
            # likely no BAM pattern at the interval's boundaries.
            Bench ${SAMPLE_ID} ~{truvari_pctsize_strict} 0 ~{truvari_refdist} 1 ${SAMPLE_ID}_ins.vcf.gz ${SAMPLE_ID}_svimasm_ins.vcf.gz ${SAMPLE_ID}_out1.vcf.gz 0
            CONCAT_STRING="${SAMPLE_ID}_out1.vcf.gz"
            if [ ~{match_ins_to_dup} -eq 1 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongMatchInsToDup ${SAMPLE_ID}_ins.vcf.gz ${SAMPLE_ID}_svimasm_dup.vcf.gz $(bcftools index --nrecords ${SAMPLE_ID}_svimasm_dup.vcf.gz) ~{truvari_pctsize_strict} ~{match_ins_to_dup_slack_bp} ~{max_read_length} | bgzip --compress-level 1 > ${SAMPLE_ID}_out2.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_out2.vcf.gz
                CONCAT_STRING="${CONCAT_STRING} ${SAMPLE_ID}_out2.vcf.gz"
            fi
            if [ ~{match_ins_to_insdup} -eq 1 -a ~{convert_ins_to_dup} -eq 1 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongMatchInsToDup ${SAMPLE_ID}_ins.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz $(bcftools index --nrecords ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz) ~{truvari_pctsize_loose} ~{match_ins_to_dup_slack_bp} ~{max_read_length} | bgzip --compress-level 1 > ${SAMPLE_ID}_out3.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_out3.vcf.gz
                CONCAT_STRING="${CONCAT_STRING} ${SAMPLE_ID}_out3.vcf.gz"
            fi
            if [ ~{match_ins_to_dup} -eq 1 -o ~{match_ins_to_insdup} -eq 1 ]; then
                # No need to call function `Concat()` here, since no ALT is
                # symbolic.
                ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type v ${CONCAT_STRING} | bcftools sort --output-type z --output ${SAMPLE_ID}_ins_training.vcf.gz
                rm -f ${SAMPLE_ID}_out*.vcf.gz* ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins_training.vcf.gz
            else
                mv ${SAMPLE_ID}_out1.vcf.gz ${SAMPLE_ID}_ins_training.vcf.gz
                mv ${SAMPLE_ID}_out1.vcf.gz.tbi ${SAMPLE_ID}_ins_training.vcf.gz.tbi
            fi

            # 5. INV
            # Remark: svim-asm does not seem to emit calls at some query INVs
            # that correspond to gaps in the dipcall BED, so it could make sense
            # to mark as true such query INVs. In practice this degrades 
            # performance, probably because it adds to the training set several
            # events that are not simple INVs.
            # Approx. 15% of all INV get marked as true by a dipcall gap.
            Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ~{truvari_refdist} 1 ${SAMPLE_ID}_inv.vcf.gz ${SAMPLE_ID}_svimasm_inv.vcf.gz ${SAMPLE_ID}_inv1.vcf.gz 0
            if [ ~{match_to_gaps} -eq 1 -a ${N_GAPS} -gt 0 ]; then
                bcftools view --header-only ${SAMPLE_ID}_inv.vcf.gz > ${SAMPLE_ID}_gaps.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBed2IntervalVcf ${SAMPLE_ID}_gaps.bed INV >> ${SAMPLE_ID}_gaps.vcf
                bgzip -@ ${N_THREADS} ${SAMPLE_ID}_gaps.vcf
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_gaps.vcf.gz
                Bench ${SAMPLE_ID} ~{truvari_pctsize_loose} ~{truvari_pctovl_loose} ~{truvari_refdist} 1 ${SAMPLE_ID}_inv.vcf.gz ${SAMPLE_ID}_gaps.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz 0
                N_INV1=$(bcftools index --nrecords ${SAMPLE_ID}_inv1.vcf.gz)
                N_INV2=$(bcftools index --nrecords ${SAMPLE_ID}_inv2.vcf.gz)
                if [ ${N_INV1} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_inv1.vcf.gz*
                    mv ${SAMPLE_ID}_inv2.vcf.gz ${SAMPLE_ID}_inv_training.vcf.gz
                    mv ${SAMPLE_ID}_inv2.vcf.gz.tbi ${SAMPLE_ID}_inv_training.vcf.gz.tbi
                elif [ ${N_INV2} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_inv2.vcf.gz*
                    mv ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv_training.vcf.gz
                    mv ${SAMPLE_ID}_inv1.vcf.gz.tbi ${SAMPLE_ID}_inv_training.vcf.gz.tbi
                else
                    Concat ${SAMPLE_ID} ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz inv_training
                fi
                rm -f ${SAMPLE_ID}_gaps.vcf.gz*
            else
                mv ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv_training.vcf.gz
                mv ${SAMPLE_ID}_inv1.vcf.gz.tbi ${SAMPLE_ID}_inv_training.vcf.gz.tbi
            fi
            
            # Uploading and deallocating the sample
            gcloud storage mv ${SAMPLE_ID}_'*_training.vcf.gz*' ~{remote_outdir}/
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
