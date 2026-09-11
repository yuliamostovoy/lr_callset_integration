version 1.0


# A simplified version of `SV_Integration_UltralongGetTrainingIntervals.wdl`.
#
workflow SV_Integration_BndGetTrainingIntervals {
    input {
        File samples_tsv
        Int remove_orientations = 1

        Int use_interval_svtypes = 0
        Int interval_svtypes_min_sv_length = 1000

        Int truvari_bnddist = 500
        
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID, ???"
        remote_indir_query: "Without final slash. Contains per-sample annotated VCFs created by `SV_Integration_BndAnnotate.wdl`."
        remote_indir_svimasm: "Without final slash. Contains per-sample canonized and filtered svim-asm VCFs."
        interval_svtypes_min_sv_length: "Only simple, non-INS SVs with SVLEN>=this are kept in the truth VCF."
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,
            remove_orientations = remove_orientations,

            use_interval_svtypes = use_interval_svtypes,
            interval_svtypes_min_sv_length = interval_svtypes_min_sv_length,

            truvari_bnddist = truvari_bnddist,

            remote_indir_query = remote_indir_query,
            remote_indir_svimasm = remote_indir_svimasm,
            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 1-core, 4GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# truvari bench                             100%       ????         ????
#
task Impl {
    input {
        File samples_tsv
        Int remove_orientations

        Int use_interval_svtypes
        Int interval_svtypes_min_sv_length

        Int truvari_bnddist
        
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir
        
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


        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2

        cat ~{samples_tsv} | tr '\t' ',' > samples.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "1" )
            if [ "${TEST}" != "1" ]; then
                continue
            fi

            # Downloading and filtering the truth VCF (skipping the sample if 
            # there is no truth).
            # Remark: we use both BND calls and ultralong calls from the truth,
            # since read-derived BNDs may describe breakpoints of variants that 
            # are too large for read-based callers but that are short enough for
            # svim-asm.
            TEST=$( gcloud storage ls ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.gz || echo "1" )
            if [ "${TEST}" = "1" ]; then
                rm -rf ${SAMPLE_ID}_*
                continue
            fi
            gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.'gz*' .
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'SVTYPE="BND"' --output-type v ${SAMPLE_ID}_canonized.vcf.gz --output ${SAMPLE_ID}_svimasm_bnd.vcf
            java -cp ~{docker_dir} BndCanonize ${SAMPLE_ID}_svimasm_bnd.vcf > ${SAMPLE_ID}_svimasm_bnd_canonized.vcf
            rm -f ${SAMPLE_ID}_svimasm_bnd.vcf
            if [ ~{remove_orientations} -eq 1 ]; then
                java -cp ~{docker_dir} BndRemoveOrientations ${SAMPLE_ID}_svimasm_bnd_canonized.vcf | bcftools sort - --output-type z > ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz
            else
                bcftools sort --output-type z ${SAMPLE_ID}_svimasm_bnd_canonized.vcf > ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz
            fi
            rm -f ${SAMPLE_ID}_svimasm_bnd_canonized.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz
            if [ ~{use_interval_svtypes} -eq 1 ]; then
                ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include '(SVTYPE="DEL" || SVTYPE="DUP" || SVTYPE="INV") && ABS(SVLEN)>='~{interval_svtypes_min_sv_length} --output-type z ${SAMPLE_ID}_canonized.vcf.gz --output ${SAMPLE_ID}_svimasm_ultralong.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ultralong.vcf.gz
                ${TIME_COMMAND} bcftools concat --allow-overlaps --output-type z ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz ${SAMPLE_ID}_svimasm_ultralong.vcf.gz --output ${SAMPLE_ID}_svimasm.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm.vcf.gz
                rm -f ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz* ${SAMPLE_ID}_svimasm_ultralong.vcf.gz*
            else
                mv ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz ${SAMPLE_ID}_svimasm.vcf.gz
                mv ${SAMPLE_ID}_svimasm_bnd_canonized.vcf.gz.tbi ${SAMPLE_ID}_svimasm.vcf.gz.tbi
            fi
            rm -f ${SAMPLE_ID}_canonized.vcf.gz*

            # Computing matches between query and truth VCF.
            # Remarks:
            # 1. Truvari bench decomposes a simple SV into its BNDs 
            #    automatically.
            # 2. We use `--pick multi` since an SV can be decomposed into
            #    multiple BNDs, in which case it would be matched only to one
            #    of them by default.
            # See https://github.com/acenglish/truvari/wiki/bench#cross-representation-matching
            gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_bnd.vcf.'gz*' .
            java -cp ~{docker_dir} BndCanonize ${SAMPLE_ID}_bnd.vcf.gz > ${SAMPLE_ID}_bnd_canonized.vcf
            rm -f ${SAMPLE_ID}_bnd.vcf.gz*
            if [ ~{remove_orientations} -eq 1 ]; then
                java -cp ~{docker_dir} BndRemoveOrientations ${SAMPLE_ID}_bnd_canonized.vcf | bcftools sort - --output-type z > ${SAMPLE_ID}_bnd_canonized.vcf.gz
            else
                bcftools sort --output-type z ${SAMPLE_ID}_bnd_canonized.vcf > ${SAMPLE_ID}_bnd_canonized.vcf.gz
            fi
            rm -f ${SAMPLE_ID}_bnd_canonized.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_bnd_canonized.vcf.gz
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm.vcf.gz -c ${SAMPLE_ID}_bnd_canonized.vcf.gz --bnddist ~{truvari_bnddist} --pick multi -o ./${SAMPLE_ID}_truvari/
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_truvari/tp-comp.vcf.gz --output ${SAMPLE_ID}_bnd_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_bnd_training.vcf.gz
            rm -rf ${SAMPLE_ID}_truvari/
            
            # Uploading and deallocating the sample
            gcloud storage mv ${SAMPLE_ID}_bnd_training.vcf.'gz*' ~{remote_outdir}/
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
