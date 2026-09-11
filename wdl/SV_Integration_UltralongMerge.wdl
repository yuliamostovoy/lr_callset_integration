version 1.0


# Merges the per-sample, per-SVTYPE VCFs produced by 
# `SV_Integration_{Ultralong,Bnd}Annotate.wdl` and 
# `SV_Integration_{Ultralong,Bnd}GetTrainingIntervals.wdl` into a single, sorted
# VCF, without doing any form of collapse, and overwriting ALT in an ad hoc, 
# non-standard way that is necessary for running VETS correctly downstream.
#
workflow SV_Integration_UltralongMerge {
    input {
        File samples_csv
        String remote_indir
        String remote_outdir
        
        String svtype
        String suffix
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            samples_csv = samples_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            svtype = svtype,
            suffix = suffix,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 4-core, 4GB VM, all HPRC+HGSVC samples:
#
# TOOL                                                CPU     RAM     TIME
# fix_sample.sh                                      
# bcftools concat                                    
#
task Impl {
    input {
        File samples_csv
        String remote_indir
        String remote_outdir
        
        String svtype
        String suffix
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 20
    }
    parameter_meta {
        samples_csv: "Format: ID, ..."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        cat << 'END' > preprocess_single_sample.sh
#!/bin/bash
set -euxo pipefail

DOCKER_DIR=$1
SVTYPE=$2
INPUT_VCF_GZ=$3

# Forcing a fixed artificial sample name
bcftools reheader --samples-list SAMPLE ${INPUT_VCF_GZ} --output ${INPUT_VCF_GZ}_reheader.vcf.gz
rm -f ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi ; mv ${INPUT_VCF_GZ}_reheader.vcf.gz ${INPUT_VCF_GZ}

# Canonizing BNDs. This is necessary for the annotated BND VCFs, since they are 
# not yet canonized. It is redundant for the TPs VCFs, since they are already
# canonized.
if [ ${SVTYPE} = "bnd" -o ${SVTYPE} = "BND" ]; then
    java -cp ${DOCKER_DIR} BndCanonize ${INPUT_VCF_GZ} | bcftools sort - --output-type z > ${INPUT_VCF_GZ}_canonized.vcf.gz
    rm -f ${INPUT_VCF_GZ} ; mv ${INPUT_VCF_GZ}_canonized.vcf.gz ${INPUT_VCF_GZ}
fi

# Creating an artificial ALT that contains ID (we assume that all IDs are unique
# both at the intra-sample level and at the inter-sample level: this is 
# enforced by the steps upstream). This is necessary to make 
# `ExtractVariantAnnotations --resource-matching-strategy START_POSITION_AND_
# GIVEN_REPRESENTATION` work downstream, and it prevents `bcftools concat` to
# collapse any record.
java -cp ${DOCKER_DIR} SetAltToID ${INPUT_VCF_GZ} | bgzip --compress-level 1 > ${INPUT_VCF_GZ}_prime
rm -f ${INPUT_VCF_GZ} ; mv ${INPUT_VCF_GZ}_prime ${INPUT_VCF_GZ}

bcftools index -f -t ${INPUT_VCF_GZ}
END
        chmod +x preprocess_single_sample.sh
        


        
        # ---------------------------- Main program ----------------------------
        
        # Concatenation without any form of collapse, since here we are only 
        # interested in preserving the original intra-sample feature values for
        # training the models downstream. The rest of a VCF record is just a key
        # used to match the resource to the input by ExtractVariantAnnotations.
        # Removing duplicates may discard information if e.g. one instance is 
        # true in one sample and false in another sample (with different feature
        # values in each).
        rm -f list.txt
        while read LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            echo ~{remote_indir}/"${SAMPLE_ID}_"~{svtype}~{suffix}'.vcf.gz' >> list.txt
            echo ~{remote_indir}/"${SAMPLE_ID}_"~{svtype}~{suffix}'.vcf.gz.tbi' >> list.txt
        done < ~{samples_csv}
        cat list.txt | gcloud storage cp -I .
        df -h 1>&2
        ls -laht 1>&2
        ls *.vcf.gz > list.txt
        ${TIME_COMMAND} xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} ./preprocess_single_sample.sh ~{docker_dir} ~{svtype}
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type v --output ~{svtype}~{suffix}_merged.vcf

        # Sorting and uploading
        ${TIME_COMMAND} bcftools sort --output-type z ~{svtype}~{suffix}_merged.vcf --output ~{svtype}~{suffix}_merged.vcf.gz
        rm -f ~{svtype}~{suffix}_merged.vcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ~{svtype}~{suffix}_merged.vcf.gz
        ls -laht 1>&2
        ${TIME_COMMAND} gcloud storage mv ~{svtype}~{suffix}'*_merged.vcf.gz*' ~{remote_outdir}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
