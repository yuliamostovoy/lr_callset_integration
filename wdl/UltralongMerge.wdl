version 1.0


# Merges the per-sample, per-SVTYPE VCFs produced by 
# `SV_Integration_UltralongAnnotate.wdl` and 
# `SV_Integration_UltralongGetTrainingIntervals.wdl` into a single, sorted VCF,
# removing only exact duplicates. The input VCFs are assumed to have standard 
# symbolic ALT fields, e.g. ALT=<DEL> with SVLEN=23456 in INFO, and not 
# ALT=<DEL-23456>.
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
# fix_sample.sh                                      300%     11M       4s
# bcftools concat                                    300%    1.5G      20s
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
        
        cat << 'END' > fix_sample.sh
#!/bin/bash
DOCKER_DIR=$1
SVTYPE=$2
INPUT_VCF_GZ=$3

bcftools reheader --samples-list SAMPLE ${INPUT_VCF_GZ} --output ${INPUT_VCF_GZ}_reheader.vcf.gz
rm -f ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi

# Adding SVLEN to symbolic ALTs, to avoid overcollapse by `bcftools concat`
# downstream.
if [ ${SVTYPE} != "ins" -a ${SVTYPE} != "INS" -a ${SVTYPE} != "bnd" -a ${SVTYPE} != "BND" ]; then
    java -cp ${DOCKER_DIR} AddSvlenToSymbolicAlt ${INPUT_VCF_GZ}_reheader.vcf.gz | bgzip --compress-level 1 > ${INPUT_VCF_GZ}
    rm -f ${INPUT_VCF_GZ}_reheader.vcf.gz
else
    mv ${INPUT_VCF_GZ}_reheader.vcf.gz ${INPUT_VCF_GZ}
fi

bcftools index -f -t ${INPUT_VCF_GZ}
END
        chmod +x fix_sample.sh
        
        
        # ---------------------------- Main program ----------------------------
        
        # Simple concatenation, with only exact duplicate removal. In the
        # future we may run truvari collapse to remove approximate duplicates.
        rm -f list.txt
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            echo ~{remote_indir}/"${SAMPLE_ID}_"~{svtype}~{suffix}'.vcf.gz' >> list.txt
            echo ~{remote_indir}/"${SAMPLE_ID}_"~{svtype}~{suffix}'.vcf.gz.tbi' >> list.txt
        done < ~{samples_csv}
        cat list.txt | gcloud storage cp -I .
        df -h 1>&2
        ls -laht 1>&2
        ls *.vcf.gz > list.txt
        ${TIME_COMMAND} xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} ./fix_sample.sh ~{docker_dir} ~{svtype}
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --file-list list.txt --output-type v --output out.vcf

        # Removing SVLEN from symbolic ALTs
        if [ ~{svtype} != "ins" -a ~{svtype} != "INS" -a ~{svtype} != "bnd" -a ~{svtype} != "BND" ]; then
            bcftools view --header-only out.vcf --output ~{svtype}~{suffix}_merged.vcf
            ${TIME_COMMAND} bcftools view --no-header out.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> ~{svtype}~{suffix}_merged.vcf
            rm -f out.vcf
        else
            mv out.vcf ~{svtype}~{suffix}_merged.vcf
        fi
        ${TIME_COMMAND} bcftools sort --output-type z ~{svtype}~{suffix}_merged.vcf --output ~{svtype}~{suffix}_merged.vcf.gz
        rm -f ~{svtype}~{suffix}_merged.vcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ~{svtype}~{suffix}_merged.vcf.gz
        ls -laht 1>&2
        
        # Uploading
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
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
