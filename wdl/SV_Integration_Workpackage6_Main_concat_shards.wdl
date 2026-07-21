version 1.0


# Concatenates the `truvari collapse` chunks. Ensures that every record in
# output has a globally unique ID (this is necessary for kanpig downstream;
# duplicated IDs may arise naturally from the previous steps of the pipeline),
# and an INFO field that counts the number of samples it was discovered in.
#
workflow SV_Integration_Workpackage6 {
    input {
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        String remote_indir
        String remote_outdir
        Int concat_all_naive = 1
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        chromosomes: "The order of the chromosomes becomes their order in the output VCF."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        concat_all_naive: "Concatenate chromosomes in a naive (1, default) or non-naive (0) way. Non-naive is necessary when different chromosomes were built by different versions of the pipeline, with slightly different code, and their headers are not exactly identical. Intra-chromosome concatenation is always performed in a naive way."
    }
    
    scatter (chr in chromosomes) {
        call SingleChromosome {
            input:
                chromosome = chr,
                remote_indir = remote_indir,
                remote_outdir = remote_outdir,
                docker_image = docker_image
        }
    }
    call AllChromosomes {
        input:
            chromosomes = chromosomes,
            out_txt = SingleChromosome.out_txt,
            remote_outdir = remote_outdir,
            naive = concat_all_naive,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, HDD:
#
# TOOL                           CPU     RAM   TIME
# download                      250%    100M    15s
# bcftools concat                30%     20M     1m
# bcftools query                100%    150M     4m
# bcftools annotate             100%    150M    10m
#
task SingleChromosome {
    input {
        String chromosome
        String remote_indir
        String remote_outdir

        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 200
        Int preemptible_number = 4
        Array[String]? upstream_signal
    }
    parameter_meta {
        upstream_signal: "Ordering-only handshake for orchestrator workflows; ignored by standalone runs."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # ---------------------------- Main program ----------------------------
        
        TEST=$( gcloud storage ls ~{remote_outdir}/~{chromosome}/~{chromosome}.done || echo "0" )
        if [ ${TEST} != "0" ]; then
            # Skipping the chromosome if it has already been processed
            :
        else
            # Localizing all chunks
            gcloud storage ls ~{remote_indir}/~{chromosome}/chunk_'*.bcf' > test.txt
            if grep -q '.bcf' test.txt ; then
                :
            else
                echo "ERROR: ~{chromosome} has no truvari collapse chunks."
                exit
            fi
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/~{chromosome}/chunk_'*.bcf*' .
            ls chunk_*.bcf | sort -V > chunk_list.txt
            cat chunk_list.txt
            df -h 1>&2
        
            # Concatenating all chunks
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type b --output out.bcf
            df -h 1>&2
            rm -rf chunk_* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            
            # Enforcing a distinct ID in every record, and annotating every
            # record with the number of samples it occurs in. Note that the
            # latter is not equal to the QUAL field in input to truvari collapse
            # upstream, so we have to recompute this number.
            CHR=~{chromosome}
            CHR=${CHR#chr}
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' in.bcf | awk -v id=${CHR} 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=sprintf("%s_%d",id,i++); print $0 }' | bgzip -c > annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 annotations.tsv.gz
            echo '##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
            ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,ID,REF,ALT,N_DISCOVERY_SAMPLES --output-type b in.bcf --output out.bcf
            df -h 1>&2
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            gcloud storage cp in.bcf ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf
            gcloud storage cp in.bcf.csi ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf.csi

            touch ~{chromosome}.done
            gcloud storage mv ~{chromosome}.done ~{remote_outdir}/~{chromosome}/
        fi
        echo "~{chromosome}" > out.txt
    >>>
    
    output {
        File out_txt = "out.txt"
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


# Performance on 12'680 samples, 15x, GRCh38, CAL_SENS<=0.999, HDD:
#
# TOOL                           CPU     RAM        TIME
# concat --naive truvari          1%     9KB        50m 
#
# SSD:
# concat truvari                200%    26KB      2h30m 
#
task AllChromosomes {
    input {
        Array[String] chromosomes
        Array[File] out_txt
        String remote_outdir
        Int naive
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 200
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
        
        
        # Localizing
        CHROMOSOMES=~{sep=',' chromosomes}
        echo ${CHROMOSOMES} | tr ',' '\n' > chr_list.txt
        rm -f file_list.txt
        while read -u 3 CHROMOSOME; do
            TEST=$( gcloud storage ls ~{remote_outdir}/${CHROMOSOME}/truvari_collapsed.bcf || echo 1 )
            if [ ${TEST} -eq 1 ]; then
                echo "ERROR: ${CHROMOSOME} has not been truvari collapsed."
                exit
            fi
            gcloud storage cp ~{remote_outdir}/${CHROMOSOME}/'*.bcf*' .
            mv truvari_collapsed.bcf ${CHROMOSOME}_truvari_collapsed.bcf
            mv truvari_collapsed.bcf.csi ${CHROMOSOME}_truvari_collapsed.bcf.csi
            echo ${CHROMOSOME}_truvari_collapsed.bcf >> file_list.txt
        done 3< chr_list.txt
        
        # Concatenating
        if [ ~{naive} -eq 1 ]; then
            CONCAT_FLAGS="--naive"
        else
            CONCAT_FLAGS=" "
        fi
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} ${CONCAT_FLAGS} --file-list file_list.txt --output-type b --output truvari_collapsed.bcf
        bcftools index --threads ${N_THREADS} -f truvari_collapsed.bcf
        
        # Uploading
        gcloud storage mv truvari_collapsed.'bcf*' ~{remote_outdir}/
    >>>

    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
