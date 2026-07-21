version 1.0


# Essentially identical to `Workpackage8.wdl`, but simplified for ultralong and
# BND VCFs. Concatenates the `truvari collapse` chunks. Ensures that every
# record in output has a globally unique ID (duplicated IDs may arise naturally
# from the previous steps of the pipeline), and an INFO field that counts the
# number of samples it was discovered in.
#
workflow SV_Integration_Workpackage15 {
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


# Performance on 12'680 samples, 15x, GRCh38, chr1, HDD:
#
# TOOL                           CPU     RAM   TIME
# bcftools concat                
# bcftools query                
# bcftools annotate             
#
task SingleChromosome {
    input {
        String chromosome
        String remote_indir
        String remote_outdir

        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 2
        Int disk_size_gb = 10
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
            N_CHUNKS=$(wc -l < chunk_list.txt)
            if [ ${N_CHUNKS} -gt 1 ]; then
                ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type b --output out.bcf
                df -h 1>&2
                rm -rf chunk_* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            else
                CHUNK_FILE=$(head -n 1 chunk_list.txt)
                mv ${CHUNK_FILE} in.bcf
                mv ${CHUNK_FILE}.csi in.bcf.csi
            fi
            
            # Enforcing a distinct ID in every record, and annotating every
            # record with the number of samples it occurs in. Note that the
            # latter is not equal to the QUAL field in input to truvari collapse
            # upstream, so we have to recompute this number.
            CHR=~{chromosome}
            CHR=${CHR#chr}
            # Enforcing a distinct ID in every record, and annotating every record
            # with the number of samples it occurs in. Note that the latter is not
            # equal to the QUAL field in input to truvari collapse upstream, so we
            # have to recompute this number.
            #
            # Both operations are keyed on record ORDER, not on CHROM/POS/REF/ALT.
            # A position-based `bcftools annotate` matches on CHROM,POS,REF,ALT and
            # cannot distinguish same-position symbolic records (identical
            # CHROM/POS/REF/ALT but different END/SVLEN); it would assign such
            # records the same ID and copy the same N_DISCOVERY_SAMPLES to all of
            # them. Matching on ID (`-c ~ID`) is not usable either (unimplemented
            # in bcftools). `bcftools query` and `bcftools view` both emit records
            # in file order, so counts.txt aligns 1:1 with the streamed records.
            ${TIME_COMMAND} bcftools query --format '%COUNT(GT="alt")\n' in.bcf > counts.txt
            ${TIME_COMMAND} bcftools view in.bcf | awk -v id=${CHR} 'BEGIN { FS="\t"; OFS="\t"; i=0; j=0; while ((getline c < "counts.txt") > 0) { cnt[j++]=c; } } /^#CHROM/ { print "##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description=\"Number of samples where the record was discovered\">"; print $0; next } /^#/ { print $0; next } { $3=sprintf("%s_%d",id,i); if ($8==".") { $8="N_DISCOVERY_SAMPLES=" cnt[i] } else { $8=$8 ";N_DISCOVERY_SAMPLES=" cnt[i] } i++; print $0 }' | bcftools view --output-type b --output out.bcf
            df -h 1>&2
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            gcloud storage cp in.bcf ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf
            gcloud storage cp in.bcf.csi ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf.csi
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


# Performance on 12'680 samples, 15x, GRCh38, HDD:
#
# TOOL                           CPU     RAM        TIME
# concat --naive truvari         
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
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f truvari_collapsed.bcf
        
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
