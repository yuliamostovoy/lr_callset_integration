version 1.0


# Runs `truvari collapse` on a small chunk of a bcftools merged cohort VCF.
# Default truvari parameters are optimal for 10k samples.
#
workflow SV_Integration_Workpackage7 {
    input {
        String remote_indir
        String chromosome_id
        File chunks_ids
        String remote_outdir
        
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Boolean use_bed
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        truvari_matching_parameters: "Truvari's definition of a match. The default string contains truvari's defaults."
    }
    
    call Impl {
        input:
            remote_indir = remote_indir,
            chromosome_id = chromosome_id,
            chunks_ids = chunks_ids,
            remote_outdir = remote_outdir,
            
            truvari_matching_parameters = truvari_matching_parameters,
            use_bed = use_bed,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, one 5 MB chunk of chr6,
# CAL_SENS<=0.999:
#
# TOOL               CPU     RAM     TIME   OUTPUT VCF
# bcftools query     100%    600M      10s  
# bcftools annotate  200%    600M      20s
# truvari collapse   100%      8G    3-30m  500M
# bcftools sort      100%      2G       1m
#
# Peak disk (old measurement with `removed.vcf`): 16G
#
task Impl {
    input {
        String remote_indir
        String chromosome_id
        File chunks_ids
        String remote_outdir
        
        String truvari_matching_parameters
        Boolean use_bed
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 16
        Int disk_size_gb = 20
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # Sets QUAL to the number of samples where a record was discovered, to
        # simulate `--keep common` (which is slow on 10k samples) with `--keep
        # maxqual` in `truvari collapse`. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # Remark: we use the number of samples rather than AC, since we don't
        # trust genotypes at this stage, and since a record being discovered
        # independently in more samples is more informative than it being
        # genotyped more times in fewer samples.
        #
        function CopyNSamplesToQual() {
            local CHUNK_ID=$1
        
            mv chunk_${CHUNK_ID}.bcf chunk_${CHUNK_ID}_in.bcf
            mv chunk_${CHUNK_ID}.bcf.csi chunk_${CHUNK_ID}_in.bcf.csi
        
            # Remark: we cannot join annotations just by ID at this stage,
            # since the IDs in the output of bcftools merge are not necessarily
            # all distinct.
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' chunk_${CHUNK_ID}_in.bcf | bgzip -c > chunk_${CHUNK_ID}_annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 chunk_${CHUNK_ID}_annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations chunk_${CHUNK_ID}_annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,QUAL --output-type z chunk_${CHUNK_ID}_in.bcf --output chunk_${CHUNK_ID}_out.vcf.gz
            rm -f chunk_${CHUNK_ID}_in.bcf* ; mv chunk_${CHUNK_ID}_out.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t chunk_${CHUNK_ID}_in.vcf.gz
            rm -f chunk_${CHUNK_ID}_annotations.tsv.gz
            
            mv chunk_${CHUNK_ID}_in.vcf.gz chunk_${CHUNK_ID}_annotated.vcf.gz
            mv chunk_${CHUNK_ID}_in.vcf.gz.tbi chunk_${CHUNK_ID}_annotated.vcf.gz.tbi
        }
        
        
        # Remark: in theory we should set `--gt all` to avoid collapsing records
        # that are present in the same sample, since we assume that intra-
        # sample merging has already done that upstream. In practice `--gt all`
        # is too slow on 10k samples.
        #
        # Remark: to further improve speed we could think of dropping genotypes
        # before running truvari collapse. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # However, this would also discard e.g. SUPP fields that were copied to
        # FORMAT upstream, so it is not correct for our setup. It would also
        # make it impossible e.g. to compare precision/recall after collapse to
        # precision/recall after cohort re-genotyping.
        #
        function Collapse() {
            local CHUNK_ID=$1
            
            mv chunk_${CHUNK_ID}_annotated.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz
            mv chunk_${CHUNK_ID}_annotated.vcf.gz.tbi chunk_${CHUNK_ID}_in.vcf.gz.tbi
            
            # Remark: we do not store `removed.vcf` since it's several GBs per
            # chunk and much bigger than the collapsed output.
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax ${INFINITY} --keep maxqual --gt off ~{truvari_matching_parameters} ${BED_FLAGS} --input chunk_${CHUNK_ID}_in.vcf.gz --output chunk_${CHUNK_ID}_out.vcf --removed-output /dev/null
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf chunk_${CHUNK_ID}_in.vcf
        
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type b chunk_${CHUNK_ID}_in.vcf --output chunk_${CHUNK_ID}_out.bcf
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf
                
            # Setting to `.` every ID written by truvari collapse, since these
            # can be very long in a large cohort and needlessly inflate output
            # size. Unique IDs genome-wide will be assigned downstream.
            ${TIME_COMMAND} bcftools annotate --remove ID --output-type b chunk_${CHUNK_ID}_in.bcf --output chunk_${CHUNK_ID}_out.bcf
            rm -f chunk_${CHUNK_ID}_in.bcf* ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf
                
            mv chunk_${CHUNK_ID}_in.bcf chunk_${CHUNK_ID}_truvari.bcf
            mv chunk_${CHUNK_ID}_in.bcf.csi chunk_${CHUNK_ID}_truvari.bcf.csi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        df -h 1>&2
        
        if ~{use_bed} ; then
            gcloud storage cp ~{remote_indir}/~{chromosome_id}/included.bed .
            BED_FLAGS="--bed included.bed"
        else 
            BED_FLAGS=" "
        fi
        while read -u 3 CHUNK_ID; do
            # Skipping the chunk if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/~{chromosome_id}/chunk_${CHUNK_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Collapsing
            gcloud storage cp ~{remote_indir}/~{chromosome_id}/chunk_${CHUNK_ID}.'bcf*' .
            CopyNSamplesToQual ${CHUNK_ID}
            Collapse ${CHUNK_ID}
            
            # Uploading
            gcloud storage mv chunk_${CHUNK_ID}_truvari.bcf'*' ~{remote_outdir}/~{chromosome_id}/
            touch chunk_${CHUNK_ID}.done
            gcloud storage mv chunk_${CHUNK_ID}.done ~{remote_outdir}/~{chromosome_id}/
            rm -rf chunk_${CHUNK_ID}*
            ls -laht 1>&2
        done 3< ~{chunks_ids}
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
