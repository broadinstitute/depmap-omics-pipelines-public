version 1.0

workflow patch_hg38_bam {
    meta {
        description: "Patch a BAM or CRAM file by realigning problematic genomic regions"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        bam: "input BAM file (provide bam+bai or cram+crai)"
        bai: "index of input BAM file"
        cram: "input CRAM file (provide bam+bai or cram+crai)"
        crai: "index of input CRAM file"
        problems_bed: "BED file of problematic regions whose overlapping reads will be realigned"
        old_ref_fasta: "reference FASTA that the input BAM/CRAM was aligned to; required to decode CRAM input"
        old_ref_fasta_index: "index of old_ref_fasta"
        ref_fasta: "corrected reference FASTA to realign problem regions against"
        ref_fasta_index: "index of ref_fasta"
        ref_fasta_amb: "BWA supporting file for ref_fasta"
        ref_fasta_ann: "BWA supporting file for ref_fasta"
        ref_fasta_bwt: "BWA supporting file for ref_fasta"
        ref_fasta_pac: "BWA supporting file for ref_fasta"
        ref_fasta_sa: "BWA supporting file for ref_fasta"

        # outputs
        patched_bam: "BAM containing only the realigned reads from the problematic regions"
        patched_bai: "index of patched_bam"
    }

    input {
        String sample_id
        File? bam
        File? bai
        File? cram
        File? crai
        File problems_bed
        File? old_ref_fasta
        File? old_ref_fasta_index
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa
    }

    call do_patch_hg38_bam {
        input:
            sample_id = sample_id,
            cram_bam = select_first([cram, bam]),
            crai_bai = select_first([crai, bai]),
            problems_bed = problems_bed,
            old_ref_fasta = old_ref_fasta,
            old_ref_fasta_index = old_ref_fasta_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_fasta_amb = ref_fasta_amb,
            ref_fasta_ann = ref_fasta_ann,
            ref_fasta_bwt = ref_fasta_bwt,
            ref_fasta_pac = ref_fasta_pac,
            ref_fasta_sa = ref_fasta_sa
    }

    output {
        File patched_bam = do_patch_hg38_bam.patched_bam
        File patched_bai = do_patch_hg38_bam.patched_bai
    }
}

task do_patch_hg38_bam {
    meta {
        description: "Patch a BAM or CRAM file by realigning problematic genomic regions"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        cram_bam: "input BAM or CRAM file"
        crai_bai: "index of cram_bam"
        problems_bed: "BED file of problematic regions whose overlapping reads will be realigned"
        old_ref_fasta: "reference FASTA that cram_bam was aligned to; required to decode CRAM input"
        old_ref_fasta_index: "index of old_ref_fasta"
        ref_fasta: "corrected reference FASTA to realign problem regions against"
        ref_fasta_index: "index of ref_fasta"
        ref_fasta_amb: "BWA supporting file for ref_fasta"
        ref_fasta_ann: "BWA supporting file for ref_fasta"
        ref_fasta_bwt: "BWA supporting file for ref_fasta"
        ref_fasta_pac: "BWA supporting file for ref_fasta"
        ref_fasta_sa: "BWA supporting file for ref_fasta"

        # outputs
        patched_bam: "BAM containing only the realigned reads from the problematic regions"
        patched_bai: "index of patched_bam"

        cram_bam: { localization_optional: true }
        crai_bai: { localization_optional: true }
        old_ref_fasta: { localization_optional: true }
        old_ref_fasta_index: { localization_optional: true }
    }

    input {
        String sample_id
        File cram_bam
        File crai_bai
        File problems_bed
        File? old_ref_fasta
        File? old_ref_fasta_index
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools_picard"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 8
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 100 + additional_disk_gb

    command <<<
        set -euo pipefail

        # set up auth for samtools streaming from GCS buckets
        export GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token)"
        export GCS_REQUESTER_PAYS_PROJECT="$(gcloud config get-value project -q)"

        echo "Subsetting input BAM"
        samtools view \
            -@ ~{cpu} \
            -M \
            -h \
            --customized-index \
            ~{'-T "' + old_ref_fasta + '"'} \
            -L "~{problems_bed}" \
            "~{cram_bam}" \
            "~{crai_bai}" \
            -o "subset.bam"

        echo "Getting read names overlapping problem regions"
        samtools view \
            -@ ~{cpu} \
            "subset.bam" \
            | cut -f1 \
            | LC_ALL=C sort -u \
            > "problem_read_names.txt"

        echo "Extracting alignments of those reads"
        samtools view \
            -@ ~{cpu} \
            -N "problem_read_names.txt" \
            -o "problem_reads.bam" \
            "subset.bam"

        echo "Converting to FASTQ (preserves qualities)"
        samtools fastq \
            -@ ~{cpu} \
            -1 "reads_1.fq.gz" \
            -2 "reads_2.fq.gz" \
            -s "singletons.fq.gz" \
            -0 /dev/null \
            -n \
            "problem_reads.bam"

        echo "Realigning paired-end reads to corrected reference"
        bwa mem -t ~{cpu} "~{ref_fasta}" "reads_1.fq.gz" "reads_2.fq.gz" -o "pe.sam"
        samtools view -b -@ ~{cpu} "pe.sam" -o "pe.bam"
        samtools sort -@ ~{cpu} "pe.bam" -o "realigned.pe.bam"
        rm "pe.sam" "pe.bam"

        if [ "$(gzip -dc "singletons.fq.gz" | wc -c)" -gt 0 ]; then
            echo "Realigning singleton reads to corrected reference"
            bwa mem -t ~{cpu} "~{ref_fasta}" "singletons.fq.gz" -o "se.sam"
            samtools view -b -@ ~{cpu} "se.sam" -o "se.bam"
            samtools sort -@ ~{cpu} "se.bam" -o "realigned.se.bam"
            rm "se.sam" "se.bam"
            samtools merge -@ ~{cpu} -f "realigned.bam" "realigned.pe.bam" "realigned.se.bam"
        else
            mv "realigned.pe.bam" "realigned.bam"
        fi

        echo "Adding arbitrary read group tags (required by GATK tools)"
        java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups \
            I="realigned.bam" \
            O="~{sample_id}.patch.bam" \
            RGID=1 \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=unit1 \
            RGSM="~{sample_id}"

        samtools index -@ ~{cpu} "~{sample_id}.patch.bam"
        mv "~{sample_id}.patch.bam.bai" "~{sample_id}.patch.bai"
    >>>

    output {
        File patched_bam = "~{sample_id}.patch.bam"
        File patched_bai = "~{sample_id}.patch.bai"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
