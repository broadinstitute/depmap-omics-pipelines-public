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
        output_cram_or_bam: "output file format; either 'BAM' or 'CRAM'"

        # outputs
        patched_bam_cram: "patched BAM or CRAM; either the full file with realigned reads merged in, or just the realigned reads from the problematic regions"
        patched_bai_crai: "index of patched_bam_cram"
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
        String output_cram_or_bam
    }

    call do_patch_hg38_bam {
        input:
            sample_id = sample_id,
            cram_bam = select_first([bam, cram]),
            crai_bai = select_first([bai, crai]),
            problems_bed = problems_bed,
            old_ref_fasta = old_ref_fasta,
            old_ref_fasta_index = old_ref_fasta_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            output_cram_or_bam = output_cram_or_bam
    }

    output {
        File patched_bam_cram = do_patch_hg38_bam.patched_bam_cram
        File patched_bai_crai = do_patch_hg38_bam.patched_bai_crai
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
        output_cram_or_bam: "output file format; either 'BAM' or 'CRAM'"

        # outputs
        patched_bam_cram: "patched BAM or CRAM; either the full file with realigned reads merged in, or just the realigned reads from the problematic regions"
        patched_bai_crai: "index of patched_bam_cram"
    }

    input {
        String sample_id
        File cram_bam
        File crai_bai
        File problems_bed
        File old_ref_fasta
        File old_ref_fasta_index
        File ref_fasta
        File ref_fasta_index
        String output_cram_or_bam

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools_picard"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            10 * size(cram_bam, "GiB")
            + size(ref_fasta, "GiB")
        ) + additional_disk_gb
    )

    Int n_threads = cpu

    parameter_meta {
        cram_bam: { localization_optional: true }
        crai_bai: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        # Get problematic subset of full input BAM
        samtools view \
            -@ ~{n_threads} \
            -M \
            -h \
            -T "~{old_ref_fasta}" \
            -L "~{problems_bed}" \
            "~{cram_bam}" \
            -o "subset.bam"

        # Get read names overlapping problem regions
        samtools view \
            -@ ~{n_threads} \
            "subset.bam" \
            | cut -f1 \
            | LC_ALL=C sort -u \
            > "problem_read_names.txt"

        # Extract full alignments of those reads
        samtools view \
            -@ ~{n_threads} \
            -N "problem_read_names.txt" \
            -o "problem_reads.bam" \
            "subset.bam"

        # Convert to FASTQ (preserves qualities)
        samtools fastq \
            -@ ~{n_threads} \
            -1 "reads_1.fq.gz" \
            -2 "reads_2.fq.gz" \
            -s "singletons.fq.gz" \
            -0 /dev/null \
            -n \
            "problem_reads.bam"

        # Realign paired-end reads to corrected reference
        bwa mem -t ~{n_threads} "~{ref_fasta}" "reads_1.fq.gz" "reads_2.fq.gz" \
          | samtools view -b -@ ~{n_threads} - \
          | samtools sort -@ ~{n_threads} -o "realigned.pe.bam"

        if [ -s "singletons.fq.gz" ]; then
            bwa mem -t ~{n_threads} "~{ref_fasta}" "singletons.fq.gz" \
              | samtools sort -@ ~{n_threads} -o "realigned.se.bam" -
            samtools merge -@ ~{n_threads} -f "realigned.bam" "realigned.pe.bam" "realigned.se.bam"
        else
            cp "realigned.pe.bam" "realigned.bam"
        fi

        # Remove problem reads from subsetted BAM/CRAM
        samtools view \
            -@ ~{n_threads} \
            -T "~{old_ref_fasta}" \
            -b \
            -N "^problem_read_names.txt" \
            -o "subset_without_problem_reads.bam" \
            "subset.bam"

        # Merge realigned reads back and sort
        samtools merge -@ ~{n_threads} -u - "subset_without_problem_reads.bam" "realigned.bam" \
          | samtools sort -@ ~{n_threads} -o "patched.bam" -

        # Add read group tags (required by GATK tools)
        java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups \
            I="patched.bam" \
            O="~{sample_id}.patch.bam" \
            RGID=1 \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=unit1 \
            RGSM="~{sample_id}"

        if [ "~{output_cram_or_bam}" = "CRAM" ]; then
            samtools view \
                -@ ~{n_threads} \
                -C \
                -T "~{ref_fasta}" \
                -o "~{sample_id}.patch.cram" \
                "~{sample_id}.patch.bam"
            samtools index -@ ~{n_threads} "~{sample_id}.patch.cram"
            mv "~{sample_id}.patch.cram.crai" "~{sample_id}.patch.crai"
        else
            samtools index -@ ~{n_threads} "~{sample_id}.patch.bam"
            mv "~{sample_id}.patch.bam.bai" "~{sample_id}.patch.bai"
        fi
    >>>

    output {
        File patched_bam_cram = if output_cram_or_bam == "CRAM"
            then "~{sample_id}.patch.cram"
            else "~{sample_id}.patch.bam"
        File patched_bai_crai = if output_cram_or_bam == "CRAM"
            then "~{sample_id}.patch.crai"
            else "~{sample_id}.patch.bai"
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
