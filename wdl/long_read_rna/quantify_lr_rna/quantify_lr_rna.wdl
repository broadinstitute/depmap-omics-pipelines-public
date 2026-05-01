version 1.0

workflow quantify_lr_rna {
    meta {
        description: "Quantify long-read RNA-seq expression with IsoQuant"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_bam: "coordinate-sorted long-read RNA-seq BAM"
        input_bai: "index of input_bam"
        ref_fasta: "reference FASTA"
        ref_annotation_db: "IsoQuant annotation database built from the reference GTF"
        data_type: "IsoQuant --data_type (e.g. 'nanopore' or 'pacbio_ccs')"
        model_construction_strategy: "IsoQuant --model_construction_strategy"
        stranded: "IsoQuant --stranded (e.g. 'forward', 'reverse', or 'none')"
        transcript_quantification: "IsoQuant --transcript_quantification method"
        gene_quantification: "IsoQuant --gene_quantification method"
        report_novel_unspliced: "IsoQuant --report_novel_unspliced setting"
        report_canonical: "IsoQuant --report_canonical setting"

        # outputs
        discovered_transcript_counts: "gzip-compressed TSV of read counts for IsoQuant-discovered transcripts"
        discovered_transcript_tpm: "gzip-compressed TSV of TPM values for IsoQuant-discovered transcripts"
        exon_counts: "gzip-compressed TSV of per-exon read counts"
        extended_annotation: "gzip-compressed GTF of the extended annotation including novel transcripts"
        gene_counts: "gzip-compressed TSV of per-gene read counts"
        gene_tpm: "gzip-compressed TSV of per-gene TPM values"
        intron_counts: "gzip-compressed TSV of per-intron read counts"
        read_assignments: "gzip-compressed TSV of per-read transcript assignments"
        transcript_counts: "gzip-compressed TSV of per-transcript read counts for reference transcripts"
        transcript_model_reads: "gzip-compressed TSV mapping transcript models to supporting reads"
        transcript_tpm: "gzip-compressed TSV of per-transcript TPM values for reference transcripts"
    }

    input {
        String sample_id
        File input_bam
        File input_bai
        File ref_fasta = "gs://ccleparams/hg38ref_no_alt/GRCh38_no_alt.fa"
        File ref_annotation_db = "gs://ccleparams/long-read/gencode/v38/gencode.v38.primary_assembly.annotation.db"

        # isoquant options
        String data_type
        String model_construction_strategy
        String stranded
        String transcript_quantification
        String gene_quantification
        String report_novel_unspliced
        String report_canonical
    }

    call run_isoquant {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            ref_annotation_db = ref_annotation_db,
            ref_fasta = ref_fasta,
            data_type = data_type,
            model_construction_strategy = model_construction_strategy,
            stranded = stranded,
            transcript_quantification = transcript_quantification,
            gene_quantification = gene_quantification,
            report_novel_unspliced = report_novel_unspliced,
            report_canonical = report_canonical
    }

    output {
        File discovered_transcript_counts = run_isoquant.discovered_transcript_counts
        File discovered_transcript_tpm = run_isoquant.discovered_transcript_tpm
        File exon_counts = run_isoquant.exon_counts
        File extended_annotation = run_isoquant.extended_annotation
        File gene_counts = run_isoquant.gene_counts
        File gene_tpm = run_isoquant.gene_tpm
        File intron_counts = run_isoquant.intron_counts
        File read_assignments = run_isoquant.read_assignments
        File transcript_counts = run_isoquant.transcript_counts
        File transcript_model_reads= run_isoquant.transcript_model_reads
        File transcript_tpm = run_isoquant.transcript_tpm
    }
}

task run_isoquant {
    meta {
        description: "Quantify long-read RNA-seq expression and discover novel transcripts using IsoQuant"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_bam: "coordinate-sorted long-read RNA-seq BAM"
        input_bai: "index of input_bam"
        ref_annotation_db: "IsoQuant annotation database built from the reference GTF"
        ref_fasta: "reference FASTA"
        data_type: "IsoQuant --data_type (e.g. 'nanopore' or 'pacbio_ccs')"
        model_construction_strategy: "IsoQuant --model_construction_strategy"
        stranded: "IsoQuant --stranded (e.g. 'forward', 'reverse', or 'none')"
        transcript_quantification: "IsoQuant --transcript_quantification method"
        gene_quantification: "IsoQuant --gene_quantification method"
        report_novel_unspliced: "IsoQuant --report_novel_unspliced setting"
        report_canonical: "IsoQuant --report_canonical setting"

        # outputs
        discovered_transcript_counts: "gzip-compressed TSV of read counts for IsoQuant-discovered transcripts"
        discovered_transcript_tpm: "gzip-compressed TSV of TPM values for IsoQuant-discovered transcripts"
        exon_counts: "gzip-compressed TSV of per-exon read counts"
        extended_annotation: "gzip-compressed GTF of the extended annotation including novel transcripts"
        gene_counts: "gzip-compressed TSV of per-gene read counts"
        gene_tpm: "gzip-compressed TSV of per-gene TPM values"
        intron_counts: "gzip-compressed TSV of per-intron read counts"
        read_assignments: "gzip-compressed TSV of per-read transcript assignments"
        transcript_counts: "gzip-compressed TSV of per-transcript read counts for reference transcripts"
        transcript_model_reads: "gzip-compressed TSV mapping transcript models to supporting reads"
        transcript_tpm: "gzip-compressed TSV of per-transcript TPM values for reference transcripts"
    }

    input {
        String sample_id
        File input_bam
        File input_bai
        File ref_annotation_db
        File ref_fasta
        String data_type
        String model_construction_strategy
        String stranded
        String transcript_quantification
        String gene_quantification
        String report_novel_unspliced
        String report_canonical

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/isoquant"
        String docker_image_hash_or_tag = ":3.7.0--hdfd78af_0"
        Int cpu = 8
        Int mem_gb = 32
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(
            size(input_bam, "GiB") + size(ref_annotation_db, "GiB")
            + size(ref_fasta, "GiB")
        ) + 20 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        # Ensure BAI is newer than BAM to avoid warnings
        touch "~{input_bai}"

        # Run IsoQuant with custom annotation database
        /usr/local/bin/isoquant.py \
            --bam "~{input_bam}" \
            --count_exons \
            --data_type "~{data_type}" \
            --gene_quantification "~{gene_quantification}" \
            --genedb "~{ref_annotation_db}" \
            --labels "~{sample_id}" \
            --model_construction_strategy "~{model_construction_strategy}" \
            --prefix "~{sample_id}" \
            --reference "~{ref_fasta}" \
            --report_canonical "~{report_canonical}" \
            --report_novel_unspliced "~{report_novel_unspliced}" \
            --stranded "~{stranded}" \
            --threads ~{cpu} \
            --transcript_quantification "~{transcript_quantification}"

        echo "Compressing output TSVS"
        find isoquant_output/~{sample_id}/ -maxdepth 1 -type f -not -name '*.gz' -exec gzip {} +
    >>>

    output {
        File discovered_transcript_counts = "isoquant_output/~{sample_id}/~{sample_id}.discovered_transcript_counts.tsv.gz"
        File discovered_transcript_tpm = "isoquant_output/~{sample_id}/~{sample_id}.discovered_transcript_tpm.tsv.gz"
        File exon_counts = "isoquant_output/~{sample_id}/~{sample_id}.exon_counts.tsv.gz"
        File extended_annotation = "isoquant_output/~{sample_id}/~{sample_id}.extended_annotation.gtf.gz"
        File gene_counts = "isoquant_output/~{sample_id}/~{sample_id}.gene_counts.tsv.gz"
        File gene_tpm = "isoquant_output/~{sample_id}/~{sample_id}.gene_tpm.tsv.gz"
        File intron_counts = "isoquant_output/~{sample_id}/~{sample_id}.intron_counts.tsv.gz"
        File read_assignments = "isoquant_output/~{sample_id}/~{sample_id}.read_assignments.tsv.gz"
        File transcript_counts = "isoquant_output/~{sample_id}/~{sample_id}.transcript_counts.tsv.gz"
        File transcript_model_reads = "isoquant_output/~{sample_id}/~{sample_id}.transcript_model_reads.tsv.gz"
        File transcript_tpm = "isoquant_output/~{sample_id}/~{sample_id}.transcript_tpm.tsv.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
