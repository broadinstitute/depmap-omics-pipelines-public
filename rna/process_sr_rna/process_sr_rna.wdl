version 1.0

import "../align_sr_rna/align_sr_rna.wdl" as align_sr_rna
import "../call_sr_rna_fusions/call_sr_rna_fusions.wdl" as call_sr_rna_fusions
import "../quantify_sr_rna/quantify_sr_rna.wdl" as quantify_sr_rna

workflow process_sr_rna {
    meta {
        description: "End-to-end short-read RNA-seq pipeline: align with STAR,call fusions with Arriba, andquantify expression with Salmon"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_file_type: "'FASTQ', 'BAM', or 'CRAM'; determines how input reads are provided to STAR"
        fastqs: "input FASTQ files; required when input_file_type is FASTQ"
        cram_bam: "input BAM or CRAM file; required when input_file_type is BAM or CRAM"
        crai_bai: "index of cram_bam"
        align_ref_fasta: "reference FASTA used to decode CRAM input for alignment"
        align_ref_fasta_index: "index of align_ref_fasta"
        ref_fasta: "reference FASTA passed to Arriba for fusion calling"
        star_index: "tar.gz archive of the STAR genome index"
        gtf: "GTF annotation file used by Arriba and Salmon"
        targets: "FASTA of transcript sequences used as Salmon quantification targets"
        stranded: "if true, Salmon also runs with automatic library type detection"
        call_fusions: "if true, call RNA fusions with Arriba"
        quantify: "if true, quantify expression with Salmon"

        # outputs
        reads_per_gene: "gzip-compressed TSV of per-gene read counts from STAR"
        junctions: "gzip-compressed TSV of splice junctions detected by STAR"
        fusions: "TSV of fusion gene candidates from Arriba"
        fusions_discarded: "TSV of discarded fusion candidates from Arriba"
        quant_transcripts_auto: "Salmon transcript-level quantification TSV with automatic library type"
        quant_genes_auto: "Salmon gene-level quantification TSV with automatic library type"
        quant_transcripts_iu: "Salmon transcript-level quantification TSV using unstranded IU library type"
        quant_genes_iu: "Salmon gene-level quantification TSV using unstranded IU library type"
    }

    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File? align_ref_fasta
        File? align_ref_fasta_index
        File ref_fasta
        File star_index
        File gtf
        File targets
        Boolean stranded

        Boolean call_fusions = true
        Boolean quantify = true
    }

    call align_sr_rna.align_sr_rna {
        input:
            sample_id = sample_id,
            input_file_type = input_file_type,
            fastqs = fastqs,
            cram_bam = cram_bam,
            crai_bai = crai_bai,
            ref_fasta = align_ref_fasta,
            ref_fasta_index = align_ref_fasta_index,
            star_index = star_index
    }

    if (call_fusions) {
        call call_sr_rna_fusions.call_sr_rna_fusions {
            input:
                sample_id = sample_id,
                analysis_ready_bam = align_sr_rna.analysis_ready_bam,
                ref_fasta = ref_fasta,
                gtf = gtf
        }
    }

    if (quantify) {
        call quantify_sr_rna.quantify_sr_rna {
            input:
                sample_id = sample_id,
                transcriptome_bam = align_sr_rna.transcriptome_bam,
                gtf = gtf,
                targets = targets,
                stranded = stranded
        }
    }

    output {
        File reads_per_gene = align_sr_rna.reads_per_gene
        File junctions = align_sr_rna.junctions
        File? fusions = call_sr_rna_fusions.fusions
        File? fusions_discarded = call_sr_rna_fusions.fusions_discarded
        File? quant_transcripts_auto = quantify_sr_rna.quant_transcripts_auto
        File? quant_genes_auto = quantify_sr_rna.quant_genes_auto
        File? quant_transcripts_iu = quantify_sr_rna.quant_transcripts_iu
        File? quant_genes_iu = quantify_sr_rna.quant_genes_iu
    }
}
