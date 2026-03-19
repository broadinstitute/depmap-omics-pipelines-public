version 1.0

import "../align_sr_rna/align_sr_rna.wdl" as align_sr_rna
import "../call_sr_rna_fusions/call_sr_rna_fusions.wdl" as call_sr_rna_fusions
import "../quantify_sr_rna/quantify_sr_rna.wdl" as quantify_sr_rna

workflow process_sr_rna {
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

    meta {
        allowNestedInputs: true
    }
}
