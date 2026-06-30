version 1.0

import "../call_structural_variants/call_structural_variants.wdl" as call_structural_variants
import "../annotate_structural_variants/annotate_structural_variants.wdl" as annotate_structural_variants
import "../select_structural_variants/select_structural_variants.wdl" as select_structural_variants

workflow call_select_svs {
    meta {
        description: "call structural variants with Manta, annotate with VEP and gene reannotation, and select somatic SVs"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        tumor_bam: "tumor BAM file"
        tumor_bai: "index of tumor_bam"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"

        # outputs
        sv_candidate_vcf: "bgzip-compressed VCF of unfiltered SV candidates"
        sv_candidate_vcf_index: "index of sv_candidate_vcf"
        sv_candidate_indel_vcf: "bgzip-compressed VCF of small indel candidates for use by SNV callers"
        sv_candidate_indel_vcf_index: "index of sv_candidate_indel_vcf"
        sv_somatic_vcf: "bgzip-compressed VCF of somatic SVs (or tumor SVs in tumor-only mode)"
        sv_somatic_vcf_index: "index of sv_somatic_vcf"
        sv_annot_vcf: "bgzip-compressed VCF with Ensembl VEP annotations added"
        sv_annot_bedpe: "BEDPE of annotated SVs converted from sv_annot_vcf"
        sv_annot_reannotated_bedpe: "BEDPE with genes reannotated by direct GTF intersection at each breakpoint"
        sv_del_annotation: "BED of genes overlapping the full span of deletion SVs"
        sv_dup_annotation: "BED of genes overlapping the full span of duplication SVs"
        sv_selected_somatic: "Parquet file of selected somatic structural variants"
    }

    input {
        String sample_id
        File tumor_bam
        File tumor_bai
        File ref_fasta
        File ref_fasta_index
    }

    call call_structural_variants.call_structural_variants {
        input:
            sample_id = sample_id,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    call annotate_structural_variants.annotate_structural_variants {
        input:
            sample_id = sample_id,
            vcf = call_structural_variants.sv_somatic_vcf
    }

    call select_structural_variants.select_structural_variants {
        input:
            sample_id = sample_id,
            input_bedpe = annotate_structural_variants.sv_annot_bedpe,
            gene_annotation = annotate_structural_variants.sv_annot_reannotated_bedpe,
            del_annotation = annotate_structural_variants.sv_del_annotation,
            dup_annotation = annotate_structural_variants.sv_dup_annotation
    }

    output {
        File sv_candidate_vcf = call_structural_variants.sv_candidate_vcf
        File sv_candidate_vcf_index = call_structural_variants.sv_candidate_vcf_index
        File sv_candidate_indel_vcf = call_structural_variants.sv_candidate_indel_vcf
        File sv_candidate_indel_vcf_index = call_structural_variants.sv_candidate_indel_vcf_index
        File sv_somatic_vcf = call_structural_variants.sv_somatic_vcf
        File sv_somatic_vcf_index = call_structural_variants.sv_somatic_vcf_index
        File sv_annot_vcf = annotate_structural_variants.sv_annot_vcf
        File sv_annot_bedpe = annotate_structural_variants.sv_annot_bedpe
        File sv_annot_reannotated_bedpe = annotate_structural_variants.sv_annot_reannotated_bedpe
        File sv_del_annotation = annotate_structural_variants.sv_del_annotation
        File sv_dup_annotation = annotate_structural_variants.sv_dup_annotation
        File sv_selected_somatic = select_structural_variants.sv_selected_somatic
    }
}
