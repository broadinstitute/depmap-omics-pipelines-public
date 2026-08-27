version 1.0

import "../../wgs/call_mutations/call_mutations.wdl" as call_mutations_wf
import "../../wgs/prep_annotations/prep_annotations.wdl" as prep_annotations_wf
import "../../wgs/annotate_mutations/annotate_mutations.wdl" as annotate_mutations_wf
import "../../wgs/annotate_mutations_merge/annotate_mutations_merge.wdl" as annotate_mutations_merge_wf
import "../../wgs/select_somatic_variants/select_somatic_variants.wdl" as select_somatic_variants_wf

workflow call_select_muts {
    meta {
        description: "call somatic variants with Mutect2, annotate, and select somatic variants"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        cram_bam: "BAM/CRAM file of tumor reads"
        crai_bai: "index of cram_bam"
        normal_reads: "optional BAM/CRAM file of matched normal reads"
        normal_reads_index: "index of normal_reads"
        scatter_count: "number of shards for Mutect2"

        # outputs
        mut_vcf: "bgzip-compressed VCF of Mutect2 filtered variants"
        mut_vcf_index: "index of mut_vcf"
        mut_annot_vcf: "bgzip-compressed VCF of annotated variants"
        mut_annot_vcf_index: "index of mut_annot_vcf"
        mut_somatic_variants: "Parquet file of selected somatic variants"
        mut_passed_variants: "Parquet file of all variants that pass quality filters"
    }

    input {
        String sample_id
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File cram_bam
        File crai_bai

        # call_mutations
        Int scatter_count

        # annotate_mutations
        String vep_pick_order
    }

    call call_mutations_wf.call_mutations {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_index,
            ref_dict = ref_dict,
            tumor_reads = cram_bam,
            tumor_reads_index = crai_bai,
            scatter_count = scatter_count
    }

    call prep_annotations_wf.prep_annotations {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_vcf = call_mutations.filtered_vcf
    }

    call annotate_mutations_wf.annotate_mutations {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_vcf = prep_annotations.mut_annot_bcftools_fixed_vcf,
            vep_pick_order = vep_pick_order
    }

    call annotate_mutations_merge_wf.annotate_mutations_merge {
        input:
            sample_id = sample_id,
            input_vcfs = select_all([
                annotate_mutations.mut_annot_bcftools_vcf,
                annotate_mutations.mut_annot_vep_vcf
            ])
    }

    call select_somatic_variants_wf.select_somatic_variants {
        input:
            sample_id = sample_id,
            duckdb = annotate_mutations_merge.mut_duckdb
    }

    output {
        File mut_vcf = call_mutations.filtered_vcf
        File mut_vcf_index = call_mutations.filtered_vcf_idx
        File mut_annot_vcf = annotate_mutations_merge.mut_annot_vcf
        File mut_annot_vcf_index = annotate_mutations_merge.mut_annot_vcf_index
        File mut_somatic_variants = select_somatic_variants.mut_somatic_variants
        File mut_passed_variants = select_somatic_variants.mut_sig_variants
    }
}
