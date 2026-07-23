version 1.0

import "../prep_annotations/prep_annotations.wdl" as prep_annotations
import "../annotate_mutations/annotate_mutations.wdl" as annotate_mutations
import "../annotate_mutations_merge/annotate_mutations_merge.wdl" as annotate_mutations_merge
import "../select_somatic_variants/select_somatic_variants.wdl" as select_somatic_variants

workflow prep_annotate_select_muts {
    meta {
        description: "prepare a VCF for annotation, annotate with multiple sources, merge annotations, and select somatic variants"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        input_vcf: "input somatic variant VCF to prepare and annotate"
        vcf_patch: "optional VCF of curated variants to substitute into patch_region_bed regions"

        # outputs
        mut_annot_bcftools_fixed_vcf: "bgzip-compressed VCF after all fixes and filters have been applied"
        mut_annot_bcftools_vcf: "bgzip-compressed VCF with bcftools-based annotations added"
        mut_annot_vep_vcf: "bgzip-compressed VCF with Ensembl VEP annotations added"
        mut_annot_vcf: "bgzip-compressed VCF with INFO fields merged from all annotated VCFs"
        mut_annot_vcf_index: "index of mut_annot_vcf"
        mut_duckdb: "SQL schema and data Parquet files produced by vcf-to-duckdb for loading into DuckDB"
        mut_somatic_variants: "Parquet file of selected somatic variants"
        mut_sig_variants: "Parquet file of unfiltered variants (for mutational signature analysis)"
    }

    input {
        String sample_id
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_vcf
        File? vcf_patch
    }

    call prep_annotations.prep_annotations {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_vcf = input_vcf,
            vcf_patch = vcf_patch
    }

    call annotate_mutations.annotate_mutations {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_vcf = prep_annotations.mut_annot_bcftools_fixed_vcf
    }

    call annotate_mutations_merge.annotate_mutations_merge {
        input:
            sample_id = sample_id,
            input_vcfs = select_all([
                annotate_mutations.mut_annot_bcftools_vcf,
                annotate_mutations.mut_annot_vep_vcf
            ])
    }

    call select_somatic_variants.select_somatic_variants {
        input:
            sample_id = sample_id,
            duckdb = annotate_mutations_merge.mut_duckdb
    }

    output {
        File mut_annot_bcftools_fixed_vcf = prep_annotations.mut_annot_bcftools_fixed_vcf
        File? mut_annot_bcftools_vcf = annotate_mutations.mut_annot_bcftools_vcf
        File? mut_annot_vep_vcf = annotate_mutations.mut_annot_vep_vcf
        File mut_annot_vcf = annotate_mutations_merge.mut_annot_vcf
        File mut_annot_vcf_index = annotate_mutations_merge.mut_annot_vcf_index
        Array[File] mut_duckdb = annotate_mutations_merge.mut_duckdb
        File mut_somatic_variants = select_somatic_variants.mut_somatic_variants
        File mut_sig_variants = select_somatic_variants.mut_sig_variants
    }
}
