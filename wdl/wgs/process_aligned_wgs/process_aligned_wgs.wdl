version 1.0

import "../call_cnvs/call_cnvs.wdl" as call_cnvs
import "../call_mutations/call_mutations.wdl" as call_mutations
import "../call_select_svs/call_select_svs.wdl" as call_select_svs
import "../compute_molecular_signatures/compute_molecular_signatures.wdl" as compute_molecular_signatures
import "../infer_msi_status/infer_msi_status.wdl" as infer_msi_status
import "../make_guide_mutation_beds/make_guide_mutation_beds.wdl" as make_guide_mutation_beds
import "../patch_hg38_bam_call_muts/patch_hg38_bam_call_muts.wdl" as patch_hg38_bam_call_muts
import "../prep_annotate_select_muts/prep_annotate_select_muts.wdl" as prep_annotate_select_muts

workflow process_aligned_wgs {
    meta {
        description: "process an aligned WGS BAM"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        bam: "analysis-ready BAM file"
        bai: "index of bam"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"

        # outputs
        cnv_bin_coverage: "gzip-compressed TSV of raw read counts per genomic bin"
        cnv_segments: "copy number segments"
        cnv_read_cov_bin: "read coverage per bin used for segmentation"
        cnv_input_params: "parameters used for CNV segmentation"
        cnv_cn_by_gene_weighted_mean: "weighted mean copy number by gene"
        cnv_segments_plot: "plot of copy number segments"
        mut_vcf: "bgzip-compressed filtered somatic variant VCF from Mutect2"
        mut_vcf_idx: "index of mut_vcf"
        mut_filtering_stats: "Mutect2 filtering stats"
        patched_bam: "BAM with problematic regions realigned"
        patched_bai: "index of patched_bam"
        mut_vcf_patch: "bgzip-compressed filtered somatic variant VCF called from the patched BAM"
        mut_vcf_patch_idx: "index of mut_vcf_patch"
        mut_annot_bcftools_fixed_vcf: "bgzip-compressed VCF after all fixes and filters have been applied"
        mut_annot_bcftools_vcf: "bgzip-compressed VCF with bcftools-based annotations added"
        mut_annot_vep_vcf: "bgzip-compressed VCF with Ensembl VEP annotations added"
        mut_annot_vcf: "bgzip-compressed VCF with INFO fields merged from all annotated VCFs"
        mut_annot_vcf_index: "index of mut_annot_vcf"
        mut_duckdb: "SQL schema and data Parquet files produced by vcf-to-duckdb for loading into DuckDB"
        mut_somatic_variants: "Parquet file of selected somatic variants"
        mut_sig_variants: "Parquet file of unfiltered variants (for mutational signature analysis)"
        mut_enriched_variants: "Parquet file of enriched variants"
        guide_bed_avana: "BED of mutations overlapping Avana guides"
        guide_bed_brunello: "BED of mutations overlapping Brunello guides"
        guide_bed_humagne: "BED of mutations overlapping Humagne guides"
        guide_bed_ky: "BED of mutations overlapping KY guides"
        guide_bed_tko: "BED of mutations overlapping TKO guides"
        mut_sig_scores: "mutational signature scores"
        msisensor2_score: "MSIsensor2 percentage of unstable microsatellite sites"
        msisensor2_output: "MSIsensor2 summary output file"
        msisensor2_output_dist: "MSIsensor2 per-locus repeat length distributions"
        msisensor2_output_somatic: "MSIsensor2 list of somatic microsatellite sites"
        sv_candidate_vcf: "candidate SV VCF"
        sv_candidate_vcf_index: "index of sv_candidate_vcf"
        sv_candidate_indel_vcf: "candidate small indel VCF"
        sv_candidate_indel_vcf_index: "index of sv_candidate_indel_vcf"
        sv_somatic_vcf: "somatic SV VCF"
        sv_somatic_vcf_index: "index of sv_somatic_vcf"
        sv_annot_vcf: "annotated SV VCF"
        sv_annot_bedpe: "annotated SV BEDPE"
        sv_annot_reannotated_bedpe: "re-annotated SV BEDPE"
        sv_del_annotation: "deletion annotations"
        sv_dup_annotation: "duplication annotations"
        sv_selected_somatic: "selected somatic SVs"
    }

    input {
        String sample_id
        File bam
        File bai
        File ref_fasta
        File ref_fasta_index
        File ref_dict
    }

    call call_cnvs.call_cnvs {
        input:
            sample_id = sample_id,
            cram_bam = bam,
            ref_fasta = ref_fasta
    }

    call call_select_svs.call_select_svs {
        input:
            sample_id = sample_id,
            tumor_bam = bam,
            tumor_bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    call infer_msi_status.infer_msi_status {
        input:
            sample_id = sample_id,
            bam = bam,
            bai = bai
    }

    call call_mutations.call_mutations {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_index,
            ref_dict = ref_dict,
            tumor_reads = bam,
            tumor_reads_index = bai
    }

    call patch_hg38_bam_call_muts.patch_hg38_bam_call_muts {
        input:
            sample_id = sample_id,
            bam = bam,
            bai = bai
    }

    call prep_annotate_select_muts.prep_annotate_select_muts {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_vcf = call_mutations.filtered_vcf,
            vcf_patch = patch_hg38_bam_call_muts.filtered_vcf
    }

    call make_guide_mutation_beds.make_guide_mutation_beds {
        input:
            sample_id = sample_id,
            vcf = prep_annotate_select_muts.mut_annot_bcftools_fixed_vcf
    }

    call compute_molecular_signatures.compute_molecular_signatures {
        input:
            sample_id = sample_id,
            mut_sig_variants = select_first([prep_annotate_select_muts.mut_sig_variants])
    }

    output {
        File cnv_bin_coverage = call_cnvs.bin_coverage
        File cnv_segments = call_cnvs.segments
        File cnv_read_cov_bin = call_cnvs.read_cov_bin
        File cnv_input_params = call_cnvs.input_params
        File cnv_cn_by_gene_weighted_mean = call_cnvs.cn_by_gene_weighted_mean
        File cnv_segments_plot = call_cnvs.segments_plot

        File sv_candidate_vcf = call_select_svs.sv_candidate_vcf
        File sv_candidate_vcf_index = call_select_svs.sv_candidate_vcf_index
        File sv_candidate_indel_vcf = call_select_svs.sv_candidate_indel_vcf
        File sv_candidate_indel_vcf_index = call_select_svs.sv_candidate_indel_vcf_index
        File sv_somatic_vcf = call_select_svs.sv_somatic_vcf
        File sv_somatic_vcf_index = call_select_svs.sv_somatic_vcf_index
        File sv_annot_vcf = call_select_svs.sv_annot_vcf
        File sv_annot_bedpe = call_select_svs.sv_annot_bedpe
        File sv_annot_reannotated_bedpe = call_select_svs.sv_annot_reannotated_bedpe
        File sv_del_annotation = call_select_svs.sv_del_annotation
        File sv_dup_annotation = call_select_svs.sv_dup_annotation
        File sv_selected_somatic = call_select_svs.sv_selected_somatic

        Float msisensor2_score = infer_msi_status.msisensor2_score
        File msisensor2_output = infer_msi_status.msisensor2_output
        File msisensor2_output_dist = infer_msi_status.msisensor2_output_dist
        File msisensor2_output_somatic = infer_msi_status.msisensor2_output_somatic

        File mut_vcf = call_mutations.filtered_vcf
        File mut_vcf_idx = call_mutations.filtered_vcf_idx
        File mut_filtering_stats = call_mutations.filtering_stats

        File patched_bam = patch_hg38_bam_call_muts.patched_bam
        File patched_bai = patch_hg38_bam_call_muts.patched_bai
        File mut_vcf_patch = patch_hg38_bam_call_muts.filtered_vcf
        File mut_vcf_patch_idx = patch_hg38_bam_call_muts.filtered_vcf_idx

        File mut_annot_bcftools_fixed_vcf = prep_annotate_select_muts.mut_annot_bcftools_fixed_vcf
        File? mut_annot_bcftools_vcf = prep_annotate_select_muts.mut_annot_bcftools_vcf
        File? mut_annot_vep_vcf = prep_annotate_select_muts.mut_annot_vep_vcf
        File mut_annot_vcf = prep_annotate_select_muts.mut_annot_vcf
        File mut_annot_vcf_index = prep_annotate_select_muts.mut_annot_vcf_index
        Array[File] mut_duckdb = prep_annotate_select_muts.mut_duckdb
        File mut_somatic_variants = prep_annotate_select_muts.mut_somatic_variants
        File? mut_sig_variants = prep_annotate_select_muts.mut_sig_variants
        File? mut_enriched_variants = prep_annotate_select_muts.mut_enriched_variants

        File guide_bed_avana = make_guide_mutation_beds.guide_bed_avana
        File guide_bed_brunello = make_guide_mutation_beds.guide_bed_brunello
        File guide_bed_humagne = make_guide_mutation_beds.guide_bed_humagne
        File guide_bed_ky = make_guide_mutation_beds.guide_bed_ky
        File guide_bed_tko = make_guide_mutation_beds.guide_bed_tko

        File mut_sig_scores = compute_molecular_signatures.mut_sig_scores
    }
}
