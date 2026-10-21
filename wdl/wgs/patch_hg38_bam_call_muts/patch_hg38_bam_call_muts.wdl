version 1.0

import "../patch_hg38_bam/patch_hg38_bam.wdl" as patch_hg38_bam
import "../call_mutations/call_mutations.wdl" as call_mutations

workflow patch_hg38_bam_call_muts {
    meta {
        description: "patch a BAM or CRAM by realigning problematic regions, then call somatic mutations with Mutect2"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        bam: "input BAM file (provide bam+bai or cram+crai)"
        bai: "index of input BAM file"
        cram: "input CRAM file (provide bam+bai or cram+crai)"
        crai: "index of input CRAM file"
        old_ref_fasta: "reference FASTA that the input BAM/CRAM was aligned to; required to decode CRAM input"
        old_ref_fasta_index: "index of old_ref_fasta"
        problems_bed: "BED file of problematic regions to realign"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_fasta_amb: "BWA .amb index file of ref_fasta"
        ref_fasta_ann: "BWA .ann index file of ref_fasta"
        ref_fasta_bwt: "BWA .bwt index file of ref_fasta"
        ref_fasta_pac: "BWA .pac index file of ref_fasta"
        ref_fasta_sa: "BWA .sa index file of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        scatter_count: "number of intervals to scatter Mutect2 over"

        # outputs
        filtered_vcf: "bgzip-compressed filtered somatic variant VCF from Mutect2"
        filtered_vcf_idx: "index of filtered_vcf"
    }

    input {
        String sample_id
        File? bam
        File? bai
        File? cram
        File? crai
        File? old_ref_fasta
        File? old_ref_fasta_index
        File problems_bed
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa
        File ref_dict
        Int scatter_count
    }

    call patch_hg38_bam.patch_hg38_bam {
        input:
            sample_id = sample_id,
            bam = bam,
            bai = bai,
            cram = cram,
            crai = crai,
            old_ref_fasta = old_ref_fasta,
            old_ref_fasta_index = old_ref_fasta_index,
            problems_bed = problems_bed,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_fasta_amb = ref_fasta_amb,
            ref_fasta_ann = ref_fasta_ann,
            ref_fasta_bwt = ref_fasta_bwt,
            ref_fasta_pac = ref_fasta_pac,
            ref_fasta_sa = ref_fasta_sa
    }

    call call_mutations.call_mutations {
        input:
            sample_id = sample_id,
            intervals = problems_bed,
            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            tumor_reads = patch_hg38_bam.patched_bam,
            tumor_reads_index = patch_hg38_bam.patched_bai
    }

    output {
        File patched_bam = patch_hg38_bam.patched_bam
        File patched_bai = patch_hg38_bam.patched_bai
        File filtered_vcf = call_mutations.filtered_vcf
        File filtered_vcf_idx = call_mutations.filtered_vcf_idx
    }
}
