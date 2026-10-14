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
    }

    call patch_hg38_bam.patch_hg38_bam {
        input:
            sample_id = sample_id,
            bam = bam,
            bai = bai,
            cram = cram,
            crai = crai,
            old_ref_fasta = old_ref_fasta,
            old_ref_fasta_index = old_ref_fasta_index
    }

    call call_mutations.call_mutations {
        input:
            sample_id = sample_id,
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
