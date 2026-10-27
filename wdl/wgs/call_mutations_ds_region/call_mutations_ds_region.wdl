version 1.0

workflow call_mutations_ds_region {
    meta {
        description: "Call somatic mutations with DeepSomatic in tumor-only mode"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        cram_bam: "tumor BAM or CRAM file"
        crai_bai: "index of cram_bam"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        regions: "genomic regions to process"

        # outputs
        vcf: "somatic variant VCF"
        vcf_index: "index of vcf"
        gvcf: "gVCF"
        gvcf_index: "index of gvcf"
        maf: "somatic variants in MAF format"
        vcf_stats_report: "DeepSomatic VCF stats HTML report"
    }

    input {
        String sample_id
        File cram_bam
        File crai_bai
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Array[String] regions
    }

    call do_call_mutations_ds {
        input:
            sample_id = sample_id,
            cram_bam = cram_bam,
            crai_bai = crai_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            regions = regions
    }

    call convert_vcf_to_maf {
        input:
            sample_id = sample_id,
            vcf = do_call_mutations_ds.vcf
    }

    output {
        File vcf_ds = do_call_mutations_ds.vcf
        File vcf_ds_index = do_call_mutations_ds.vcf_index
        File gvcf = do_call_mutations_ds.gvcf
        File gvcf_index = do_call_mutations_ds.gvcf_index
        File maf = convert_vcf_to_maf.maf
    }
}

task do_call_mutations_ds {
    meta {
        description: "Run DeepSomatic on a single interval shard in tumor-only mode"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        cram_bam: "tumor BAM or CRAM file"
        crai_bai: "index of cram_bam"
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        regions: "genomic regions to process"

        # outputs
        vcf: "somatic variant VCF"
        gvcf: "gVCF"
    }

    input {
        String sample_id
        File cram_bam
        File crai_bai
        File ref_fasta
        File ref_fasta_index
        Array[String] regions

        String docker_image = "google/deepsomatic"
        String docker_image_hash_or_tag = ":1.10.0"
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 100 + additional_disk_gb # need a lot of space for intermediates

    Int n_shards = cpu
    Int seq2_arg = cpu - 1

    command <<<
        set -euo pipefail

        mkdir -p ./intermediates

        seq 0 ~{seq2_arg} | parallel -q \
            --halt 2 \
            --line-buffer /opt/deepvariant/bin/make_examples_somatic \
            --mode calling \
            --ref "~{ref_fasta}" \
            --reads_tumor "~{cram_bam}" \
            --examples "./intermediates/make_examples_somatic.tfrecord@~{n_shards}.gz" \
            --checkpoint "/opt/models/deepsomatic/wgs_tumor_only" \
            --checkpoint_json "/opt/models/deepsomatic/wgs_tumor_only/model.example_info.json" \
            --gvcf "./intermediates/gvcf.tfrecord@~{n_shards}.gz" \
            --regions "~{sep=" " regions}" \
            --sample_name_tumor "~{sample_id}" \
            --task {}

        /opt/deepvariant/bin/call_variants \
            --outfile "./intermediates/call_variants_output.tfrecord.gz" \
            --examples "./intermediates/make_examples_somatic.tfrecord@~{n_shards}.gz" \
            --checkpoint "/opt/models/deepsomatic/wgs_tumor_only"

        /opt/deepvariant/bin/postprocess_variants \
            --ref "~{ref_fasta}" \
            --infile "./intermediates/call_variants_output@1.tfrecord.gz" \
            --outfile "~{sample_id}.vcf.gz" \
            --process_somatic=true \
            --pon_filtering "/opt/models/deepsomatic/pons/PON_dbsnp138_gnomad_ILMN1000g_pon.vcf.gz" \
            --gvcf_outfile "~{sample_id}.g.vcf.gz" \
            --nonvariant_site_tfrecord_path "./intermediates/gvcf.tfrecord@~{n_shards}.gz"

        bcftools index "~{sample_id}.vcf.gz" --output="~{sample_id}.vcf.gz.csi"
        bcftools index "~{sample_id}.g.vcf.gz" --output="~{sample_id}.g.vcf.gz.csi"
    >>>

    output {
        File vcf = "~{sample_id}.vcf.gz"
        File vcf_index = "~{sample_id}.vcf.gz.csi"
        File gvcf = "~{sample_id}.g.vcf.gz"
        File gvcf_index = "~{sample_id}.g.vcf.gz.csi"
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

task convert_vcf_to_maf {
    meta {
        description: "Convert a somatic VCF to MAF format"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        vcf: "somatic variant VCF to convert"

        # outputs
        maf: "somatic variants in MAF format"
    }

    input {
        String sample_id
        File vcf

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/mafsmith"
        String docker_image_hash_or_tag = ":0.1.0"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        mafsmith vcf2maf \
            -i "~{vcf}" \
            -o "~{sample_id}.maf" \
            --tumor-id "~{sample_id}" \
            --skip-annotation
    >>>

    output {
        File maf = "~{sample_id}.maf"
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
