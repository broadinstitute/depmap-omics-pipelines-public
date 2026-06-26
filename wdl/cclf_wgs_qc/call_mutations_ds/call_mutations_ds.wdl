version 1.0

workflow call_mutations_ds {
    input {
        # basic inputs
        String sample_id
        File cram_bam
        File crai_bai
        File ref_fasta
        File ref_fasta_index
    }

    call do_call_mutations_ds {
        input:
            sample_id = sample_id,
            cram_bam = cram_bam,
            crai_bai = crai_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    call convert_vcf_to_maf {
        input:
            sample_id = sample_id,
            vcf = do_call_mutations_ds.vcf
    }

    output {
        File vcf = do_call_mutations_ds.vcf
        File vcf_index = do_call_mutations_ds.vcf_index
        File vcf_stats_report = do_call_mutations_ds.vcf_stats_report
        File maf = convert_vcf_to_maf.maf
    }
}

task do_call_mutations_ds {
    input {
        String sample_id
        File cram_bam
        File crai_bai
        File ref_fasta
        File ref_fasta_index

        String docker_image = "google/deepsomatic"
        String docker_image_hash_or_tag = ":1.10.0-gpu"
        Int mem_gb = 64
        Int cpu = 16
        Int preemptible = 0
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 200 + additional_disk_gb # need a lot of space for intermediates

    command <<<
        set -euo pipefail

        mkdir -p ./intermediates
        mkdir -p ./logs

        run_deepsomatic \
            --model_type="WGS_TUMOR_ONLY" \
            --ref="~{ref_fasta}" \
            --reads_tumor="~{cram_bam}" \
            --output_vcf="~{sample_id}.vcf.gz" \
            --sample_name_tumor="~{sample_id}" \
            --num_shards=~{cpu} \
            --intermediate_results_dir="./intermediates" \
            --logging_dir="./logs" \
            --use_default_pon_filtering="true" \
            --vcf_stats_report="true"

        bcftools index "~{sample_id}.vcf.gz" --output "~{sample_id}.vcf.gz.csi"
    >>>

    output {
        File vcf = "~{sample_id}.vcf.gz"
        File vcf_index = "~{sample_id}.vcf.gz.csi"
        File vcf_stats_report = "~{sample_id}.visual_report.html"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
    }

    meta {
        allowNestedInputs: true
    }
}

task convert_vcf_to_maf {
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

    Int disk_space = ceil(2 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        mafsmith vcf2maf \
            --input "~{vcf}" \
            --output "~{sample_id}.maf" \
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

    meta {
        allowNestedInputs: true
    }
}
