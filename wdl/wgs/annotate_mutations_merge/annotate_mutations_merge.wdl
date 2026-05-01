version 1.0

workflow annotate_mutations_merge {
    meta {
        description: "Merge INFO fields from multiple annotated VCFs and convert to DuckDB"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_vcfs: "array of annotated VCFs whose INFO fields will be merged into a single VCF"
        vcf_to_duckdb_config: "optional JSON config file controlling vcf-to-duckdb field selection and types"

        # outputs
        mut_annot_vcf: "bgzip-compressed VCF with INFO fields merged from all input VCFs"
        mut_annot_vcf_index: "index of mut_annot_vcf"
        mut_duckdb: "SQL schema and data Parquet files produced by vcf-to-duckdb for loading into DuckDB"
    }

    input {
        String sample_id
        Array[File] input_vcfs
        File? vcf_to_duckdb_config
    }

    call merge_info {
        input:
            vcfs = select_all(input_vcfs),
            output_file_base_name = sample_id + "_info_merged",
    }

    call vcf_to_duckdb {
        input:
            vcf = merge_info.vcf_info_merged,
            vcf_index = merge_info.vcf_info_merged_index,
            config = vcf_to_duckdb_config
    }

    output {
        File mut_annot_vcf = merge_info.vcf_info_merged
        File mut_annot_vcf_index = merge_info.vcf_info_merged_index
        Array[File] mut_duckdb = vcf_to_duckdb.duckdb
    }
}

task merge_info {
    meta {
        description: "Merge INFO fields from multiple annotated VCFs into a single VCF"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        vcfs: "array of annotated VCFs to merge"
        output_file_base_name: "base name for the output file (without extension)"

        # outputs
        vcf_info_merged: "bgzip-compressed VCF with INFO fields merged from all input VCFs"
        vcf_info_merged_index: "index of vcf_info_merged"
    }

    input {
        Array[File] vcfs
        String output_file_base_name

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/vcf-info-merger"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(vcfs, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        python -m vcf_info_merger merge \
            --vcf ~{sep=" --vcf " vcfs} \
            --out="~{output_file_base_name}.vcf.gz"

        bcftools index "~{output_file_base_name}.vcf.gz" \
            --output="~{output_file_base_name}.vcf.gz.csi"
    >>>

    output {
        File vcf_info_merged = "~{output_file_base_name}.vcf.gz"
        File vcf_info_merged_index = "~{output_file_base_name}.vcf.gz.csi"
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

task vcf_to_duckdb {
    meta {
        description: "Convert an annotated VCF to DuckDB-compatible Parquet files"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        vcf: "bgzip-compressed annotated VCF to convert"
        vcf_index: "index of vcf"
        config: "optional JSON config file controlling field selection and types"
        batch_size: "number of variants to process per batch"

        # outputs
        duckdb: "SQL schema and data Parquet files produced by vcf-to-duckdb for loading into DuckDB"
    }

    input {
        File vcf
        File vcf_index
        File? config
        Int batch_size = 1000000

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/vcf-to-duckdb"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 16
        Int cpu = 4
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(vcf, "GiB")) + 30 + additional_disk_gb
    Int duck_db_memory_limit = floor(mem_gb * 0.75)

    command <<<
        set -euo pipefail

        mkdir -p "./data/tmp"

        python -m vcf_to_duckdb convert \
            --vcf="~{vcf}" \
            --tab="tmp.tsv" \
            --db="tmp.duckdb" \
            --parquet-dir="./parq" \
            --no-multiallelics \
            ~{"--config " + config} \
            --batch-size=~{batch_size} \
            --compression-level=15 \
            --tmp-dir="./data/tmp" \
            --memory-limit="~{duck_db_memory_limit}GB"
    >>>

    output {
        Array[File] duckdb = glob("parq/*.*")
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
