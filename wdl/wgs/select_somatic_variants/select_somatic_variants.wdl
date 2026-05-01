version 1.0

workflow select_somatic_variants {
    meta {
        description: "Select somatic variants from an annotated DuckDB Parquet dataset"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        duckdb: "Parquet files produced by vcf-to-duckdb"
        min_af: "minimum allele frequency for a variant to be retained"
        min_depth: "minimum read depth for a variant to be retained"
        max_pop_af: "maximum population allele frequency"
        save_enriched_variants: "if true, also save the intermediate enriched variants data as Parquet"

        # outputs
        mut_somatic_variants: "Parquet file of selected somatic variants"
        mut_enriched_variants: "Parquet file of unfiltered variants (with all output and debugging columns)"
        mut_sig_variants: "Parquet file of unfiltered variants (for mutational signature analysis)"
    }

    input {
        String sample_id
        Array[File] duckdb
        Float min_af = 0.15
        Int min_depth = 5
        Float max_pop_af = 0.00001
        Boolean save_enriched_variants = false
    }

    call do_select_somatic_variants {
        input:
            sample_id = sample_id,
            duckdb = duckdb,
            min_af = min_af,
            min_depth = min_depth,
            max_pop_af = max_pop_af,
            save_enriched_variants = save_enriched_variants
    }

    output {
        File mut_somatic_variants = do_select_somatic_variants.somatic_variants
        File? mut_enriched_variants = do_select_somatic_variants.enriched_variants
        File mut_sig_variants = do_select_somatic_variants.sig_variants
    }
}

task do_select_somatic_variants {
    meta {
        description: "Select somatic variants from an annotated DuckDB Parquet dataset"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        duckdb: "Parquet files produced by vcf-to-duckdb"
        min_af: "minimum allele frequency for a variant to be retained"
        min_depth: "minimum read depth for a variant to be retained"
        max_pop_af: "maximum population allele frequency"
        save_enriched_variants: "if true, also save the intermediate enriched variants data as Parquet"
        batch_size: "number of variants to process per batch (default: 300000)"

        # outputs
        somatic_variants: "Parquet file of selected somatic variants"
        enriched_variants: "Parquet file of unfiltered variants (with all output and debugging columns)"
        sig_variants: "Parquet file of unfiltered variants (for mutational signature analysis)"
    }

    input {
        String sample_id
        Array[File] duckdb
        Float min_af
        Int min_depth
        Float max_pop_af
        Boolean save_enriched_variants
        Int batch_size = 300000

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/select-somatic-variants"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 32
        Int cpu = 8
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(10 * size(duckdb, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        mkdir -p parq
        mv ~{sep=" " duckdb} parq/

        mkdir -p db_tmp

        if [ "~{save_enriched_variants}" = "true" ]
        then
            variants_enriched_line="--variants-enriched-out-file=~{sample_id}.enriched_variants.parquet"
        else
            variants_enriched_line=""
        fi

        python -m select_somatic_variants \
            --db="tmp.duckdb" \
            --parquet-dir="./parq" \
            --somatic-variants-out-file="~{sample_id}.somatic_variants.parquet" \
            $variants_enriched_line \
            --mut-sig-out-file="~{sample_id}.mut_sig_variants.parquet" \
            --db-tmp-dir-path="./db_tmp/" \
            --sample-id="~{sample_id}" \
            --min-af=~{min_af} \
            --min-depth=~{min_depth} \
            --max-pop-af=~{max_pop_af} \
            --batch-size=~{batch_size}
    >>>

    output {
        File somatic_variants = "~{sample_id}.somatic_variants.parquet"
        File? enriched_variants = "~{sample_id}.enriched_variants.parquet"
        File sig_variants = "~{sample_id}.mut_sig_variants.parquet"
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
