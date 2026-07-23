version 1.0

workflow plot_all_cnvs {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    input {
        String sample_set_id
        Array[String] cell_line_names
        Array[File] cnv_read_cov_bins
    }

    call make_cnv_heat_map {
        input:
            sample_set_id = sample_set_id,
            cell_line_names = cell_line_names,
            cnv_read_cov_bins = cnv_read_cov_bins
    }

    output {
        File cnv_heat_map = make_cnv_heat_map.cnv_heat_map
    }
}

task make_cnv_heat_map {
    meta {
        description: "TODO"
        allowNestedInputs: true
    }

    input {
        String sample_set_id
        Array[String] cell_line_names
        Array[File] cnv_read_cov_bins

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/call_cnvs"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 1
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
        Int additional_mem_gb = 0
    }

    Int mem_gb = 4 + ceil(size(cnv_read_cov_bins, "GiB") * 3) + additional_mem_gb
    Int disk_space = ceil(size(cnv_read_cov_bins, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        names=(~{sep=' ' cell_line_names})
        files=(~{sep=' ' cnv_read_cov_bins})

        mkdir ./renamed

        for i in "${!names[@]}"; do
            cp "${files[$i]}" "./renamed/${names[$i]}.tsv.gz"
        done

        Rscript /app/plot_heat_map.R \
            "~{sample_set_id}.png" \
            ./renamed/*.tsv.gz
    >>>

    output {
        File cnv_heat_map = "~{sample_set_id}.png"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
