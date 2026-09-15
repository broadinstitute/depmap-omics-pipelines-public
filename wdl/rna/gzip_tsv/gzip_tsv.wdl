version 1.0

workflow gzip_tsv {
    input {
        File in_tsv
    }

    call gzip_file {
        input:
            in_tsv = in_tsv
    }

    output {
        File out_tsv = gzip_file.out_tsv
    }
}

task gzip_file {
    input {
        File in_tsv

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
        Int additional_mem_gb = 0
        Int disk_space = 10
    }

    command <<<
        gzip -c "~{in_tsv}" > "~{basename(in_tsv)}.gz"
    >>>

    output {
        File out_tsv = "~{basename(in_tsv)}.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
