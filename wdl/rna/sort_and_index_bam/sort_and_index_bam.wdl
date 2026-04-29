version 1.0

workflow sort_and_index_bam {
    meta {
        description: "Sort and index a BAM file"
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_bam: "input BAM file"

        # outputs
        output_bam: "sorted BAM"
        output_bai: "index of output_bam"
    }

    input {
        String sample_id
        File input_bam
    }

    call do_sort_and_index_bam {
        input:
            sample_id = sample_id,
            input_bam = input_bam
    }

    output {
        File output_bam = do_sort_and_index_bam.output_bam
        File output_bai = do_sort_and_index_bam.output_bai
    }
}

task do_sort_and_index_bam {
    meta {
        description: "Coordinate-sort a BAM file and index it using samtools"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_bam: "input BAM file to sort and index"

        # outputs
        output_bam: "coordinate-sorted BAM"
        output_bai: "index of sorted_bam"
    }

    input {
        String sample_id
        File input_bam

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 8
        Int preemptible = 1
        Int max_retries = 0
        Int additional_memory_gb = 0
        Int additional_disk_gb = 0
    }

    Int mem_gb = cpu + 4 + additional_memory_gb
    Int disk_space = ceil(3 * size(input_bam, "GiB")) + 20 + additional_disk_gb
    Int mem_per_thread = ((mem_gb - 4) * 1024) / (cpu + 1)
    String output_bam = sample_id + ".bam"
    String output_bai = sample_id + ".bai"

    command <<<
        set -euo pipefail

        echo "Sorting BAM"
        samtools sort \
            -@ ~{cpu} \
            -o ~{output_bam} \
            -T "tmp_srt" \
            -m ~{mem_per_thread}M \
            ~{input_bam}

        echo "Indexing BAM"
        samtools index \
            -@ ~{cpu} \
            -b \
            "~{sample_id}.bam" \
            "~{sample_id}.bai"
    >>>

    output {
        File output_bam = "~{sample_id}.bam"
        File output_bai = "~{sample_id}.bai"
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
