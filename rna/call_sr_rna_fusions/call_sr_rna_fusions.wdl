version 1.0

workflow call_sr_rna_fusions {
    input {
        String sample_id
        File analysis_ready_bam
        File ref_fasta
        File gtf
    }

    call call_fusions_with_arriba {
        input:
            sample_id = sample_id,
            analysis_ready_bam = analysis_ready_bam,
            ref_fasta = ref_fasta,
            gtf = gtf
    }

    output {
        File fusions = call_fusions_with_arriba.fusions
        File fusions_discarded = call_fusions_with_arriba.fusions_discarded
    }
}

task call_fusions_with_arriba {
    input {
        String sample_id
        File analysis_ready_bam
        File ref_fasta
        File gtf

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/star_arriba"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 16
        Int cpu = 1
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
        Int additional_mem_gb = 0
    }

    Int disk_space = ceil(
        size(analysis_ready_bam, "GiB") + size(ref_fasta, "GiB") + size(gtf, "GiB")
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        arriba \
            -x "~{analysis_ready_bam}" \
            -o "~{sample_id}.arriba_out_fusions.tsv" \
            -O "~{sample_id}.arriba_out_fusions.discarded.tsv" \
            -a "~{ref_fasta}" \
            -g "~{gtf}" \
            -b /arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz
    >>>

    output {
        File fusions = "~{sample_id}.arriba_out_fusions.tsv"
        File fusions_discarded = "~{sample_id}.arriba_out_fusions.discarded.tsv"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_space + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    meta {
        allowNestedInputs: true
    }
}
