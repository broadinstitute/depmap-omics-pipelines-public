version 1.0

workflow quantify_sr_rna {
    input {
        String sample_id
        File transcriptome_bam
        File targets
        File gtf
        Boolean stranded
    }

    if (stranded) {
        call quantify_with_salmon as quantify_with_salmon_auto {
            input:
                sample_id = sample_id,
                transcriptome_bam = transcriptome_bam,
                gtf = gtf,
                targets = targets,
                lib_type = "A"
        }
    }

    call quantify_with_salmon as quantify_with_salmon_iu {
        input:
            sample_id = sample_id,
            transcriptome_bam = transcriptome_bam,
            gtf = gtf,
            targets = targets,
            lib_type = "IU"
    }

    output {
        File? quant_transcripts_auto = quantify_with_salmon_auto.quant_transcripts
        File? quant_genes_auto = quantify_with_salmon_auto.quant_genes
        File quant_genes_iu = quantify_with_salmon_iu.quant_genes
        File quant_transcripts_iu = quantify_with_salmon_iu.quant_transcripts
    }
}

task quantify_with_salmon {
    input {
        String sample_id
        File transcriptome_bam
        File targets
        File gtf
        String lib_type

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/salmon"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 16
        Int cpu = 8
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
        Int additional_mem_gb = 0
    }

    Int disk_space = ceil(
        size(transcriptome_bam, "GiB") + size(targets, "GiB") + size(gtf, "GiB")
    ) + 10 + additional_disk_gb

    Int n_threads = cpu - 1

    command <<<
        set -euo pipefail

        salmon quant \
            --threads ~{n_threads} \
            --targets "~{targets}" \
            --libType "~{lib_type}" \
            --alignments "~{transcriptome_bam}" \
            --geneMap "~{gtf}" \
            --numBootstraps 20 \
            --seqBias \
            --gcBias \
            --posBias \
            --dumpEq \
            --useEM \
            --output "~{sample_id}.salmon_quant"

        mv "~{sample_id}.salmon_quant/quant.sf" "~{sample_id}.salmon_quant_transcripts.tsv"
        mv "~{sample_id}.salmon_quant/quant.genes.sf" "~{sample_id}.salmon_quant_genes.tsv"
    >>>

    output {
        File quant_transcripts = "~{sample_id}.salmon_quant_transcripts.tsv"
        File quant_genes = "~{sample_id}.salmon_quant_genes.tsv"
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
