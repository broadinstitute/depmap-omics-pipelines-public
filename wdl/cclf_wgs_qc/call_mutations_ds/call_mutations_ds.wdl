version 1.0

workflow call_mutations_ds {
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
        scatter_count: "number of intervals to scatter DeepSomatic over"

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
        Int scatter_count = 1
    }

    call split_intervals {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = scatter_count
    }

    scatter (regions in split_intervals.intervals) {
        call do_call_mutations_ds {
            input:
                sample_id = sample_id,
                cram_bam = cram_bam,
                crai_bai = crai_bai,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                regions = regions
        }
    }

    call merge_vcfs {
        input:
            input_vcfs = do_call_mutations_ds.vcf,
            output_basename = "~{sample_id}"
    }

    call merge_vcfs as merge_gvcfs {
        input:
            input_vcfs = do_call_mutations_ds.gvcf,
            output_basename = "~{sample_id}.g"
    }

    call make_stats_report {
        input:
            sample_id = sample_id,
            vcf = merge_vcfs.merged_vcf
    }

    call convert_vcf_to_maf {
        input:
            sample_id = sample_id,
            vcf = merge_vcfs.merged_vcf
    }

    output {
        File vcf = merge_vcfs.merged_vcf
        File vcf_index = merge_vcfs.merged_vcf_index
        File gvcf = merge_gvcfs.merged_vcf
        File gvcf_index = merge_gvcfs.merged_vcf_index
        File maf = convert_vcf_to_maf.maf
        File vcf_stats_report = make_stats_report.vcf_stats_report
    }
}

task split_intervals {
    meta {
        description: "Split genomic intervals for scattering DeepSomatic"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        scatter_count: "number of output interval groups to produce"

        # outputs
        intervals: "scattered intervals as tab-separated region strings, one group per row"
    }

    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Int scatter_count

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk"
        String docker_image_hash_or_tag = ":4.6.1.0"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = 12 + additional_disk_gb

    command <<<
        set -euo pipefail

        mkdir -p interval-files

        # split intervals (of at least size(chr21))
        java -jar /root/gatk.jar SplitIntervals \
            -R "~{ref_fasta}" \
            -scatter ~{scatter_count} \
            -O interval-files \
            --min-contig-size 46709983

        for f in interval-files/*-scattered.interval_list; do
            awk '
                BEGIN { first = 1 }
                    !/^@/ {
                        if (!first) printf "\t";
                            printf "%s:%s-%s", $1, $2, $3;
                            first = 0;
                        }
                END { printf "\n" }
            ' "$f"
        done > intervals.tsv
    >>>

    output {
        Array[Array[String]] intervals = read_tsv("intervals.tsv")
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
        regions: "genomic regions to process for this shard"

        # outputs
        vcf: "somatic variant VCF for this shard"
        gvcf: "gVCF for this shard"
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
    >>>

    output {
        File vcf = "~{sample_id}.vcf.gz"
        File gvcf = "~{sample_id}.g.vcf.gz"
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

task merge_vcfs {
    meta {
        description: "Merge scattered VCF shards into a single VCF"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_vcfs: "VCF shards to merge"
        output_basename: "base name for the merged output VCF"

        # outputs
        merged_vcf: "merged VCF"
        merged_vcf_index: "index of merged_vcf"
    }

    input {
        Array[File] input_vcfs
        String output_basename

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk"
        String docker_image_hash_or_tag = ":4.6.1.0"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(3 * size(input_vcfs, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        java -jar /root/gatk.jar MergeVcfs \
            -I ~{sep=' -I ' input_vcfs} \
            -O "~{output_basename}.vcf.gz" \
            --CREATE_INDEX
    >>>

    output {
        File merged_vcf = "~{output_basename}.vcf.gz"
        File merged_vcf_index = "~{output_basename}.vcf.gz.tbi"
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

task make_stats_report {
    meta {
        description: "Generate a DeepSomatic VCF stats HTML report"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        vcf: "somatic variant VCF to report on"

        # outputs
        vcf_stats_report: "HTML report of VCF variant statistics"
    }

    input {
        String sample_id
        File vcf

        String docker_image = "google/deepsomatic"
        String docker_image_hash_or_tag = ":1.10.0"
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 1
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(size(vcf, "GiB")) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        /opt/deepvariant/bin/vcf_stats_report \
            --input_vcf "~{vcf}" \
            --outfile_base "~{sample_id}"
    >>>

    output {
        File vcf_stats_report = "~{sample_id}.visual_report.html"
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
