version 1.0

# Adapted from https://github.com/broadinstitute/gatk/blob/4.5.0.0/scripts/mutect2_wdl/mutect2.wdl
#
# Modifications:
#   - Remove tasks and options related to mutect3 training data
#   - Use SSDs instead of HDDs by default
#
# Use this file as a base and manually implement changes in the upstream code in order
# to retain these modifications.

struct Runtime {
    String gatk_docker
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
}

workflow call_mutations {
    meta {
        description: "Call mutations with Mutect2"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        ref_fasta: "reference FASTA"
        ref_fai: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        tumor_reads: "tumor BAM or CRAM file"
        tumor_reads_index: "index of tumor_reads"
        pon: "panel of normals VCF for filtering artifacts"
        pon_idx: "index of pon"
        gnomad: "gnomAD VCF for germline resource"
        gnomad_idx: "index of gnomad"
        m2_extra_args: "additional arguments to pass to Mutect2"
        m2_extra_filtering_args: "additional arguments to pass to FilterMutectCalls"
        run_orientation_bias_mixture_model_filter: "if true, run the orientation bias mixture model filter"
        compress_vcfs: "if true, output gzip-compressed VCFs"
        gatk_docker: "GATK Docker image"
        scatter_count: "number of intervals to scatter Mutect2 over"
        preemptible: "number of preemptible retries"
        max_retries: "number of non-preemptible retries"
        small_task_cpu: "CPU for small tasks"
        small_task_mem: "memory in GB for small tasks"
        small_task_disk: "disk in GB for small tasks"
        learn_read_orientation_mem: "memory in MB for LearnReadOrientationModel"
        emergency_extra_disk: "additional disk added to every task as a last resort"

        # outputs
        filtered_vcf: "filtered somatic variant VCF"
        filtered_vcf_idx: "index of filtered_vcf"
        filtering_stats: "FilterMutectCalls filtering statistics"
    }

    input {
        # basic inputs
        String sample_id
        File ref_fasta
        File ref_fai
        File ref_dict
        File tumor_reads
        File tumor_reads_index

        # optional but usually recommended resources
        File? pon
        File? pon_idx
        File? gnomad
        File? gnomad_idx

        # extra arguments
        String? m2_extra_args
        String? m2_extra_filtering_args

        # additional modes and outputs
        Boolean run_orientation_bias_mixture_model_filter = true
        Boolean compress_vcfs = true

        # runtime
        String gatk_docker = "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk:4.6.1.0"
        Int scatter_count
        Int preemptible = 1
        Int max_retries = 1
        Int small_task_cpu = 1
        Int small_task_mem = 4
        Int small_task_disk = 15
        Int learn_read_orientation_mem = 8000

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int emergency_extra_disk = 0
    }

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
    Int tumor_reads_size = ceil(size(tumor_reads, "GB") + size(tumor_reads_index, "GB"))
    Int gnomad_vcf_size = if defined(gnomad) then ceil(size(gnomad, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + emergency_extra_disk

    Runtime standard_runtime = {"gatk_docker": gatk_docker,
            "max_retries": max_retries, "preemptible": preemptible, "cpu": small_task_cpu,
            "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
            "disk": small_task_disk + disk_pad}

    Int m2_output_size = tumor_reads_size / scatter_count
    Int m2_per_scatter_size = ref_size + gnomad_vcf_size + m2_output_size + disk_pad

    call SplitIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            runtime_params = standard_runtime
    }

    scatter (subintervals in SplitIntervals.interval_files) {
        call M2 {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_reads = tumor_reads,
                tumor_reads_index = tumor_reads_index,
                pon = pon,
                pon_idx = pon_idx,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                mem = small_task_mem,
                preemptible = 2,
                max_retries = max_retries,
                m2_extra_args = m2_extra_args,
                run_ob_filter = run_orientation_bias_mixture_model_filter,
                compress_vcfs = compress_vcfs,
                gatk_docker = gatk_docker,
                disk_space = m2_per_scatter_size
        }
    }

    if (run_orientation_bias_mixture_model_filter) {
        call LearnReadOrientationModel {
            input:
                f1r2_tar_gz = M2.f1r2_counts,
                runtime_params = standard_runtime,
                mem = learn_read_orientation_mem
        }
    }

    call MergeVCFs {
        input:
            input_vcfs = M2.unfiltered_vcf,
            input_vcf_indices = M2.unfiltered_vcf_idx,
            compress_vcfs = compress_vcfs,
            runtime_params = standard_runtime
    }

    call MergeStats {
        input:
            stats = M2.stats,
            runtime_params = standard_runtime
    }

    call Filter {
        input:
            sample_id = sample_id,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            unfiltered_vcf = MergeVCFs.merged_vcf,
            unfiltered_vcf_idx = MergeVCFs.merged_vcf_idx,
            compress_vcfs = compress_vcfs,
            mutect_stats = MergeStats.merged_stats,
            artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
            m2_extra_filtering_args = m2_extra_filtering_args,
            runtime_params = standard_runtime,
            disk_space = ceil(size(MergeVCFs.merged_vcf, "GB") * 4) + disk_pad
    }

    output {
        File filtered_vcf = Filter.filtered_vcf
        File filtered_vcf_idx = Filter.filtered_vcf_idx
        File filtering_stats = Filter.filtering_stats
    }
}

task SplitIntervals {
    meta {
        description: "Split genomic intervals for scattering Mutect2"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        intervals: "optional interval list to restrict splitting"
        ref_fasta: "reference FASTA"
        ref_fai: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        scatter_count: "number of output interval files to produce"
        runtime_params: "runtime parameters struct"

        # outputs
        interval_files: "scattered interval list files"
    }

    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count

      # runtime
      Runtime runtime_params
    }

    command <<<
        set -eou pipefail

        mkdir interval-files

        java -Xmx~{runtime_params.command_mem}m -jar /root/gatk.jar SplitIntervals \
            -R ~{ref_fasta} \
            ~{"-L " + intervals} \
            -scatter ~{scatter_count} \
            -O interval-files

        cp interval-files/*.interval_list .
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    meta {
        description: "Run Mutect2 on a single interval shard"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        intervals: "interval list for this scatter shard"
        ref_fasta: "reference FASTA"
        ref_fai: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        tumor_reads: "tumor BAM or CRAM file"
        tumor_reads_index: "index of tumor_reads"
        pon: "panel of normals VCF for filtering artifacts"
        pon_idx: "index of pon"
        gnomad: "gnomAD VCF for germline resource"
        gnomad_idx: "index of gnomad"
        m2_extra_args: "additional arguments to pass to Mutect2"
        run_ob_filter: "if true, collect F1R2 counts for orientation bias filtering"
        compress_vcfs: "if true, output gzip-compressed VCFs"
        gatk_docker: "GATK Docker image"
        mem: "memory override in GB"
        preemptible: "number of preemptible retries"
        max_retries: "number of non-preemptible retries"
        disk_space: "disk space override in GB"
        cpu: "number of CPUs"

        # outputs
        unfiltered_vcf: "unfiltered Mutect2 VCF for this shard"
        unfiltered_vcf_idx: "index of unfiltered_vcf"
        tumor_sample: "tumor sample name extracted from BAM header"
        normal_sample: "normal sample name (empty for tumor-only)"
        stats: "Mutect2 stats file for this shard"
        f1r2_counts: "F1R2 counts tar archive for orientation bias model"
    }

    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File tumor_reads
        File tumor_reads_index
        File? pon
        File? pon_idx
        File? gnomad
        File? gnomad_idx
        String? m2_extra_args
        Boolean? run_ob_filter
        Boolean compress_vcfs

        # runtime
        String gatk_docker
        Int? mem
        Int? preemptible
        Int? max_retries
        Int? disk_space
        Int? cpu
    }

    String output_vcf = "output" + if compress_vcfs then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then select_first([mem]) * 1000 else 3500
    Int command_mem = machine_mem - 500

    parameter_meta {
        intervals: { localization_optional: true }
        ref_fasta: { localization_optional: true }
        ref_fai: { localization_optional: true }
        ref_dict: { localization_optional: true }
        tumor_reads: { localization_optional: true }
        tumor_reads_index: { localization_optional: true }
        pon: { localization_optional: true }
        pon_idx: { localization_optional: true }
        gnomad: { localization_optional: true }
        gnomad_idx: { localization_optional: true }
    }

    command <<<
        set -eou pipefail

        if project="$(gcloud config get-value project -q 2>/dev/null)"; then
            export GCS_REQUESTER_PAYS_PROJECT="$project"
        fi

        # We need to create these files regardless, even if they stay empty
        touch f1r2.tar.gz
        touch dataset.txt
        echo "" > normal_name.txt

        java -Xmx~{command_mem}m -jar /root/gatk.jar GetSampleName \
            -R "~{ref_fasta}" \
            -I "~{tumor_reads}" \
            -O tumor_name.txt \
            -encode \
            --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"

        java -Xmx~{command_mem}m -jar /root/gatk.jar Mutect2 \
            -R "~{ref_fasta}" \
            -I "~{tumor_reads}" \
            -tumor "$(cat tumor_name.txt)" \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            -O "~{output_vcf}" \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{m2_extra_args} \
            --gcs-project-for-requester-pays "${GCS_REQUESTER_PAYS_PROJECT}"
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " SSD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
    }
}

task MergeVCFs {
    meta {
        description: "Merge scattered Mutect2 VCF shards into a single VCF"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_vcfs: "VCF shards to merge"
        input_vcf_indices: "indices of input_vcfs"
        compress_vcfs: "if true, output gzip-compressed VCF"
        runtime_params: "runtime parameters struct"

        # outputs
        merged_vcf: "merged VCF"
        merged_vcf_idx: "index of merged_vcf"
    }

    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      Boolean compress_vcfs
      Runtime runtime_params
    }

    String output_vcf = if compress_vcfs then "merged.vcf.gz" else "merged.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
    command <<<
        set -eou pipefail

        java -Xmx~{runtime_params.command_mem}m -jar /root/gatk.jar MergeVcfs \
            -I ~{sep=' -I ' input_vcfs} \
            -O ~{output_vcf}
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task MergeStats {
    meta {
        description: "Merge per-shard Mutect2 stats files"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        stats: "per-shard Mutect2 stats files to merge"
        runtime_params: "runtime parameters struct"

        # outputs
        merged_stats: "merged Mutect2 stats file"
    }

    input {
      Array[File] stats
      Runtime runtime_params
    }

    command <<<
        set -eou pipefail

        java -Xmx~{runtime_params.command_mem}m -jar /root/gatk.jar MergeMutectStats \
            -stats ~{sep=" -stats " stats} \
            -O merged.stats
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_stats = "merged.stats"
    }
}

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter as of September 2018
task LearnReadOrientationModel {
    meta {
        description: "Learn read orientation model from F1R2 counts for orientation bias filtering"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        f1r2_tar_gz: "F1R2 counts tar archives from Mutect2 scatter shards"
        runtime_params: "runtime parameters struct"
        mem: "memory override in MB"

        # outputs
        artifact_prior_table: "artifact prior probabilities tar archive for use with FilterMutectCalls"
    }

    input {
      Array[File] f1r2_tar_gz
      Runtime runtime_params
      Int? mem  #override memory
    }

    Int machine_mem = select_first([mem, runtime_params.machine_mem])
    Int command_mem = machine_mem - 1000

    command <<<
        set -eou pipefail

        java -Xmx~{command_mem}m -jar /root/gatk.jar LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }

}

task Filter {
    meta {
        description: "Filter Mutect2 calls with FilterMutectCalls"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        ref_fasta: { localization_optional: true, description: "reference FASTA" }
        ref_fai: { localization_optional: true, description: "index of ref_fasta" }
        ref_dict: { localization_optional: true, description: "sequence dictionary for ref_fasta" }
        unfiltered_vcf: "merged unfiltered Mutect2 VCF"
        unfiltered_vcf_idx: "index of unfiltered_vcf"
        compress_vcfs: "if true, output gzip-compressed VCF"
        mutect_stats: "merged Mutect2 stats file"
        artifact_priors_tar_gz: "orientation bias artifact priors from LearnReadOrientationModel"
        contamination_table: "contamination estimates table"
        maf_segments: "tumor segmentation file for contamination filtering"
        m2_extra_filtering_args: "additional arguments to pass to FilterMutectCalls"
        runtime_params: "runtime parameters struct"
        disk_space: "disk space override in GB"

        # outputs
        filtered_vcf: "filtered somatic variant VCF"
        filtered_vcf_idx: "index of filtered_vcf"
        filtering_stats: "FilterMutectCalls filtering statistics"
    }

    input {
        String sample_id
        File ref_fasta
        File ref_fai
        File ref_dict
        File unfiltered_vcf
        File unfiltered_vcf_idx
        Boolean compress_vcfs
        File? mutect_stats
        File? artifact_priors_tar_gz
        File? contamination_table
        File? maf_segments
        String? m2_extra_filtering_args

        Runtime runtime_params
        Int? disk_space
    }

    String output_vcf = if compress_vcfs then "~{sample_id}.filtered.vcf.gz" else "~{sample_id}.filtered.vcf"
    String output_vcf_idx = output_vcf + if compress_vcfs then ".tbi" else ".idx"

    command <<<
        set -eou pipefail

        java -Xmx~{runtime_params.command_mem}m -jar /root/gatk.jar FilterMutectCalls \
            -V ~{unfiltered_vcf} \
            -R ~{ref_fasta} \
            -O "~{output_vcf}" \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            ~{"--ob-priors " + artifact_priors_tar_gz} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats "~{sample_id}.filtering.stats" \
            ~{m2_extra_filtering_args}
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " SSD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "~{sample_id}.filtering.stats"
    }
}
