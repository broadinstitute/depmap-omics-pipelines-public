version 1.0

# Adapted from https://github.com/broadinstitute/warp/blob/GDCWholeGenomeSomaticSingleSample_v1.3.1/pipelines/broad/dna_seq/somatic/single_sample/wgs/gdc_genome/GDCWholeGenomeSomaticSingleSample.wdl
# and https://github.com/broadinstitute/warp/blob/GDCWholeGenomeSomaticSingleSample_v1.3.1/pipelines/broad/reprocessing/cram_to_unmapped_bams/CramToUnmappedBams.wdl

workflow align_wgs_sample {
    meta {
        description: "Produce an analysis-ready, optionally BQSR-recalibrated BAM from CRAM, BAM, or FASTQ input"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        input_type: "type of input data: CRAM, BAM, or FASTQ"
        input_cram_bam: "input CRAM or BAM file (required when input_type is CRAM or BAM)"
        input_crai_bai: "index of input_cram_bam"
        fastqs_pe_forward: "forward (R1) FASTQ files for paired-end input (required when input_type is FASTQ)"
        fastqs_pe_reverse: "reverse (R2) FASTQ files for paired-end input, in matching order with fastqs_pe_forward"
        fastqs_pe_readgroup: "@RG header lines for each paired-end FASTQ pair"
        fastqs_pe_readgroup_id: "read group IDs for each paired-end FASTQ pair"
        fastqs_se: "single-end FASTQ files for FASTQ input"
        fastqs_se_readgroup: "@RG header lines for each single-end FASTQ"
        fastqs_se_readgroup_id: "read group IDs for each single-end FASTQ"
        delivery_ref_fasta: "reference FASTA the input CRAM/BAM was aligned to, used to decode it (required when input_type is CRAM)"
        delivery_ref_fasta_index: "index of delivery_ref_fasta"
        gb_per_bwa_bam: "target GiB of unmapped BAM per bwa shard, used to size the within-read-group FASTQ split"
        run_bqsr: "whether to run base quality score recalibration"
        rerun_bqsr: "whether to re-run BQSR even when the input was already recalibrated (requires OQ tags)"
        dbsnp_vcf: "dbSNP VCF of known sites for BQSR"
        dbsnp_vcf_index: "index of dbsnp_vcf"
        ref_fasta: "reference FASTA to align to"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary of ref_fasta"
        ref_amb: "bwa .amb index file"
        ref_ann: "bwa .ann index file"
        ref_bwt: "bwa .bwt index file"
        ref_pac: "bwa .pac index file"
        ref_sa: "bwa .sa index file"
        ref_alt: "bwa .alt file of alternate contigs (optional)"

        # outputs
        analysis_ready_bam: "analysis-ready aligned, duplicate-marked, and optionally recalibrated BAM"
        analysis_ready_bai: "index of analysis_ready_bam"
    }

    input {
        String? wandb_api_key

        String sample_id
        String input_type # "CRAM", "BAM", or "FASTQ"
        File? input_cram_bam
        File? input_crai_bai
        Array[File]? fastqs_pe_forward
        Array[File]? fastqs_pe_reverse
        Array[String]? fastqs_pe_readgroup
        Array[String]? fastqs_pe_readgroup_id
        Array[File]? fastqs_se
        Array[String]? fastqs_se_readgroup
        Array[String]? fastqs_se_readgroup_id
        File? delivery_ref_fasta
        File? delivery_ref_fasta_index

        Int gb_per_bwa_bam = 4

        Boolean run_bqsr = true
        Boolean rerun_bqsr = true
        File dbsnp_vcf
        File dbsnp_vcf_index

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File? ref_alt
    }

    if (input_type == "BAM" || input_type == "CRAM") {
        if (input_type == "CRAM") {
            File cram_ref_fasta = select_first([delivery_ref_fasta])
            File cram_ref_fasta_index = select_first([delivery_ref_fasta_index])
        }

        call split_bam_by_read_group {
            input:
                wandb_api_key = wandb_api_key,
                input_cram_bam = select_first([input_cram_bam]),
                ref_fasta = cram_ref_fasta,
                ref_fasta_index = cram_ref_fasta_index
        }

        scatter (ubam in split_bam_by_read_group.output_bams) {
            call revert_and_validate_sam {
                input:
                    wandb_api_key = wandb_api_key,
                    input_bam = ubam,
                    output_bam_filename = basename(ubam),
                    ref_fasta = select_first([delivery_ref_fasta])
            }
        }

        scatter (ubam in revert_and_validate_sam.output_bam) {
            Int n_shards = ceil(size(ubam, "GiB") / gb_per_bwa_bam)

            call get_rg {
                input:
                    wandb_api_key = wandb_api_key,
                    bam = ubam
            }

            call biobambam_bamtofastq {
                input:
                    wandb_api_key = wandb_api_key,
                    filename = ubam,
                    n_shards = n_shards
            }

            Array[File] fastq1 = biobambam_bamtofastq.output_fastq1
            Array[File] fastq2 = biobambam_bamtofastq.output_fastq2
            Array[File] fastq_o1 = biobambam_bamtofastq.output_fastq_o1
            Array[File] fastq_o2 = biobambam_bamtofastq.output_fastq_o2
            Array[File] fastq_s = biobambam_bamtofastq.output_fastq_s

            Int pe_count = length(fastq1)
            Int o1_count = length(fastq_o1)
            Int o2_count = length(fastq_o2)
            Int s_count = length(fastq_s)

            if (pe_count > 0) {
                scatter (pe_record in zip(fastq1, fastq2)) {
                    call bwa_pe {
                        input:
                            wandb_api_key = wandb_api_key,
                            fastq1 = pe_record.left,
                            fastq2 = pe_record.right,
                            rg_header = get_rg.rg_header,
                            ref_fasta = ref_fasta,
                            ref_fasta_index = ref_fasta_index,
                            ref_dict = ref_dict,
                            ref_amb = ref_amb,
                            ref_ann = ref_ann,
                            ref_bwt = ref_bwt,
                            ref_pac = ref_pac,
                            ref_sa = ref_sa,
                            ref_alt = ref_alt
                    }
                }
            }

            if (o1_count + o2_count + s_count > 0) {
                Array[File] fastq_single_records = flatten([fastq_o1, fastq_o2, fastq_s])

                scatter (fastq in fastq_single_records) {
                    call bwa_se {
                        input:
                            wandb_api_key = wandb_api_key,
                            fastq = fastq,
                            rg_header = get_rg.rg_header,
                            ref_fasta = ref_fasta,
                            ref_fasta_index = ref_fasta_index,
                            ref_dict = ref_dict,
                            ref_amb = ref_amb,
                            ref_ann = ref_ann,
                            ref_bwt = ref_bwt,
                            ref_pac = ref_pac,
                            ref_sa = ref_sa,
                            ref_alt = ref_alt
                    }
                }
            }
        }

        Array[File] pe_aligned_bams = flatten(select_all(bwa_pe.bam))
        Array[File] se_aligned_bams = flatten(select_all(bwa_se.bam))
        Array[File] rg_aligned_bams = flatten([pe_aligned_bams, se_aligned_bams])
    }


    if (input_type == "FASTQ") {
        Int fastqs_pe_count = length(select_first([fastqs_pe_forward, []]))

        if (fastqs_pe_count > 0) {
            Array[File] pe_forward = select_first([fastqs_pe_forward])
            Array[File] pe_reverse = select_first([fastqs_pe_reverse])
            Array[String] pe_readgroup = select_first([fastqs_pe_readgroup])
            Array[String] pe_readgroup_id = select_first([fastqs_pe_readgroup_id])

            scatter (i in range(fastqs_pe_count)) {
                call bwa_pe as fastq_bwa_pe {
                    input:
                        wandb_api_key = wandb_api_key,
                        fastq1 = pe_forward[i],
                        fastq2 = pe_reverse[i],
                        rg_header = pe_readgroup[i],
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_amb = ref_amb,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_sa = ref_sa,
                        ref_alt = ref_alt
                }
            }
        }

        Int fastqs_se_count = length(select_first([fastqs_se, []]))

        if (fastqs_se_count > 0) {
            Array[File] se_fastqs = select_first([fastqs_se])
            Array[String] se_readgroup = select_first([fastqs_se_readgroup])
            Array[String] se_readgroup_id = select_first([fastqs_se_readgroup_id])

            scatter (i in range(fastqs_se_count)) {
                call bwa_se as fastq_bwa_se {
                    input:
                        wandb_api_key = wandb_api_key,
                        fastq = se_fastqs[i],
                        rg_header = se_readgroup[i],
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        ref_amb = ref_amb,
                        ref_ann = ref_ann,
                        ref_bwt = ref_bwt,
                        ref_pac = ref_pac,
                        ref_sa = ref_sa,
                        ref_alt = ref_alt
                }
            }
        }
    }

    # Require at least one aligned BAM; the `+` makes miniwdl assert nonemptiness
    # at runtime. The coercion can't be proven nonempty statically, so suppress lint.
    Array[File]+ aligned_bams = flatten([
        select_first([rg_aligned_bams, []]),
        select_first([fastq_bwa_pe.bam, []]),
        select_first([fastq_bwa_se.bam, []]),
    ])  # !NonemptyCoercion

    call picard_markduplicates {
        input:
            wandb_api_key = wandb_api_key,
            bams = aligned_bams,
            outbam = sample_id + ".mrkdup.bam"
    }

    call sort_and_index_markdup_bam {
        input:
            wandb_api_key = wandb_api_key,
            input_bam = picard_markduplicates.bam,
            output_bam_basename = sample_id + ".analysis_ready"
    }

    if (run_bqsr) {
        if (input_type == "BAM" || input_type == "CRAM") {
            call check_bqsr_and_oq_tags {
                input:
                    wandb_api_key = wandb_api_key,
                    bam = select_first([input_cram_bam]),
                    bai = select_first([input_crai_bai]),
                    ref_fasta = ref_fasta
            }
        }

        if (input_type == "FASTQ" || !select_first([check_bqsr_and_oq_tags.bqsr_performed, false]) || (rerun_bqsr && select_first([check_bqsr_and_oq_tags.has_oq_tags, false]))) {
            call create_sequence_grouping_tsv {
                input:
                    wandb_api_key = wandb_api_key,
                    ref_dict = ref_dict
            }

            Int n_bqsr_splits = length(create_sequence_grouping_tsv.sequence_grouping_with_unmapped)

            scatter (subgroup in create_sequence_grouping_tsv.sequence_grouping) {
                call base_recalibrator {
                    input:
                        wandb_api_key = wandb_api_key,
                        bqsr_regions = subgroup,
                        n_bqsr_splits = n_bqsr_splits,
                        bam = sort_and_index_markdup_bam.output_bam,
                        bam_index = sort_and_index_markdup_bam.output_bai,
                        ref_fasta = ref_fasta,
                        ref_fasta_index = ref_fasta_index,
                        ref_dict = ref_dict,
                        dbsnp_vcf = dbsnp_vcf,
                        dbsnp_vcf_index = dbsnp_vcf_index
                }
            }

            call gather_bqsr_reports {
                input:
                    wandb_api_key = wandb_api_key,
                    input_bqsr_reports = base_recalibrator.bqsr_recal_file,
                    output_report_filename = sample_id + ".bqsr.grp"
            }

            scatter (subgroup in create_sequence_grouping_tsv.sequence_grouping_with_unmapped) {
                call apply_bqsr {
                    input:
                        wandb_api_key = wandb_api_key,
                        bqsr_regions = subgroup,
                        n_bqsr_splits = n_bqsr_splits,
                        bam = sort_and_index_markdup_bam.output_bam,
                        bam_index = sort_and_index_markdup_bam.output_bai,
                        bqsr_recal_file = gather_bqsr_reports.output_bqsr_report
                }
            }

            call gather_bam_files as bqsr_gather_bam_files {
                input:
                    wandb_api_key = wandb_api_key,
                    input_bams = apply_bqsr.recalibrated_bam,
                    output_bam_basename = sample_id + ".analysis_ready"
            }
        }
    }

    File final_bam = select_first([bqsr_gather_bam_files.output_bam, sort_and_index_markdup_bam.output_bam])
    File final_bai = select_first([bqsr_gather_bam_files.output_bai, sort_and_index_markdup_bam.output_bai])

    output {
        File analysis_ready_bam = final_bam
        File analysis_ready_bai = final_bai
    }
}

task split_bam_by_read_group {
    meta {
        description: "Split a CRAM/BAM into one BAM per read group"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_cram_bam: "input CRAM or BAM file to split by read group"
        ref_fasta: "reference FASTA needed to decode a CRAM input"
        ref_fasta_index: "index of ref_fasta"

        # outputs
        output_bams: "one BAM per read group found in the input"
    }

    input {
        String? wandb_api_key

        File input_cram_bam
        File? ref_fasta
        File? ref_fasta_index

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 4
        Int mem_gb = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_gb = ceil(size(input_cram_bam, "GiB") * 10) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        mkdir -p split

        samtools split \
            -@ ~{cpu} \
            -f 'split/%!.bam' \
            ~{"--reference " + ref_fasta} \
            ~{input_cram_bam}
    >>>

    output {
        Array[File] output_bams = glob("split/*.bam")
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }
}

task revert_and_validate_sam {
    meta {
        description: "Revert an aligned BAM to unmapped, coordinate-sorted reads and validate resulting BAM"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_bam: "input aligned BAM to revert"
        ref_fasta: "reference FASTA"
        output_bam_filename: "filename for the reverted BAM"

        # outputs
        output_bam: "reverted, validated BAM"
    }

    input {
        String? wandb_api_key

        File input_bam
        File ref_fasta
        String output_bam_filename

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/picard"
        String docker_image_hash_or_tag = ":3.4.0"
        Int mem_mb = 8192
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_gb = ceil(size(input_bam, "GiB") * 6 + size(ref_fasta, "GiB")) + 20 + additional_disk_gb
    Int java_mem_mb = mem_mb - 1000
    Int max_heap_mb = mem_mb - 500

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        java -Xms~{java_mem_mb}m -Xmx~{max_heap_mb}m -jar /usr/picard/picard.jar \
            RevertSam \
            --INPUT ~{input_bam} \
            --OUTPUT ~{output_bam_filename} \
            --VALIDATION_STRINGENCY LENIENT \
            --ATTRIBUTE_TO_CLEAR FT \
            --ATTRIBUTE_TO_CLEAR CO \
            --ATTRIBUTE_TO_CLEAR PA \
            --ATTRIBUTE_TO_CLEAR OA \
            --ATTRIBUTE_TO_CLEAR XA \
            --RESTORE_HARDCLIPS false \
            --SORT_ORDER coordinate

        java -Xms~{java_mem_mb}m -Xmx~{max_heap_mb}m -jar /usr/picard/picard.jar \
            ValidateSamFile \
            --INPUT ~{output_bam_filename} \
            --OUTPUT validation_report.txt \
            --MODE SUMMARY \
            --IS_BISULFITE_SEQUENCED false

        grep -q "No errors found" validation_report.txt || {
            echo "VALIDATE_SAM failed:"
            cat validation_report.txt
            exit 1
        }
    >>>

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_mb} MB"
        preemptible: preemptible
    }

    output {
        File output_bam = output_bam_filename
    }
}

task get_rg {
    meta {
        description: "Read the single @RG header line from a read-group-split BAM and extract its read group ID"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        bam: { description: "read-group-split BAM (expected to contain exactly one @RG)", localization_optional: true }

        # outputs
        rg_header: "the @RG header line"
        rg_id: "the read group ID parsed from the @RG line"
    }

    input {
        String? wandb_api_key

        File bam

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int disk_gb = 12
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
    }

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        # Set up auth for samtools streaming from GCS buckets. When running locally
        # (e.g. a plain Docker run with local file inputs), gcloud may be present but
        # unauthenticated; in that case leave these unset and let samtools read the
        # localized files directly rather than failing under `set -e`.
        if token="$(gcloud auth application-default print-access-token 2>/dev/null)"; then
            export GCS_OAUTH_TOKEN="$token"
        fi
        if project="$(gcloud config get-value project -q 2>/dev/null)"; then
            export GCS_REQUESTER_PAYS_PROJECT="$project"
        fi

        # extract the @RG header line
        samtools view -H ~{bam} | grep '^@RG' > rg.txt

        # since we already split by read group, there should be only one
        [ "$(wc -l < rg.txt)" -eq 1 ] || { echo "multiple @RG values in RG-split BAM"; exit 1; }

        # extract the ID of this read group
        awk -F'\t' '
            {
                for (i = 1; i <= NF; i++) {
                    if ($i ~ /^ID:/) {
                        sub(/^ID:/, "", $i)
                        print $i
                        exit
                    }
                }
            }' rg.txt > rg_id.txt
    >>>

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }

    output {
        String rg_header = read_string("rg.txt")
        String rg_id = read_string("rg_id.txt")
    }
}

task biobambam_bamtofastq {
    meta {
        description: "Convert a single-read-group unmapped BAM to FASTQ and shard the paired-end FASTQ file pairs"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        exclude: "comma-separated read flags to exclude from output"
        filename: "input unmapped BAM to convert to FASTQ"
        level: "gzip compression level for output FASTQ"
        tryoq: "if 1, use OQ tag quality scores when present"
        n_shards: "number of chunks to split each read group's paired reads into"

        # outputs
        output_fastq1: "forward (R1) FASTQ chunks for paired reads"
        output_fastq2: "reverse (R2) FASTQ chunks for paired reads, in matching order with output_fastq1"
        output_fastq_o1: "orphaned first-of-pair reads"
        output_fastq_o2: "orphaned second-of-pair reads"
        output_fastq_s: "single-end reads"
    }

    input {
        String? wandb_api_key

        String exclude = "QCFAIL,SECONDARY,SUPPLEMENTARY"
        File filename
        Int level = 5
        Int tryoq = 1
        Int n_shards

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/biobambam2"
        String docker_image_hash_or_tag = ":2.0.87-release-20180301132713"
        Int cpu = 2
        Int mem_gb = 8
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_gb = ceil(size(filename, "GiB") * 4) + 20 + additional_disk_gb

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        echo "Converting to FASTQ"
        /usr/local/bin/bamtofastq \
            T=tempfq \
            collate=1 \
            exclude=~{exclude} \
            filename=~{filename} \
            gz=1 \
            inputformat="bam" \
            level=~{level} \
            outputdir=. \
            outputperreadgroup=1 \
            outputperreadgroupsuffixF="_1.fq.gz" \
            outputperreadgroupsuffixF2="_2.fq.gz" \
            outputperreadgroupsuffixO="_o1.fq.gz" \
            outputperreadgroupsuffixO2="_o2.fq.gz" \
            outputperreadgroupsuffixS="_s.fq.gz" \
            tryoq=~{tryoq}

        n_shards=~{n_shards}

        for r1 in *_1.fq.gz; do
            [ -f "$r1" ] || continue
            rgid="${r1%_1.fq.gz}"
            r2="${rgid}_2.fq.gz"

            echo "Splitting $rgid"

            total_lines=$(zcat "$r1" | wc -l)
            chunk_reads=$(( (total_lines / 4 + n_shards - 1) / n_shards ))
            chunk_lines=$(( chunk_reads * 4 ))
            [ "$chunk_lines" -gt 0 ] || chunk_lines=4

            echo "Using $chunk_lines lines per chunk"

            # $FILE is expanded by split's --filter subshell, not by this shell,
            # so single quotes are intentional here.
            # shellcheck disable=SC2016
            zcat "$r1" | split -l "$chunk_lines" -d -a 3 \
                --additional-suffix=_1.fq.gz \
                --filter='gzip -5 > "$FILE"' \
                - "${rgid}__chunk__"
            rm "$r1"

            # shellcheck disable=SC2016
            zcat "$r2" | split -l "$chunk_lines" -d -a 3 \
                --additional-suffix=_2.fq.gz \
                --filter='gzip -5 > "$FILE"' \
                - "${rgid}__chunk__"
            rm "$r2"
        done

        echo "FASTQ files created:"
        ls -l ./*.fq.gz

        # Fail loudly if any output FASTQ is a truncated/incomplete gzip stream
        # rather than letting bwa emit a malformed SAM record downstream.
        echo "Checking .fq.gz files"
        for fq in *.fq.gz; do
            [ -e "$fq" ] || continue
            gzip -t "$fq" || { echo "ERROR: ${fq} is a truncated/corrupt gzip" 1>&2; exit 1; }
        done
    >>>

    output {
        Array[File] output_fastq1 = glob("*__chunk__*_1.fq.gz")
        Array[File] output_fastq2 = glob("*__chunk__*_2.fq.gz")
        Array[File] output_fastq_o1 = glob("*_o1.fq.gz")
        Array[File] output_fastq_o2 = glob("*_o2.fq.gz")
        Array[File] output_fastq_s = glob("*_s.fq.gz")
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }
}

task bwa_pe {
    meta {
        description: "Align a paired-end FASTQ shard with bwa mem and write an unsorted BAM"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        fastq1: "forward (R1) FASTQ"
        fastq2: "reverse (R2) FASTQ, in matching order with fastq1"
        rg_header: "@RG header line to attach to aligned reads"
        ref_fasta: "reference FASTA to align to"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary of ref_fasta"
        ref_amb: "bwa .amb index file"
        ref_ann: "bwa .ann index file"
        ref_bwt: "bwa .bwt index file"
        ref_pac: "bwa .pac index file"
        ref_sa: "bwa .sa index file"
        ref_alt: "bwa .alt file of alternate contigs (optional)"
        bwa_batch_size: "bwa mem -K: input bases per batch; fixing it makes output thread-count-independent and bounds peak memory"

        # outputs
        bam: "unsorted aligned BAM"
    }

    input {
        String? wandb_api_key

        File fastq1
        File fastq2
        String rg_header
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File? ref_alt # !UnusedDeclaration

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bwa"
        String docker_image_hash_or_tag = ":0.7.15-r1142-dirty"
        Int cpu = 16
        Int preemptible = 1
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0

        # bwa mem -K: number of input bases processed per batch. Fixing this makes
        # output independent of thread count (reproducibility). It is also the main
        # driver of variable peak RAM: a batch of this many bases plus its alignment
        # results is held in memory at once. Lowering it (e.g. 50000000) cuts peak
        # memory for very large shards at a small throughput/reproducibility cost.
        Int bwa_batch_size = 100000000
    }

    Float ref_size = size([
        ref_fasta,
        ref_dict,
        ref_amb,
        ref_ann,
        ref_bwt,
        ref_pac,
        ref_sa,
        ref_fasta_index,
        ref_alt
    ], "GiB")

    Int disk_gb = ceil((size([fastq1, fastq2], "GiB") * 10) + ref_size) + 20 + additional_disk_gb

    # Peak RSS ~= the BWA index resident in RAM (~6 GiB for hg38) + the per-batch
    # working set + samtools view compression buffers + OS headroom. The variable
    # part scales with thread count (each thread holds part of the in-flight batch
    # and its alignment results), so it is modeled as cpu * 750 MiB rather than purely
    # from FASTQ size. A size-proportional term is kept as empirical headroom since
    # large shards (each FASTQ 15-30 GiB) were OOMing under the previous 24 GiB cap.
    # Floor at 24000 MiB (sane minimum machine) and cap at 49152 MiB (48 GiB) so a
    # 30+30 GiB shard requests ~42 GiB instead of being clamped to 24 GiB.
    Int computed_mem_mb = 12000 + (cpu * 750) + ceil(size([fastq1, fastq2], "MiB") * 0.3) + additional_memory_mb
    Int mem_mb = if computed_mem_mb < 24000 then 24000 else if computed_mem_mb > 49152 then 49152 else computed_mem_mb

    String out_basename = sub(basename(fastq1, ".fq.gz"), "_1$", "")

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        bwa mem -K ~{bwa_batch_size} \
            -t ~{cpu} \
            -T 0 \
            -R "~{rg_header}" \
            ~{ref_fasta} \
            ~{fastq1} \
            ~{fastq2} \
        | samtools view \
            -@ ~{cpu} \
            -Shb \
            -o "~{out_basename}.bam" \
            -
    >>>

    output {
        File bam = "~{out_basename}.bam"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_mb} MiB"
        preemptible: preemptible
    }
}

task bwa_se {
    meta {
        description: "Align a single-end FASTQ shard with bwa mem and write an unsorted BAM"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        fastq: "single-end FASTQ to align"
        rg_header: "@RG header line to attach to aligned reads"
        ref_fasta: "reference FASTA to align to"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary of ref_fasta"
        ref_amb: "bwa .amb index file"
        ref_ann: "bwa .ann index file"
        ref_bwt: "bwa .bwt index file"
        ref_pac: "bwa .pac index file"
        ref_sa: "bwa .sa index file"
        ref_alt: "bwa .alt file of alternate contigs (optional)"
        bwa_batch_size: "bwa mem -K: input bases per batch; fixing it makes output thread-count-independent and bounds peak memory"

        # outputs
        bam: "unsorted aligned BAM"
    }

    input {
        String? wandb_api_key

        File fastq
        String rg_header
        File ref_fasta
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fasta_index
        File? ref_alt

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bwa"
        String docker_image_hash_or_tag = ":0.7.15-r1142-dirty"
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0

        # bwa mem -K: number of input bases processed per batch. Fixing this makes
        # output independent of thread count (reproducibility). It is also the main
        # driver of variable peak RAM: a batch of this many bases plus its alignment
        # results is held in memory at once. Lowering it (e.g. 50000000) cuts peak
        # memory for very large shards at a small throughput/reproducibility cost.
        Int bwa_batch_size = 100000000
    }

    Float ref_size = size([
        ref_fasta,
        ref_dict,
        ref_amb,
        ref_ann,
        ref_bwt,
        ref_pac,
        ref_sa,
        ref_fasta_index,
        ref_alt
    ], "GiB")

    Int disk_gb = ceil((size(fastq, "GiB") * 4) + ref_size) + 20 + additional_disk_gb

    # Same model as bwa_pe: BWA index (~6 GiB for hg38) + a per-batch working set that
    # scales with thread count (cpu * 750 MiB) + a size-proportional headroom term.
    # Floor/cap are lower than bwa_pe since SE handles few, small singleton reads, but
    # the cap still allows headroom for the rare large shard.
    Int computed_mem_mb = 12000 + (cpu * 750) + ceil(size(fastq, "MiB") * 0.3) + additional_memory_mb
    Int mem_mb = if computed_mem_mb < 16000 then 16000 else if computed_mem_mb > 32768 then 32768 else computed_mem_mb

    String out_basename = basename(fastq, ".fq.gz")

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        bwa mem -K ~{bwa_batch_size} \
            -t ~{cpu} \
            -T 0 \
            -R "~{rg_header}" \
            ~{ref_fasta} \
            ~{fastq} \
        | samtools view \
            -@ ~{cpu} \
            -Shb \
            -o "~{out_basename}.bam" \
            -
    >>>

    output {
        File bam = "~{out_basename}.bam"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_mb} MiB"
        preemptible: preemptible
    }
}

task picard_markduplicates {
    meta {
        description: "Mark duplicate reads across the per-shard aligned BAMs with Picard MarkDuplicates, merging them into one BAM"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        bams: "aligned BAM shards to merge and mark duplicates over"
        outbam: "filename for the merged, duplicate-marked BAM"
        compression_level: "compression level for the output BAM"
        validation_stringency: "Picard validation stringency"
        assume_sort_order: "sort order Picard should assume for the inputs"
        read_name_regex: "regex for parsing read names to detect optical duplicates; set to 'null' to disable optical duplicate detection"
        sorting_collection_size_ratio: "Picard SORTING_COLLECTION_SIZE_RATIO override"

        # outputs
        bam: "merged, duplicate-marked BAM"
        metrics: "duplication metrics file"
    }

    input {
        String? wandb_api_key

        Array[File]+ bams
        String outbam

        Int compression_level = 2
        String validation_stringency = "SILENT"
        String assume_sort_order = "queryname"

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
        # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
        String? read_name_regex
        Int additional_disk_gb = 20

        Float? sorting_collection_size_ratio

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/picard"
        String docker_image_hash_or_tag = ":3.4.0"
        Int mem_gb = 8
        Int cpu = 1
        Int preemptible = 1
        Int max_retries = 1
    }

    Int disk_gb = 3 * ceil(size(bams, "GiB")) + 20 + additional_disk_gb

    # Leave ~2 GiB of headroom for JVM off-heap (metaspace, thread stacks, GC),
    # htsjdk I/O buffers, and the OS so the container isn't OOM-killed. Xms == Xmx
    # for a stable, predictable heap footprint.
    Int java_heap_gb = mem_gb - 2

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        java \
            -Xms~{java_heap_gb}g -Xmx~{java_heap_gb}g \
            -Dsamjdk.compression_level=~{compression_level} \
            -jar /usr/picard/picard.jar \
            MarkDuplicates \
            INPUT=~{sep=' INPUT=' bams} \
            OUTPUT=~{outbam} \
            METRICS_FILE="~{outbam}.metrics" \
            VALIDATION_STRINGENCY=~{validation_stringency} \
            ASSUME_SORT_ORDER=~{assume_sort_order} \
            ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
            ~{"READ_NAME_REGEX=" + read_name_regex}
    >>>

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GiB"
        preemptible: preemptible
    }

    output {
        File bam = "~{outbam}"
        File metrics = "~{outbam}.metrics"
    }
}

task sort_and_index_markdup_bam {
    meta {
        description: "Coordinate-sort and index a duplicate-marked BAM with samtools"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_bam: "duplicate-marked BAM to sort"
        output_bam_basename: "basename for the sorted BAM and its index"
        tmp_prefix: "prefix for samtools sort temporary files"

        # outputs
        output_bam: "coordinate-sorted BAM"
        output_bai: "index of output_bam"
    }

    input {
        String? wandb_api_key

        File input_bam
        String output_bam_basename
        String tmp_prefix = "tmp_srt"

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int computed_mem_mb = ceil(size(input_bam, "MiB") * 0.25) + 12000 + additional_memory_mb
    Int mem_mb = if computed_mem_mb < 48000 then computed_mem_mb else 48000
    Int disk_gb = ceil(size(input_bam, "GiB") * 3) + 10 + additional_disk_gb
    Int mem_mb_per_thread = floor(mem_mb / cpu * 0.85)
    Int index_threads = cpu
    String output_bam_name = output_bam_basename + ".bam"
    String output_bai_name = output_bam_basename + ".bai"

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        samtools sort \
            -@ ~{cpu} \
            -o ~{output_bam_name} \
            -T ~{tmp_prefix} \
            -m ~{mem_mb_per_thread}M \
            ~{input_bam}

        samtools index \
            -@ ~{index_threads} \
            -b \
            ~{output_bam_name} \
            ~{output_bai_name}
    >>>

    output {
        File output_bam = "~{output_bam_name}"
        File output_bai = "~{output_bai_name}"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_mb} MiB"
        preemptible: preemptible
    }
}

task check_bqsr_and_oq_tags {
    meta {
        description: "Inspect a BAM header and a read to determine whether BQSR was already applied and whether reads carry OQ (original quality) tags"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        bam: { description: "BAM to inspect (streamed, not localized)", localization_optional: true }
        bai: { description: "index of bam", localization_optional: true }
        ref_fasta: { description: "reference FASTA", localization_optional: true }

        # outputs
        bqsr_performed: "true if the header shows BQSR was already applied"
        has_oq_tags: "true if reads carry OQ original-quality tags"
    }

    input {
        String? wandb_api_key

        File bam
        File bai
        File ref_fasta

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_gb = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        # Set up auth for samtools streaming from GCS buckets. When running locally
        # (e.g. a plain Docker run with local file inputs), gcloud may be present but
        # unauthenticated; in that case leave these unset and let samtools read the
        # localized files directly rather than failing under `set -e`.
        if token="$(gcloud auth application-default print-access-token 2>/dev/null)"; then
            export GCS_OAUTH_TOKEN="$token"
        fi
        if project="$(gcloud config get-value project -q 2>/dev/null)"; then
            export GCS_REQUESTER_PAYS_PROJECT="$project"
        fi

        # get header (and single read) to check for previous BQSR recalibration
        samtools head "~{bam}" -n 1 --reference "~{ref_fasta}" > head.txt

        BQSR_PERFORMED=$(grep -q "apply_bqsr" head.txt && echo "true" || echo "false") && \
            echo "${BQSR_PERFORMED}" > bqsr_performed.txt

        if [ "${BQSR_PERFORMED}" = "true" ]; then
            echo "BQSR has already been performed on ~{bam}"
        else
            echo "BQSR has not been performed on ~{bam} yet"
        fi

        HAS_OQ_TAGS=$(tail -n 1 head.txt | awk '
            BEGIN {
                FS="\t"
                found = 0
            }

            {
                for (i=1; i<=NF; i++) {
                    if ($i ~ /^OQ:/) {
                        found = 1
                        print "true"
                        exit
                    }
                }
            }

            END {
                if (!found) {
                    print "false"
                }
            }
        ') && echo "${HAS_OQ_TAGS}" > has_oq_tags.txt

        if [ "${HAS_OQ_TAGS}" = "true" ]; then
            echo "Found OQ tags; BQSR can be performed again with --use-original-qualities"
        else
            echo "Reads are missing OQ tag; BQSR cannot be performed again"
        fi
    >>>

    output {
        Boolean bqsr_performed = read_boolean("bqsr_performed.txt")
        Boolean has_oq_tags = read_boolean("has_oq_tags.txt")
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }
}

# BQSR scatter-gather tasks from
# https://github.com/broadinstitute/warp/blob/9942eed7d6c8aba87ddf096fc77cddc6f5052e54/tasks/broad/Utilities.wdl
task create_sequence_grouping_tsv {
    meta {
        description: "Group reference contigs into balanced interval lists for scattering BQSR, from the reference sequence dictionary"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        ref_dict: "sequence dictionary of the reference FASTA"

        # outputs
        sequence_grouping: "balanced groups of contigs for scattering BaseRecalibrator"
        sequence_grouping_with_unmapped: "sequence_grouping plus a final group for unmapped reads, for scattering ApplyBQSR"
    }

    input {
        String? wandb_api_key

        File ref_dict

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/python_extra"
        String docker_image_hash_or_tag = ":3.12.4-slim"
        Int mem_gb = 1
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_gb = 1 + additional_disk_gb

    # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
    # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
    # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        python3 <<CODE
        with open("~{ref_dict}", "r") as ref_dict_file:
            sequence_tuple_list = []
            longest_sequence = 0

            for line in ref_dict_file:
                if line.startswith("@SQ"):
                    line_split = line.split("\t")
                    # (Sequence_Name, Sequence_Length)
                    sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))

            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

        # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
        # the last element after a :, so we add this as a sacrificial element.
        hg38_protection_tag = ":1+"
        # initialize the tsv string with the first sequence
        tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
        temp_size = sequence_tuple_list[0][1]

        for sequence_tuple in sequence_tuple_list[1:]:
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
            else:
                tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
                temp_size = sequence_tuple[1]

        # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
        with open("sequence_grouping.txt","w") as tsv_file:
            tsv_file.write(tsv_string)
            tsv_file.close()

        tsv_string += '\n' + "unmapped"

        with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
            tsv_file_with_unmapped.write(tsv_string)
            tsv_file_with_unmapped.close()
        CODE
    >>>

    output {
        Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
        Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }
}

task base_recalibrator {
    meta {
        description: "Compute a BQSR recalibration table over a group of intervals with GATK BaseRecalibrator"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        bqsr_regions: "intervals this shard recalibrates over"
        n_bqsr_splits: "total number of BQSR shards, used to size disk and memory"
        bam: { description: "sorted, duplicate-marked BAM (streamed, not localized)", localization_optional: true }
        bam_index: { description: "index of bam", localization_optional: true }
        dbsnp_vcf: { description: "dbSNP VCF of known sites", localization_optional: true }
        dbsnp_vcf_index: { description: "index of dbsnp_vcf", localization_optional: true }
        ref_dict: { description: "sequence dictionary of ref_fasta", localization_optional: true }
        ref_fasta: { description: "reference FASTA", localization_optional: true }
        ref_fasta_index: { description: "index of ref_fasta", localization_optional: true }

        # outputs
        bqsr_recal_file: "BQSR recalibration table for this interval group"
    }

    input {
        String? wandb_api_key

        Array[String] bqsr_regions
        Int n_bqsr_splits
        File bam
        File bam_index
        File dbsnp_vcf
        File dbsnp_vcf_index
        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk"
        String docker_image_hash_or_tag = ":4.6.1.0"
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    String output_grp = basename(bam, ".bam") + "_bqsr.grp"

    Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Float dbsnp_size = size(dbsnp_vcf, "GiB")
    Int disk_gb = ceil(
        (size(bam, "GiB") / n_bqsr_splits) + ref_size + dbsnp_size
    ) + additional_disk_gb + 20

    Int mem_mb = ceil(size(bam, "MiB") / n_bqsr_splits) + 6000 + additional_memory_mb
    Int jvm_mem_mb = mem_mb - 1000
    Int max_heap_mb = mem_mb - 500

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        java \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -XX:+PrintFlagsFinal \
            -Xlog:gc=debug:file=gc_log.log \
            -Xms~{jvm_mem_mb}m -Xmx~{max_heap_mb}m \
            -jar /root/gatk.jar \
            BaseRecalibrator \
            --input ~{bam} \
            --use-original-qualities \
            --known-sites ~{dbsnp_vcf} \
            --reference ~{ref_fasta} \
            --tmp-dir . \
            --output ~{output_grp} \
            --intervals ~{sep=" --intervals " bqsr_regions}
    >>>

    output {
        File bqsr_recal_file = "~{output_grp}"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_mb} MiB"
        preemptible: preemptible
    }
}

task gather_bqsr_reports {
    meta {
        description: "Merge per-interval BQSR recalibration tables into one with GATK GatherBQSRReports"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_bqsr_reports: "per-shard BQSR recalibration tables to merge"
        output_report_filename: "filename for the merged recalibration table"

        # outputs
        output_bqsr_report: "merged BQSR recalibration table"
    }

    input {
        String? wandb_api_key

        Array[File] input_bqsr_reports
        String output_report_filename

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk"
        String docker_image_hash_or_tag = ":4.6.1.0"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_gb = 3 * ceil(size(input_bqsr_reports, "GiB")) + additional_disk_gb

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        java -Xms3000m -Xmx3000m -jar /root/gatk.jar \
            GatherBQSRReports \
            -I ~{sep=' -I ' input_bqsr_reports} \
            -O ~{output_report_filename}
    >>>

    output {
        File output_bqsr_report = "~{output_report_filename}"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }
}

task apply_bqsr {
    meta {
        description: "Apply a BQSR recalibration table to a group of intervals with GATK ApplyBQSR"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        bqsr_regions: "intervals this shard recalibrates over"
        n_bqsr_splits: "total number of BQSR shards, used to size disk and memory"
        bam: { description: "sorted, duplicate-marked BAM (streamed, not localized)", localization_optional: true }
        bam_index: { description: "index of bam", localization_optional: true }
        bqsr_recal_file: "merged BQSR recalibration table to apply"
        emit_original_quals: "whether to write original qualities to OQ tags in the output"

        # outputs
        recalibrated_bam: "recalibrated BAM for this interval group"
    }

    input {
        String? wandb_api_key

        Array[String] bqsr_regions
        Int n_bqsr_splits
        File bam
        File bam_index
        File bqsr_recal_file
        Boolean emit_original_quals = false

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/gatk"
        String docker_image_hash_or_tag = ":4.6.1.0"
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }

    Int mem_mb = ceil(size(bam, "MiB") / n_bqsr_splits) + 4000 + additional_memory_mb
    Int jvm_mem_mb = mem_mb - 1000
    Int max_heap_mb = mem_mb - 500
    Int disk_gb = ceil((size(bam, "GiB") * 3) / n_bqsr_splits) + 20 + additional_disk_gb

    String output_bam = basename(bam)

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        java \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -XX:+PrintFlagsFinal \
            -Xlog:gc=debug:file=gc_log.log \
            -Xms~{jvm_mem_mb}m -Xmx~{max_heap_mb}m \
            -jar /root/gatk.jar ApplyBQSR \
            --input ~{bam} \
            --bqsr-recal-file ~{bqsr_recal_file} \
            --emit-original-quals ~{emit_original_quals} \
            --use-original-qualities \
            --add-output-sam-program-record \
            --tmp-dir . \
            --output ~{output_bam} \
            --intervals ~{sep=" --intervals " bqsr_regions}
    >>>

    output {
        File recalibrated_bam = "~{output_bam}"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_mb} MiB"
        preemptible: preemptible
    }
}

task gather_bam_files {
    meta {
        description: "Concatenate per-interval recalibrated BAMs into one indexed BAM with Picard GatherBamFiles"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        input_bams: "per-shard recalibrated BAMs to concatenate, in order"
        output_bam_basename: "basename for the gathered BAM and its index"

        # outputs
        output_bam: "concatenated recalibrated BAM"
        output_bai: "index of output_bam"
    }

    input {
        String? wandb_api_key

        Array[File] input_bams
        String output_bam_basename

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/picard"
        String docker_image_hash_or_tag = ":3.4.0"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int java_max_memory_mb = ceil(mem_gb * 1000) - 500
    Int java_inital_memory_mb = (mem_gb * 1000) - 1000
    Int disk_gb = 2 * ceil(size(input_bams, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        ~{"export WANDB_API_KEY=" + wandb_api_key}
        # shellcheck disable=SC1090
        source <(curl -fsSL https://raw.githubusercontent.com/broadinstitute/magicwand/refs/heads/github-override/install.sh)
        magicwand init

        java \
            -Xms~{java_inital_memory_mb}m \
            -Xmx~{java_max_memory_mb}m \
            -jar /usr/picard/picard.jar \
            GatherBamFiles \
            INPUT=~{sep=' INPUT=' input_bams} \
            OUTPUT=~{output_bam_basename}.bam \
            CREATE_INDEX=true
    >>>

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bai = "~{output_bam_basename}.bai"
    }

    runtime {
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        maxRetries: max_retries
        memory: "~{mem_gb} GB"
        preemptible: preemptible
    }
}
