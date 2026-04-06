version 1.0

workflow prep_annotations {
    meta {
        description: "Prepare a somatic variant VCF for annotation"
    }

    parameter_meta {
        # inputs
        ref_fasta: "reference FASTA"
        ref_fasta_index: "index of ref_fasta"
        ref_dict: "sequence dictionary for ref_fasta"
        sample_id: "ID of this sample"
        input_vcf: "input VCF to prepare; will be bgzip-compressed if not already"
        vcf_patch: "optional VCF of curated variants to substitute into patch_region_bed regions"
        patch_region_bed: "optional BED file of regions whose variants will be replaced by vcf_patch"
        fix_ploidy: "if true, correct multi-allelic genotypes where all reads support the alt "
        filter_vcf: "if true, apply bcftools exclude filter using exclude_string "
        normalize_indels: "if true, split multi-allelic sites and left-align indels "
        exclude_string: "bcftools exclude expression used when filter_vcf is true"

        # outputs
        mut_annot_bcftools_fixed_vcf: "bgzip-compressed VCF after all fixes and filters have been applied"
    }

    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String sample_id
        File input_vcf
        File? vcf_patch
        File? patch_region_bed

        # bcftools-based fixes and filters
        Boolean fix_ploidy = true
        Boolean filter_vcf = true
        Boolean normalize_indels = true
        String? exclude_string
    }

    Boolean needs_compression = sub(
        basename(input_vcf), "\\.vcf\\.gz$", ""
    ) == basename(input_vcf)

    if (needs_compression) {
        call compress_vcf {
            input:
                vcf = input_vcf,
                output_file_base_name = sample_id
        }
    }

    if (defined(vcf_patch) && defined(patch_region_bed)) {
        call patch_vcf {
            input:
                sample_id = sample_id,
                vcf = select_first([compress_vcf.vcf_compressed, input_vcf]),
                vcf_patch = select_first([vcf_patch]),
                patch_region_bed = select_first([patch_region_bed]),
                output_file_base_name = sample_id + "_patched"
        }
    }

    call fix_with_bcftools {
        input:
            vcf = select_first([
                patch_vcf.vcf_patched, compress_vcf.vcf_compressed, input_vcf
            ]),
            output_file_base_name = sample_id,
            fix_ploidy = fix_ploidy,
            filter_vcf = filter_vcf,
            exclude_string = exclude_string,
            normalize_indels = normalize_indels,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File mut_annot_bcftools_fixed_vcf = fix_with_bcftools.vcf_fixed
    }
}

task compress_vcf {
    meta {
        description: "bgzip-compress a VCF using bcftools"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        vcf: "input VCF to compress"
        output_file_base_name: "base name for the output file (without extension)"

        # outputs
        vcf_compressed: "bgzip-compressed output VCF"
    }

    input {
        File vcf
        String output_file_base_name

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(2 * size(vcf, "GiB")) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        bcftools view \
            "~{vcf}" \
            --output="~{output_file_base_name}.vcf.gz" \
            --no-version
    >>>

    output {
        File vcf_compressed = "~{output_file_base_name}.vcf.gz"
    }

    runtime {
        docker: "~{docker_image}~{docker_image_hash_or_tag}"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
        continueOnReturnCode: [0, 9] # FilterMutectCalls generates weird vcf.gz files
    }
}

task patch_vcf {
    meta {
        description: "Replace variants in specified genomic regions using a patch VCF"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        sample_id: "ID of this sample"
        vcf: "base VCF whose variants in patch_region_bed will be replaced"
        vcf_patch: "VCF of curated variants to substitute into patch_region_bed regions"
        patch_region_bed: "BED file of regions whose variants will be replaced"
        output_file_base_name: "base name for the output file (without extension)"

        # outputs
        vcf_patched: "bgzip-compressed VCF with patch regions replaced"
    }

    input {
        String sample_id
        File vcf
        File vcf_patch
        File patch_region_bed
        String output_file_base_name

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = ceil(
        2 * (size(vcf, "GiB") + size(vcf_patch, "GiB"))
    ) + 10 + additional_disk_gb

    command <<<
        set -euo pipefail

        # original VCF gets its sample name from the BAM, while the patch uses CDS-*
        echo "Reheadering VCFs to match sample names"

        echo "~{sample_id}" > sample_name.txt
        bcftools reheader --samples sample_name.txt "~{vcf}" \
            --output vcf_reheadered.vcf.gz \
            && rm "~{vcf}"
        bcftools reheader --samples sample_name.txt "~{vcf_patch}" \
            --output vcf_patch_reheadered.vcf.gz \
            && rm "~{vcf_patch}"

        echo "Indexing reheadered VCFs"
        bcftools index vcf_reheadered.vcf.gz
        bcftools index vcf_patch_reheadered.vcf.gz

        echo "Removing variants in patch regions from the base VCF"
        bcftools view vcf_reheadered.vcf.gz \
            --targets-file "^~{patch_region_bed}" \
            --output base_filtered.vcf.gz \
            && rm vcf_reheadered.vcf.gz

        bcftools index base_filtered.vcf.gz

        echo "Concatenating base (with patch regions removed) and patch VCFs"
        bcftools concat \
            --allow-overlaps \
            base_filtered.vcf.gz \
            vcf_patch_reheadered.vcf.gz \
            --output-type z \
            | bcftools sort --output "~{output_file_base_name}.vcf.gz"
    >>>

    output {
        File vcf_patched = "~{output_file_base_name}.vcf.gz"
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

task fix_with_bcftools {
    meta {
        description: "Apply miscellaneous fixes and filters to a somatic variant VCF"
        allowNestedInputs: true
    }

    parameter_meta {
        # inputs
        vcf: "input VCF to fix and filter"
        output_file_base_name: "base name for the output file (without extension)"
        fix_ploidy: "if true, correct multi-allelic genotypes where all reads support the alt"
        filter_vcf: "if true, apply bcftools exclude filter using exclude_string"
        normalize_indels: "if true, split multi-allelic sites and left-align indels"
        ref_fasta: "reference FASTA; required when normalize_indels is true"
        ref_fasta_index: "index of ref_fasta"
        exclude_string: "bcftools exclude expression used when filter_vcf is true"

        # outputs
        vcf_fixed: "bgzip-compressed VCF after all fixes and filters have been applied"
    }

    input {
        File vcf
        String output_file_base_name
        Boolean fix_ploidy
        Boolean filter_vcf
        Boolean normalize_indels
        File? ref_fasta
        File? ref_fasta_index
        String? exclude_string

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools"
        String docker_image_hash_or_tag = ":production"
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 0
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        ceil(size(vcf, "GiB") + 2 * 10 * size(vcf, "GiB")) + 10 + additional_disk_gb
    )

    command <<<
        set -euo pipefail

        TMP_VCF="~{basename(vcf, '.vcf.gz')}_tmp.vcf.gz"

        # fix for https://github.com/broadinstitute/gatk/issues/6857
        # temporarily convert to text VCF
        echo "Fixing AS_FilterStatus allele delimiter"
        bcftools view "~{vcf}" -o tmp.vcf && rm "~{vcf}"

        awk '{
            # find the filter status field
            match($0, /\tAS_FilterStatus=[^;]+;/);

            if (RSTART != 0) {
                # the line matches
                before = substr($0, 1, RSTART-1);
                match_str = substr($0, RSTART, RLENGTH);
                after = substr($0, RSTART+RLENGTH);

                # temp replace "|" with a non-printing char to swap "|" and "," chars
                gsub(/\|/, "\x1e", match_str);
                gsub(/,/, "|", match_str);
                gsub(/\x1e/, ",", match_str);

                # print modified line
                print before match_str after;
            } else {
                # no match
                print $0;
            }
        }' tmp.vcf > swapped.vcf && rm tmp.vcf

        bcftools view swapped.vcf -o "~{vcf}" && rm swapped.vcf

        if ~{fix_ploidy}; then
            echo "Fixing ploidy"
            # set the genotype annotation to be homozygous when we have no ref reads and at
            # least 3 alt reads, e.g., 0/1/2 --> 1/2, 0/1/2/3 --> 1/2/3, etc.
            bcftools +setGT "~{vcf}" \
                -- -t q -i'INFO/DP>8 & AF>0.9' -n c:'m|m' \
                > "${TMP_VCF}"
            bcftools +setGT "${TMP_VCF}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2"' -n c:'1/2' \
                > "~{vcf}"
            bcftools +setGT "~{vcf}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3"' -n c:'1/2/3' \
                > "${TMP_VCF}"
            bcftools +setGT "${TMP_VCF}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4"' -n c:'1/2/3/4' \
                > "~{vcf}"
            bcftools +setGT "~{vcf}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5"' -n c:'1/2/3/4/5' \
                > "${TMP_VCF}"
            bcftools +setGT "${TMP_VCF}" \
                -- -t q -i'AD[*:0]=0 & INFO/DP>8 & GT="0/1/2/3/4/5/6"' -n c:'1/2/3/4/5/6' \
                > "~{vcf}"
        fi

        if ~{filter_vcf}; then
            echo "Filtering VCF: ~{exclude_string}"
            bcftools view \
                "~{vcf}" \
                --exclude='~{exclude_string}' \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        if ~{normalize_indels}; then
            echo "Normalizing indels"
            bcftools norm "~{vcf}" \
                --multiallelics=- \
                --site-win=10000 \
                --fasta-ref="~{ref_fasta}" \
                --output="${TMP_VCF}"
            rm "~{vcf}" && mv "${TMP_VCF}" "~{vcf}"
        fi

        # ensure it's bgzipped
        bcftools view \
            "~{vcf}" \
            --output="~{output_file_base_name}.vcf.gz" \
            --no-version
    >>>

    output {
        File vcf_fixed = "~{output_file_base_name}.vcf.gz"
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
