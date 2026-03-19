version 1.0

workflow align_sr_rna {
    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        File star_index
        String? align_insertion_flush
        Int? align_intron_max
        Int? align_mates_gap_max
        Int? align_sjdb_overhang_min
        String? align_sj_stitch_mismatch_nmax
        Int? align_spliced_mate_map_lmin
        Int? align_spliced_mate_map_lmin_over_lmate
        Int? chim_junction_overhang_min
        Int? chim_multimap_nmax
        Int? chim_multimap_score_range
        Int? chim_nonchim_score_drop_min
        Int? chim_out_junction_format
        String? chim_out_type
        Int? chim_score_junction_non_gtag
        Int? chim_segment_min
        String? genome_dir
        String? genome_load
        Int? limit_sjdb_insert_nsj
        String? out_filter_intron_motifs
        Int? out_filter_match_nmin
        Float? out_filter_match_nmin_over_lread
        Int? out_filter_mismatch_nmax
        Float? out_filter_mismatch_nover_lmax
        Int? out_filter_multimap_nmax
        Float? out_filter_score_min_over_lread
        String? out_filter_type
        String? out_reads_unmapped
        String? out_sam_attr_rg_line
        String? out_sam_attributes
        String? out_sam_strand_field
        String? out_sam_type
        String? out_sam_unmapped
        Float? pe_overlap_mmp
        Int? pe_overlap_nbases_min
        String? quant_mode
        String? read_files_command
        String? twopass_mode
        Int? win_anchor_multimap_nmax
        String? extra_star_args
    }

    call align_with_star {
        input:
            sample_id = sample_id,
            input_file_type = input_file_type,
            fastqs = fastqs,
            cram_bam = cram_bam,
            crai_bai = crai_bai,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            star_index = star_index,
            align_insertion_flush = align_insertion_flush,
            align_intron_max = align_intron_max,
            align_mates_gap_max = align_mates_gap_max,
            align_sjdb_overhang_min = align_sjdb_overhang_min,
            align_sj_stitch_mismatch_nmax = align_sj_stitch_mismatch_nmax,
            align_spliced_mate_map_lmin = align_spliced_mate_map_lmin,
            align_spliced_mate_map_lmin_over_lmate = align_spliced_mate_map_lmin_over_lmate,
            chim_junction_overhang_min = chim_junction_overhang_min,
            chim_multimap_nmax = chim_multimap_nmax,
            chim_multimap_score_range = chim_multimap_score_range,
            chim_nonchim_score_drop_min = chim_nonchim_score_drop_min,
            chim_out_junction_format = chim_out_junction_format,
            chim_out_type = chim_out_type,
            chim_score_junction_non_gtag = chim_score_junction_non_gtag,
            chim_segment_min = chim_segment_min,
            genome_dir = genome_dir,
            genome_load = genome_load,
            limit_sjdb_insert_nsj = limit_sjdb_insert_nsj,
            out_filter_intron_motifs = out_filter_intron_motifs,
            out_filter_match_nmin = out_filter_match_nmin,
            out_filter_match_nmin_over_lread = out_filter_match_nmin_over_lread,
            out_filter_mismatch_nmax = out_filter_mismatch_nmax,
            out_filter_mismatch_nover_lmax = out_filter_mismatch_nover_lmax,
            out_filter_multimap_nmax = out_filter_multimap_nmax,
            out_filter_score_min_over_lread = out_filter_score_min_over_lread,
            out_filter_type = out_filter_type,
            out_reads_unmapped = out_reads_unmapped,
            out_sam_attr_rg_line = out_sam_attr_rg_line,
            out_sam_attributes = out_sam_attributes,
            out_sam_strand_field = out_sam_strand_field,
            out_sam_type = out_sam_type,
            out_sam_unmapped = out_sam_unmapped,
            pe_overlap_mmp = pe_overlap_mmp,
            pe_overlap_nbases_min = pe_overlap_nbases_min,
            quant_mode = quant_mode,
            read_files_command = read_files_command,
            twopass_mode = twopass_mode,
            win_anchor_multimap_nmax = win_anchor_multimap_nmax,
            extra_star_args = extra_star_args
    }

    output {
        File analysis_ready_bam = align_with_star.analysis_ready_bam
        File junctions = align_with_star.junctions
        File transcriptome_bam = align_with_star.transcriptome_bam
        File reads_per_gene = align_with_star.reads_per_gene
    }
}

task align_with_star {
    input {
        String sample_id
        String input_file_type
        Array[File]? fastqs
        File? cram_bam
        File? crai_bai
        File? ref_fasta
        File? ref_fasta_index
        File star_index
        String? align_insertion_flush
        Int? align_intron_max
        Int? align_mates_gap_max
        Int? align_sjdb_overhang_min
        String? align_sj_stitch_mismatch_nmax
        Int? align_spliced_mate_map_lmin
        Int? align_spliced_mate_map_lmin_over_lmate
        Int? chim_junction_overhang_min
        Int? chim_multimap_nmax
        Int? chim_multimap_score_range
        Int? chim_nonchim_score_drop_min
        Int? chim_out_junction_format
        String? chim_out_type
        Int? chim_score_junction_non_gtag
        Int? chim_segment_min
        String? genome_dir
        String? genome_load
        Int? limit_sjdb_insert_nsj
        String? out_filter_intron_motifs
        Int? out_filter_match_nmin
        Float? out_filter_match_nmin_over_lread
        Int? out_filter_mismatch_nmax
        Float? out_filter_mismatch_nover_lmax
        Int? out_filter_multimap_nmax
        Float? out_filter_score_min_over_lread
        String? out_filter_type
        String? out_reads_unmapped
        String? out_sam_attr_rg_line
        String? out_sam_attributes
        String? out_sam_strand_field
        String? out_sam_type
        String? out_sam_unmapped
        Float? pe_overlap_mmp
        Int? pe_overlap_nbases_min
        String? quant_mode
        String? read_files_command
        String? twopass_mode
        Int? win_anchor_multimap_nmax
        String? extra_star_args

        String docker_image = "us-central1-docker.pkg.dev/depmap-omics/terra-images/star_arriba"
        String docker_image_hash_or_tag = ":production"
        Int cpu = 16
        Int mem_gb = 48
        Int preemptible = 1
        Int max_retries = 1
        Int additional_disk_gb = 0
    }

    Int disk_space = (
        if input_file_type == "FASTQ" then
            ceil(size(select_first([fastqs, []]), "GiB"))
        else (
            if input_file_type == "CRAM" then
                ceil(size(select_first([cram_bam, "/dev/null"]), "GiB") * 8)
            else # BAM
                ceil(size(select_first([cram_bam, "/dev/null"]), "GiB") * 2)
        )
    ) + ceil(size(star_index, "GiB") * 5) + 20 + additional_disk_gb

    Int n_threads = cpu - 1

    command <<<
        set -euo pipefail

        if [[ "~{input_file_type}" == "CRAM" ]]; then
            echo "Converting CRAM to FASTQ"

            samtools sort -n \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                "~{cram_bam}" \
            | samtools fastq \
                -@ ~{n_threads} \
                --reference "~{ref_fasta}" \
                -1 "~{sample_id}.1.fq.gz" \
                -2 "~{sample_id}.2.fq.gz" -

            READ_FILES_IN="~{sample_id}.1.fq.gz ~{sample_id}.2.fq.gz"
        elif [[ "~{input_file_type}" == "BAM" ]]; then
            echo "Converting BAM to FASTQ"

            samtools sort -n \
                -@ ~{n_threads} \
                "~{cram_bam}" \
            | samtools fastq \
                -@ ~{n_threads} \
                -1 "~{sample_id}.1.fq.gz" \
                -2 "~{sample_id}.2.fq.gz" -

            READ_FILES_IN="~{sample_id}.1.fq.gz ~{sample_id}.2.fq.gz"
        elif [[ "~{input_file_type}" == "FASTQ" ]]; then
            READ_FILES_IN="~{sep=' ' fastqs}"
        fi

        echo "Extracting STAR index"
        mkdir -p "star_index"
        tar -C "star_index" -xzf "~{star_index}" --strip-components=1

        echo "Running STAR"
        STAR \
            --readFilesIn $READ_FILES_IN \
            --outFileNamePrefix "~{sample_id}." \
            --runThreadN ~{n_threads} \
            ~{"--alignInsertionFlush " + align_insertion_flush} \
            ~{"--alignIntronMax " + align_intron_max} \
            ~{"--alignMatesGapMax " + align_mates_gap_max} \
            ~{"--alignSJDBoverhangMin " + align_sjdb_overhang_min} \
            ~{"--alignSJstitchMismatchNmax " + align_sj_stitch_mismatch_nmax} \
            ~{"--alignSplicedMateMapLmin " + align_spliced_mate_map_lmin} \
            ~{"--alignSplicedMateMapLminOverLmate " + align_spliced_mate_map_lmin_over_lmate} \
            ~{"--chimJunctionOverhangMin " + chim_junction_overhang_min} \
            ~{"--chimMultimapNmax " + chim_multimap_nmax} \
            ~{"--chimMultimapScoreRange " + chim_multimap_score_range} \
            ~{"--chimNonchimScoreDropMin " + chim_nonchim_score_drop_min} \
            ~{"--chimOutJunctionFormat " + chim_out_junction_format} \
            ~{"--chimOutType " + chim_out_type} \
            ~{"--chimScoreJunctionNonGTAG " + chim_out_type} \
            ~{"--chimSegmentMin " + chim_segment_min} \
            ~{"--genomeDir " + genome_dir} \
            ~{"--genomeLoad " + genome_load} \
            ~{"--limitSjdbInsertNsj " + limit_sjdb_insert_nsj} \
            ~{"--outFilterIntronMotifs " + out_filter_intron_motifs} \
            ~{"--outFilterMatchNmin " + out_filter_match_nmin} \
            ~{"--outFilterMatchNminOverLread " + out_filter_match_nmin_over_lread} \
            ~{"--outFilterMismatchNmax " + out_filter_mismatch_nmax} \
            ~{"--outFilterMismatchNoverLmax " + out_filter_mismatch_nover_lmax} \
            ~{"--outFilterMultimapNmax " + out_filter_multimap_nmax} \
            ~{"--outFilterScoreMinOverLread " + out_filter_score_min_over_lread} \
            ~{"--outFilterType " + out_filter_type} \
            ~{"--outReadsUnmapped " + out_reads_unmapped} \
            ~{"--outSAMattrRGline " + out_sam_attr_rg_line} \
            ~{"--outSAMattributes " + out_sam_attributes} \
            ~{"--outSAMstrandField " + out_sam_strand_field} \
            ~{"--outSAMtype " + out_sam_type} \
            ~{"--outSAMunmapped " + out_sam_unmapped} \
            ~{"--peOverlapMMp " + pe_overlap_mmp} \
            ~{"--peOverlapNbasesMin " + pe_overlap_nbases_min} \
            ~{"--quantMode " + quant_mode} \
            ~{"--readFilesCommand " + read_files_command} \
            ~{"--twopassMode " + twopass_mode} \
            ~{"--winAnchorMultimapNmax " + win_anchor_multimap_nmax} \
            ~{extra_star_args}

        gzip -c "~{sample_id}.SJ.out.tab" > "~{sample_id}.junctions.tsv.gz"
        gzip -c "~{sample_id}.ReadsPerGene.out.tab" > "~{sample_id}.reads_per_gene.tsv.gz"
    >>>

    output {
        File analysis_ready_bam = "~{sample_id}.Aligned.out.bam"
        File junctions = "~{sample_id}.junctions.tsv.gz"
        File reads_per_gene = "~{sample_id}.reads_per_gene.tsv.gz"
        File transcriptome_bam = "~{sample_id}.Aligned.toTranscriptome.out.bam"
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
