version 1.0

import "subworkflows/align.wdl" as align
import "subworkflows/bucketize.wdl" as bucketize
import "subworkflows/cross.wdl" as cross
import "subworkflows/complex.wdl" as complex
import "subworkflows/blast.wdl" as blast
import "subworkflows/locate.wdl" as locate
import "subworkflows/metrics.wdl" as metrics
import "tasks/trimmomatic.wdl" as trimmomatic
import "tasks/quality.wdl" as quality
import "tasks/structs/compute.wdl"

workflow donor_recipient {

    input {
        File donor_ref_genome
        File recipient_ref_genome
        File reads_fastq_1
        File reads_fastq_2

        File blastdb
    }

    # Compute resources
    Compute server = read_json("../config/sizes.json")

    # Quality Control metrics
    call quality.fast_qc as quality_before_trim {
        input:
            fastq_1 = reads_fastq_1,
            fastq_2 = reads_fastq_2,
            resources = server.size["local_server"]
    }

    # Quality trim the reads
    call trimmomatic.trim as trim {
        input:
            fastq_1 = reads_fastq_1,
            fastq_2 = reads_fastq_2,
            adapter = "all_adapters.fa",
            seed_mismatches = 2,
            paired_clip_threshold = 30,
            unpaired_clip_threshold = 10,
            leading = 3,
            trailing = 3,
            sliding_window_quality = 20,
            sliding_window_length = 4,
            min_length = 50,
            is_phred33 = true,
            resources = server.size["local_server"]
    }

    # Quality control metrics after trimmomatic
    call quality.fast_qc as quality_after_trim {
        input:
            fastq_1 = trim.fastq_1_paired,
            fastq_2 = trim.fastq_2_paired,
            resources = server.size["local_server"]
    }

    # Compare quality control metrics before and after trimmomatic
    call quality.multi_qc as trimmomatic_metrics {
        input:
            quality_files = flatten(
                [
                quality_before_trim.files,
                quality_after_trim.files
                ]),
            report_name = "multiqc_before_and_after_trimmomatic_report.html",
            resources = server.size["local_server"]
    }

    # Donor Reference Genome
    ## Align
    call align.main as align_donor {
        input:
            ref_genome = donor_ref_genome,
            fastq_1 = trim.fastq_1_paired,
            fastq_2 = trim.fastq_2_paired,
            resources = server.size["local_server"],
            ignore_matches_shorted_than = 19,
            ignore_gaps_longer_than = 100,
            discard_if_repeated_in_ref_genome_more_than = 10000,
            matching_score = 1,
            mismatch_penalty = 4,
            gap_open_penalty = 6,
            gap_extension_penalty = 1,
            output_all_found_alignments = true,
            base_filename = "reads-to-donor"
    }

    ## Bucketize
    call bucketize.main as bucketize_donor {
        input:
            bam = align_donor.bam,
            resources = server.size["local_server"],
            base_filename = "reads-to-donor"
    }

    # Recipient Reference Genome
    ## Align 
    call align.main as align_recipient {
        input:
            ref_genome = recipient_ref_genome,
            fastq_1 = trim.fastq_1_paired,
            fastq_2 = trim.fastq_2_paired,
            resources = server.size["local_server"],
            ignore_matches_shorted_than = 19,
            ignore_gaps_longer_than = 100,
            discard_if_repeated_in_ref_genome_more_than = 10000,
            matching_score = 1,
            mismatch_penalty = 4,
            gap_open_penalty = 6,
            gap_extension_penalty = 1,
            output_all_found_alignments = true,
            base_filename = "reads-to-recipient"
    }

    ## Bucketize
    call bucketize.main as bucketize_recipient {
        input:
            bam = align_recipient.bam,
            resources = server.size["local_server"],
            base_filename = "reads-to-recipient"
    }

    # Cross BAM files buckets
    call cross.main as crossing {
        input:
            donor_MM = bucketize_donor.MM,
            donor_MU = bucketize_donor.MU,
            donor_UM = bucketize_donor.UM,
            donor_UU = bucketize_donor.UU,
            recipient_MM = bucketize_recipient.MM,
            recipient_MU = bucketize_recipient.MU,
            recipient_UM = bucketize_recipient.UM,
            recipient_UU = bucketize_recipient.UU,
            resources = server.size["local_server"]
    }

    # Filter out low complexity sequences
    call complex.main as complex_only {
        input:
            bam_files = crossing.bams,
            filter_shorter_than = 5,
            filter_longer_than = 10000,
            filter_if_gc_content_lower_than = 10,
            filter_if_gc_content_higher_than = 90,
            filter_if_avg_quality_below = 20,
            low_complexity_method = 'dust',
            low_complexity_threshold = '7',
            resources = server.size["local_server"]
    }

    # Locate sequences in reference genome
    call locate.main as locate_donor {
        input:
            bams = complex_only.bams,
            ref_genome = donor_ref_genome,
            resources = server.size["local_server"]
    }

    call locate.main as locate_recipient {
        input:
            bams = complex_only.bams,
            ref_genome = recipient_ref_genome,
            resources = server.size["local_server"]
    }

    # Blastn search over all crossings of interest (see tasks/blast.wdl)
    call blast.main as blaster {
        input:
            fastas = complex_only.fastas,
            blastdb = blastdb,
            resources = server.size["local_server"],
            evalue = 1
    }

    # Calculate metrics
    call metrics.main as metrics_donor_bucketized {
        input:
            bams = bucketize_donor.bams,
            resources = server.size["local_server"]
    }

    call metrics.main as metrics_recipient_bucketized {
        input:
            bams = bucketize_recipient.bams,
            resources = server.size["local_server"]
    }

    call metrics.main as metrics_cross {
        input:
            bams = crossing.bams,
            resources = server.size["local_server"]
    }
    
    call metrics.main as metrics_complex {
        input:
            bams = complex_only.bams,
            resources = server.size["local_server"]
    }

    # Compare quality control metrics for crossed bams
    call quality.multi_qc as crossed_multiqc_metrics {
        input:
            quality_files = flatten(
                [
                quality_before_trim.files,
                quality_after_trim.files
                ]),
            report_name = "multiqc_before_and_after_trimmomatic_report.html",
            resources = server.size["local_server"]
    }


    output {
        File out_pre_fastq_1_zip = quality_before_trim.fastq_1_zip
        File out_pre_fastq_2_zip = quality_before_trim.fastq_2_zip
        File out_pre_fastq_1_html = quality_before_trim.fastq_1_html
        File out_pre_fastq_2_html = quality_before_trim.fastq_2_html

        File out_post_fastq_1_zip = quality_after_trim.fastq_1_zip
        File out_post_fastq_2_zip = quality_after_trim.fastq_2_zip
        File out_post_fastq_1_html = quality_after_trim.fastq_1_html
        File out_post_fastq_2_html = quality_after_trim.fastq_2_html

        File out_fastq_1_paired = trim.fastq_1_paired
        File out_fastq_1_unpaired = trim.fastq_1_unpaired
        File out_fastq_2_paired = trim.fastq_2_paired
        File out_fastq_2_unpaired = trim.fastq_2_unpaired

        File out_multiqc_before_and_after_trim_report = trimmomatic_metrics.out

        Array[File] out_donor_bams = bucketize_donor.bams
        Array[File] out_recipient_bams = bucketize_recipient.bams
        Array[File] out_cross_bams = crossing.bams
        Array[File] out_complex_bams = complex_only.bams
        Array[File] out_fastas = complex_only.fastas
        
        Array[File] out_donor_mpileups = locate_donor.mpileups
        Array[File] out_recipient_mpileups = locate_recipient.mpileups

        File out_donor_MMd_MUr = blaster.donor_MMd_MUr
        File out_donor_MUd_UMr = blaster.donor_MUd_UMr
        File out_donor_UMd_MUr = blaster.donor_UMd_MUr
        File out_recipient_MMd_MUr = blaster.recipient_MMd_MUr
        File out_recipient_MUd_UMr = blaster.recipient_MUd_UMr
        File out_recipient_UMd_MUr = blaster.recipient_UMd_MUr

        Array[File] out_donor_stats = metrics_donor_bucketized.stats
        Array[File] out_donor_flagstats = metrics_donor_bucketized.flagstats
        Array[File] out_donor_count = metrics_donor_bucketized.count

        Array[File] out_recipient_stats = metrics_recipient_bucketized.stats
        Array[File] out_recipient_flagstats = metrics_recipient_bucketized.flagstats
        Array[File] out_recipient_count = metrics_recipient_bucketized.count

        Array[File] out_cross_stats = metrics_cross.stats
        Array[File] out_cross_flagstats = metrics_cross.flagstats
        Array[File] out_cross_count = metrics_cross.count

        Array[File] out_complex_stats = metrics_complex.stats
        Array[File] out_complex_flagstats = metrics_complex.flagstats
        Array[File] out_complex_count = metrics_complex.count

        File out_crossed_multiqc_report = crossed_multiqc_metrics.out

    }

}



