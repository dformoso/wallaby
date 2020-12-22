version 1.0

import "subworkflows/align.wdl" as align
import "subworkflows/bucketize.wdl" as bucketize
import "subworkflows/cross.wdl" as cross
import "subworkflows/complex.wdl" as complex
#import "subworkflows/blast.wdl" as blast
#import "subworkflows/locate.wdl" as locate
import "subworkflows/metrics.wdl" as metrics
import "tasks/trimmomatic.wdl" as trimmomatic
import "tasks/quality.wdl" as quality
import "tasks/samtools.wdl" as samtools
import "tasks/structs/compute.wdl"

workflow main {

    input {
        String donor_name
        File donor_ref_genome
        String recipient_name
        File recipient_ref_genome
        String srr_name
        File srr_fastq_1
        File srr_fastq_2
        File blastdb
    }

    # Compute resources
    Compute server = read_json("../config/sizes.json")

    # Quality Control metrics
    call quality.fast_qc as srr_fastqc_before_trim {
        input:
            fastq_1 = srr_fastq_1,
            fastq_2 = srr_fastq_2,
            resources = server.size["2cpu_8mem_100disk"]
    }

    # Quality trim the reads
    call trimmomatic.trim as srr_trim_adapters {
        input:
            fastq_1 = srr_fastq_1,
            fastq_2 = srr_fastq_2,
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
            resources = server.size["8cpu_32mem_512disk"]
    }

    # Quality control metrics after trimmomatic
    call quality.fast_qc as srr_fastqc_after_trim {
        input:
            fastq_1 = srr_trim_adapters.fastq_1_paired,
            fastq_2 = srr_trim_adapters.fastq_2_paired,
            resources = server.size["2cpu_8mem_100disk"]
    }

    # Compare quality control metrics before and after trimmomatic
    call quality.multi_qc as srr_multiqc_after_trim {
        input:
            quality_files = flatten(
                [
                srr_fastqc_before_trim.files,
                srr_fastqc_after_trim.files
                ]),
            report_name = "${srr_name}_multiqc_trim_report.html",
            include = "../inputs/*",            
            resources = server.size["8cpu_32mem_512disk"]
    }

    # Donor Reference Genome
    ## Align
    call align.main as donor_align {
        input:
            ref_genome = donor_ref_genome,
            fastq_1 = srr_trim_adapters.fastq_1_paired,
            fastq_2 = srr_trim_adapters.fastq_2_paired,
            ignore_matches_shorted_than = 19,
            ignore_gaps_longer_than = 100,
            discard_if_repeated_in_ref_genome_more_than = 10000,
            matching_score = 1,
            mismatch_penalty = 4,
            gap_open_penalty = 6,
            gap_extension_penalty = 1,
            output_all_found_alignments = true,
            base_filename = "${srr_name}-to-${donor_name}",
            resources = server.size["32cpu_64mem_512disk"]
    }

    ## Bucketize
    call bucketize.main as donor_bucketize {
        input:
            bam = donor_align.bam,
            base_filename = "${srr_name}-to-${donor_name}",
            resources = server.size["32cpu_64mem_512disk"]
    }

    # Recipient Reference Genome
    ## Align 
    call align.main as recipient_align {
        input:
            ref_genome = recipient_ref_genome,
            fastq_1 = srr_trim_adapters.fastq_1_paired,
            fastq_2 = srr_trim_adapters.fastq_2_paired,
            ignore_matches_shorted_than = 19,
            ignore_gaps_longer_than = 100,
            discard_if_repeated_in_ref_genome_more_than = 10000,
            matching_score = 1,
            mismatch_penalty = 4,
            gap_open_penalty = 6,
            gap_extension_penalty = 1,
            output_all_found_alignments = true,
            base_filename = "${srr_name}-to-${recipient_name}",
            resources = server.size["32cpu_64mem_512disk"]
    }

    ## Bucketize
    call bucketize.main as recipient_bucketize {
        input:
            bam = recipient_align.bam,
            base_filename = "${srr_name}-to-${recipient_name}",
            resources = server.size["32cpu_64mem_512disk"]
    }

    # Cross BAM files buckets
    call cross.main as crossing {
        input:
            donor_MM = donor_bucketize.MM,
            donor_MU = donor_bucketize.MU,
            donor_UM = donor_bucketize.UM,
            donor_UU = donor_bucketize.UU,
            recipient_MM = recipient_bucketize.MM,
            recipient_MU = recipient_bucketize.MU,
            recipient_UM = recipient_bucketize.UM,
            recipient_UU = recipient_bucketize.UU,
            donor_name = "${donor_name}",
            recipient_name = "${recipient_name}",
            srr_name = "${srr_name}",
            resources = server.size["32cpu_64mem_512disk"]
    }

    # Filter out low complexity sequences
    call complex.main as crossed_filtered {
        input:
            bam_files = crossing.bams,
            filter_shorter_than = 5,
            filter_longer_than = 10000,
            filter_if_gc_content_lower_than = 10,
            filter_if_gc_content_higher_than = 90,
            filter_if_avg_quality_below = 20,
            low_complexity_method = 'dust',
            low_complexity_threshold = '7',
            resources = server.size["8cpu_32mem_512disk"]
    }

    # Create indexes (BAI files) for all crossed_filtered BAM files
    scatter (bam in crossed_filtered.bams) {
        call samtools.index as indexing_bams {
            input:
                file = bam,
                resources = server.size["4cpu_32mem_100disk"]
        }
    }

    # Create BED files for all crossed_filtered BAM files
    scatter (bam in crossed_filtered.bams) {
        call samtools.bam_to_bed as bams_to_beds {
            input:
                file = bam,
                resources = server.size["4cpu_32mem_100disk"]
        }
    }

    # Locate sequences in reference genome
#    call locate.main as donor_locate {
#        input:
#            bams = crossed_filtered.bams,
#            ref_genome = donor_ref_genome,
#            resources = server.size["4cpu_32mem_100disk"]
#    }
#
#    call locate.main as recipient_locate {
#        input:
#            bams = crossed_filtered.bams,
#            ref_genome = recipient_ref_genome,
#            resources = server.size["4cpu_32mem_100disk"]
#    }

    # Blastn search over all crossings of interest (see tasks/blast.wdl)
#    call blast.main as blaster {
#        input:
#            fastas = crossed_filtered.fastas,
#            blastdb = blastdb,
#            resources = server.size["32cpu_64mem_200disk"],
#            evalue = 1
#    }

    # Calculate metrics
    call metrics.main as donor_bucketized_metrics {
        input:
            bams = donor_bucketize.bams,
            resources = server.size["4cpu_32mem_100disk"]
    }

    call metrics.main as recipient_bucketized_metrics {
        input:
            bams = recipient_bucketize.bams,
            resources = server.size["4cpu_32mem_100disk"]
    }
    
    call metrics.main as crossed_filtered_metrics {
        input:
            bams = crossed_filtered.bams,
            resources = server.size["4cpu_32mem_100disk"]
    }

    # Compare quality control for all donor files
    call quality.multi_qc as donor_crossed_filtered_multiqc {
        input:
            quality_files = flatten(
                [
                crossed_filtered.bams,
                crossed_filtered_metrics.stats,
                crossed_filtered_metrics.flagstats
                ]),
            report_name = "${srr_name}-to-${donor_name}_multiqc_metrics.html",
            enable_fullnames = false,
            include = "../inputs/*/*${donor_name}*",
            resources = server.size["4cpu_32mem_100disk"]
    }

    # Compare quality control for all recipient files
    call quality.multi_qc as recipient_crossed_filtered_multiqc {
        input:
            quality_files = flatten(
                [
                crossed_filtered.bams,
                crossed_filtered_metrics.stats,
                crossed_filtered_metrics.flagstats
                ]),
            report_name = "${srr_name}-to-${recipient_name}_multiqc_metrics.html",
            enable_fullnames = false,
            include = "../inputs/*/*${recipient_name}*",
            resources = server.size["4cpu_32mem_100disk"]
    }

    output {
        File out_pre_fastq_1_zip = srr_fastqc_before_trim.fastq_1_zip
        File out_pre_fastq_2_zip = srr_fastqc_before_trim.fastq_2_zip
        File out_pre_fastq_1_html = srr_fastqc_before_trim.fastq_1_html
        File out_pre_fastq_2_html = srr_fastqc_before_trim.fastq_2_html

        File out_post_fastq_1_zip = srr_fastqc_after_trim.fastq_1_zip
        File out_post_fastq_2_zip = srr_fastqc_after_trim.fastq_2_zip
        File out_post_fastq_1_html = srr_fastqc_after_trim.fastq_1_html
        File out_post_fastq_2_html = srr_fastqc_after_trim.fastq_2_html

        File out_fastq_1_paired = srr_trim_adapters.fastq_1_paired
        File out_fastq_1_unpaired = srr_trim_adapters.fastq_1_unpaired
        File out_fastq_2_paired = srr_trim_adapters.fastq_2_paired
        File out_fastq_2_unpaired = srr_trim_adapters.fastq_2_unpaired

        File out_multiqc_before_and_after_trim_report = srr_multiqc_after_trim.out

        File out_donor_bam = donor_align.bam
        File out_donor_bai = donor_align.bai
        File out_recipient_bam = recipient_align.bam
        File out_recipient_bai = recipient_align.bai

        Array[File] out_bucket_donor_bams = donor_bucketize.bams
        Array[File] out_bucket_donor_bais = donor_bucketize.bais
        Array[File] out_bucket_recipient_bams = recipient_bucketize.bams
        Array[File] out_bucket_recipient_bais = recipient_bucketize.bais

        Array[File] out_crossed_filtered_bams = crossed_filtered.bams
#        Array[File] out_crossed_filtered_fastas = crossed_filtered.fastas
        Array[File] out_crossed_filtered_bais = indexing_bams.out
        Array[File] out_crossed_filtered_beds = bams_to_beds.out

#        Array[File] out_donor_mpileups = donor_locate.mpileups
#        Array[File] out_recipient_mpileups = recipient_locate.mpileups

#        File out_donor_MMd_MUr = blaster.donor_MMd_MUr
#        File out_donor_MUd_UMr = blaster.donor_MUd_UMr
#        File out_donor_UMd_MUr = blaster.donor_UMd_MUr
#        File out_recipient_MMd_MUr = blaster.recipient_MMd_MUr
#        File out_recipient_MUd_UMr = blaster.recipient_MUd_UMr
#        File out_recipient_UMd_MUr = blaster.recipient_UMd_MUr

        Array[File] out_donor_stats = donor_bucketized_metrics.stats
        Array[File] out_donor_flagstats = donor_bucketized_metrics.flagstats

        Array[File] out_recipient_stats = recipient_bucketized_metrics.stats
        Array[File] out_recipient_flagstats = recipient_bucketized_metrics.flagstats

        Array[File] out_crossed_filtered_stats = crossed_filtered_metrics.stats
        Array[File] out_crossed_filtered_flagstats = crossed_filtered_metrics.flagstats

        File out_multiqc_donor_crossed_filtered_multiqc_report = donor_crossed_filtered_multiqc.out
        File out_multiqc_recipient_crossed_filtered_multiqc_report = recipient_crossed_filtered_multiqc.out

    }
}



