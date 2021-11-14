version 1.0

import "subworkflows/bucketize.wdl" as bucketize
import "subworkflows/cross.wdl" as cross
import "subworkflows/metrics.wdl" as metrics
import "tasks/align.wdl" as align
import "tasks/quality.wdl" as quality
import "tasks/samtools.wdl" as samtools
import "tasks/filter.wdl" as filter
import "tasks/overlaps.wdl" as overlaps
import "tasks/structs/compute.wdl"

workflow main {

    input {
        String donor_name
        Array[File] donor_index
        File donor_ref_genome

        String recipient_name
        Array[File] recipient_index
        File recipient_ref_genome

        String srr_name
        File fastq_1
        File fastq_2
    }

    # Compute resources
    Compute server = read_json("../inputs/sizes.json")

    # Donor Reference Genome
    
    ## Align
    call align.align_convert_index as donor_align {
        input:
            index_object = donor_index,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            base_filename = "${srr_name}-to-${donor_name}",
            resources = server.size["local_instance"]
    }

    ## Bucketize
    call bucketize.main as donor_bucketize {
        input:
            bam = donor_align.bam,
            base_filename = "${srr_name}-to-${donor_name}",
            resources = server.size["local_instance"]
    }

    # Recipient Reference Genome

    ## Align 
    call align.align_convert_index as recipient_align {
        input:
            index_object = recipient_index,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            base_filename = "${srr_name}-to-${recipient_name}",
            resources = server.size["local_instance"]
    }

    ## Bucketize
    call bucketize.main as recipient_bucketize {
        input:
            bam = recipient_align.bam,
            base_filename = "${srr_name}-to-${recipient_name}",
            resources = server.size["local_instance"]
    }

    # Cross BAM files buckets
    call cross.main as crossing {
        input:
            donor_MM = donor_bucketize.MM,
            donor_MU = donor_bucketize.MU,
            donor_UM = donor_bucketize.UM,
            recipient_MM = recipient_bucketize.MM,
            recipient_MU = recipient_bucketize.MU,
            recipient_UM = recipient_bucketize.UM,
            donor_name = "${donor_name}",
            recipient_name = "${recipient_name}",
            srr_name = "${srr_name}",
            resources = server.size["local_instance"]
    }

    # Filter out low complexity sequences
    scatter (bam in crossing.bams) {
        call filter.bam_reads as filtered {
            input:
                bam = bam,
                filter = "true",
                validation_stringency = "SILENT",
                filter_shorter_than = 5,
                filter_longer_than = 10000,
                filter_if_avg_quality_below = 20,
                filter_if_gc_content_lower_than = 10,
                filter_if_gc_content_higher_than = 90,
                low_complexity_method = 'dust',
                low_complexity_threshold = '7',
                filter_type = "includeReadList",
                resources = server.size["local_instance"]
        }
    }

    # Create indexes (BAI files) for all crossed_filtered BAM files
    scatter (bam in filtered.bams) {
        call samtools.index as filtered_bam_to_bai {
            input:
                file = bam,
                resources = server.size["local_instance"]
        }
    }

    # Create BED files for all crossed_filtered BAM files
    scatter (bam in filtered.bams) {
        call samtools.bam_to_bed as filtered_bam_to_bed {
            input:
                file = bam,
                resources = server.size["local_instance"]
        }
    }

    # Calculate metrics for crossed and filtered files
    call metrics.main as crossing_metrics {
        input:
            bams = crossing.bams,
            resources = server.size["local_instance"]
    }
    
    call metrics.main as filtered_metrics {
        input:
            bams = filtered.bams,
            resources = server.size["local_instance"]
    }

    # Compare quality control between crossed and filtered donor files
    call quality.multi_qc as donor_multiqc {
        input:
            quality_files = flatten(
                [
                crossing.bams,
                crossing_metrics.stats,
                crossing_metrics.flagstats,
                filtered.bams,
                filtered_metrics.stats,
                filtered_metrics.flagstats
                ]),
            report_name = "${srr_name}-to-${donor_name}_multiqc_report.html",
            enable_fullnames = false,
            include = "../inputs/*/*${donor_name}*",
            resources = server.size["local_instance"]
    }

    # Compare quality control between crossed and filtered recipient files
    call quality.multi_qc as recipient_multiqc {
        input:
            quality_files = flatten(
                [
                crossing.bams,
                crossing_metrics.stats,
                crossing_metrics.flagstats,
                filtered.bams,
                filtered_metrics.stats,
                filtered_metrics.flagstats
                ]),
            report_name = "${srr_name}-to-${recipient_name}_multiqc_report.html",
            enable_fullnames = false,
            include = "../inputs/*/*${recipient_name}*",
            resources = server.size["local_instance"]
    }

    # Create putative insertion table and bam files for loci of interest
    call overlaps.putative_insertions as overlap_loci {
        input:
            bams = filtered_bams,
            bais = filtered_bam_to_bai.out,
            beds = filtered_bam_to_bed.out,
            srr_name = srr_name,
            donor_name = donor_name,
            donor_ref_genome = donor_ref_genome,
            recipient_name = recipient_name,
            recipient_ref_genome = recipient_ref_genome,
            min_num_crossings = "1",
            min_num_reads = "5",
            resources = server.size["local_instance"]
    }

    # Create indexes (BAI files) for all overlap loci BAM files
    scatter (bam in overlap_loci.out) {
        call samtools.index as overlap_bam_to_bai {
            input:
                file = bam,
                resources = server.size["local_instance"]
        }
    }

    output {
        Array[File] filtered_bams = filtered.bams
        Array[File] filtered_bais = filtered_bam_to_bai.out
        Array[File] filtered_beds = filtered_bam_to_bed.out
        Array[File] overlap_loci_bams = overlap_loci.out
        Array[File] overlap_loci_bais = overlap_bam_to_bai.out
        File overlap_loci_csv = overlap_loci.csv

        File? multiqc_donor_html = donor_multiqc.html
        File? multiqc_donor_zip = donor_multiqc.zip
        File? multiqc_recipient_html = recipient_multiqc.html
        File? multiqc_recipient_zip = recipient_multiqc.zip
    }
}