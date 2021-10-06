version 1.0

#import "subworkflows/align.wdl" as align
import "subworkflows/bucketize.wdl" as bucketize
import "subworkflows/cross.wdl" as cross
#import "subworkflows/filter.wdl" as filter
import "subworkflows/metrics.wdl" as metrics
import "tasks/align.wdl" as align
import "tasks/quality.wdl" as quality
import "tasks/samtools.wdl" as samtools
import "tasks/structs/compute.wdl"

workflow main {

    input {
        String donor_name
        Array[File] donor_index

        String recipient_name
        Array[File] recipient_index

        String srr_name
        File fastq_1
        File fastq_2
    }

    # Compute resources
    Compute server = read_json("../config/sizes.json")

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
            #donor_UU = donor_bucketize.UU,
            recipient_MM = recipient_bucketize.MM,
            recipient_MU = recipient_bucketize.MU,
            recipient_UM = recipient_bucketize.UM,
            #recipient_UU = recipient_bucketize.UU,
            donor_name = "${donor_name}",
            recipient_name = "${recipient_name}",
            srr_name = "${srr_name}",
            resources = server.size["local_instance"]
    }

    # Filter out low complexity sequences
#    call filter.main as filtered {
#        input:
#            bam_files = crossing.bams,
#            filter_shorter_than = 5,
#            filter_longer_than = 10000,
#            filter_if_avg_quality_below = 20,
#            filter_if_gc_content_lower_than = 10,
#            filter_if_gc_content_higher_than = 90,
#            low_complexity_method = 'dust',
#            low_complexity_threshold = '7',
#            resources = server.size["local_instance"]
#    }

    # Create indexes (BAI files) for all crossed_filtered BAM files
    scatter (bam in crossing.bams) {
        call samtools.index as indexing_bams {
            input:
                file = bam,
                resources = server.size["local_instance"]
        }
    }

    # Create BED files for all crossed_filtered BAM files
    scatter (bam in crossing.bams) {
        call samtools.bam_to_bed as bams_to_beds {
            input:
                file = bam,
                resources = server.size["local_instance"]
        }
    }

    # Calculate metrics
    call metrics.main as donor_bucketized_metrics {
        input:
            bams = donor_bucketize.bams,
            resources = server.size["local_instance"]
    }

    call metrics.main as recipient_bucketized_metrics {
        input:
            bams = recipient_bucketize.bams,
            resources = server.size["local_instance"]
    }
    
    call metrics.main as crossing_metrics {
        input:
            bams = crossing.bams,
            resources = server.size["local_instance"]
    }

    # Compare quality control for all donor files
    call quality.multi_qc as donor_crossing_multiqc {
        input:
            quality_files = flatten(
                [
                crossing.bams,
                crossing_metrics.stats,
                crossing_metrics.flagstats
                ]),
            report_name = "${srr_name}-to-${donor_name}_multiqc_metrics.html",
            enable_fullnames = false,
            include = "../inputs/*/*${donor_name}*",
            resources = server.size["local_instance"]
    }

    # Compare quality control for all recipient files
    call quality.multi_qc as recipient_crossing_multiqc {
        input:
            quality_files = flatten(
                [
                crossing.bams,
                crossing_metrics.stats,
                crossing_metrics.flagstats
                ]),
            report_name = "${srr_name}-to-${recipient_name}_multiqc_metrics.html",
            enable_fullnames = false,
            include = "../inputs/*/*${recipient_name}*",
            resources = server.size["local_instance"]
    }

    output {
        Array[File] out_bucket_donor_bams = donor_bucketize.bams
        Array[File] out_bucket_donor_bais = donor_bucketize.bais
        Array[File] out_bucket_recipient_bams = recipient_bucketize.bams
        Array[File] out_bucket_recipient_bais = recipient_bucketize.bais

        Array[File] out_crossing_bams = crossing.bams
        Array[File] out_crossing_bais = indexing_bams.out
        Array[File] out_crossing_beds = bams_to_beds.out

        Array[File] out_donor_stats = donor_bucketized_metrics.stats
        Array[File] out_donor_flagstats = donor_bucketized_metrics.flagstats

        Array[File] out_recipient_stats = recipient_bucketized_metrics.stats
        Array[File] out_recipient_flagstats = recipient_bucketized_metrics.flagstats

        Array[File] out_crossing_stats = crossing_metrics.stats
        Array[File] out_crossing_flagstats = crossing_metrics.flagstats

        File? out_multiqc_donor_crossing_multiqc_report = donor_crossing_multiqc.out
        File? out_multiqc_recipient_crossing_multiqc_report = recipient_crossing_multiqc.out
    }
}