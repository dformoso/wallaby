version 1.0

import "paired_donor_recipient.wdl"
import "tasks/align.wdl" as align
import "tasks/download.wdl" as download
import "tasks/trimmomatic.wdl" as trimmomatic
import "tasks/quality.wdl" as quality
import "tasks/structs/compute.wdl"

workflow multi_donor_recipient {

    input {
        String donor_name
        File donor_ref_genome
        String recipient_name
        File recipient_ref_genome
        File srr_list
        String aligner_type
    }

    # Compute resources
    Compute server = read_json("../config/sizes.json")
    
    # Download all FASTQ files from the given SRR numbers
    Array[String] srrs = read_lines(srr_list)

    # Indexing Donor and Recipient Reference Genomes
    call align.index as donor_index { 
        input: 
            fasta = donor_ref_genome,
            aligner_type = aligner_type,
            resources = server.size["local_instance"]
    }

    call align.index as recipient_index { 
        input: 
            fasta = recipient_ref_genome,
            aligner_type = aligner_type,
            resources = server.size["local_instance"]
    }

    scatter (srr_name in srrs) {

        # Download SSRs
        call download.srr as downloaded_srr {
            input:
                srr = srr_name,
                resources = server.size["local_instance"]
        }

        # SRRs Quality Control metrics
        call quality.fast_qc as srr_fastqc_before_trim {
            input:
                fastq_1 = downloaded_srr.out_1,
                fastq_2 = downloaded_srr.out_2,
                resources = server.size["local_instance"]
        }

        # Quality trim the SRR reads
        call trimmomatic.trim as srr_trim_adapters {
            input:
                fastq_1 = downloaded_srr.out_1,
                fastq_2 = downloaded_srr.out_2,
                adapter = "all_adapters.fa",
                seed_mismatches = 2,
                paired_clip_threshold = 30,
                unpaired_clip_threshold = 10,
                leading = 3,
                trailing = 3,
                sliding_window_quality = 20,
                sliding_window_length = 4,
                min_length = 10,
                is_phred33 = true,
                resources = server.size["local_instance"]
        }

        # SRRs Quality control metrics after trimmomatic
        call quality.fast_qc as srr_fastqc_after_trim {
            input:
                fastq_1 = srr_trim_adapters.fastq_1_paired,
                fastq_2 = srr_trim_adapters.fastq_2_paired,
                resources = server.size["local_instance"]
        }

        # Compare quality control SRRs metrics before and after trimmomatic
        call quality.multi_qc as srr_multiqc_after_trim {
            input:
                quality_files = flatten(
                    [
                    srr_fastqc_before_trim.files,
                    srr_fastqc_after_trim.files
                    ]),
                report_name = "${srr_name}_multiqc_trim_report.html",
                include = "../inputs/*",            
                resources = server.size["local_instance"]
        }

        # Launch Donor to Recipient main workflow 
        call paired_donor_recipient.main as donor_recipient {
            input:
                donor_name = donor_name,
                donor_index = donor_index.index_object,
                recipient_name = recipient_name,
                recipient_index = recipient_index.index_object,
                srr_name = srr_name,
                fastq_1 = srr_trim_adapters.fastq_1_paired,
                fastq_2 = srr_trim_adapters.fastq_2_paired,
                aligner_type = aligner_type
        }
    }

    output {        
        Array[File] out_pre_fastq_1_zip = srr_fastqc_before_trim.fastq_1_zip
        Array[File] out_pre_fastq_2_zip = srr_fastqc_before_trim.fastq_2_zip
        Array[File] out_pre_fastq_1_html = srr_fastqc_before_trim.fastq_1_html
        Array[File] out_pre_fastq_2_html = srr_fastqc_before_trim.fastq_2_html

        Array[File] out_post_fastq_1_zip = srr_fastqc_after_trim.fastq_1_zip
        Array[File] out_post_fastq_2_zip = srr_fastqc_after_trim.fastq_2_zip
        Array[File] out_post_fastq_1_html = srr_fastqc_after_trim.fastq_1_html
        Array[File] out_post_fastq_2_html = srr_fastqc_after_trim.fastq_2_html

        Array[Array[File]] out_bucket_donor_bams = donor_recipient.out_bucket_donor_bams
        Array[Array[File]] out_bucket_donor_bais = donor_recipient.out_bucket_donor_bais
        Array[Array[File]] out_bucket_recipient_bams = donor_recipient.out_bucket_recipient_bams
        Array[Array[File]] out_bucket_recipient_bais = donor_recipient.out_bucket_recipient_bais

        Array[Array[File]] out_crossing_bams = donor_recipient.out_crossing_bams
        Array[Array[File]] out_crossing_bais = donor_recipient.out_crossing_bais
        Array[Array[File]] out_crossing_beds = donor_recipient.out_crossing_beds

        Array[File?] out_multiqc_before_and_after_trim_report = select_all(srr_multiqc_after_trim.out)
        Array[File?] out_multiqc_all_donor_metrics_reports = select_all(donor_recipient.out_multiqc_donor_crossing_multiqc_report)
        Array[File?] out_multiqc_all_recipient_metrics_reports = select_all(donor_recipient.out_multiqc_recipient_crossing_multiqc_report)
    }
}