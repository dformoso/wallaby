version 1.0

import "donor_recipient.wdl"
import "tasks/download.wdl" as download

workflow multi_donor_recipient {

    input {
        String donor_name
        File donor_ref_genome
        String recipient_name
        File recipient_ref_genome
        File srr_list
        File blastdb
    }

    # Compute resources
    Compute server = read_json("../config/sizes.json")
    
    # Download all FASTQ files from the given SRR numbers
    Array[String] srrs = read_lines(srr_list)

    scatter (srr_name in srrs) {
        call download.srr as downloaded_srr {
            input:
                srr = srr_name,
                resources = server.size["2cpu_8mem_100disk"]
        }

        call donor_recipient.main as single_donor_recipient {
            input:
                donor_name = donor_name,
                donor_ref_genome = donor_ref_genome,
                recipient_name = recipient_name,
                recipient_ref_genome = recipient_ref_genome,
                srr_name = srr_name,
                srr_fastq_1 = downloaded_srr.out_1,
                srr_fastq_2 = downloaded_srr.out_2,
                blastdb = blastdb
        }
    }

    output {
        Array[File] out_donor_bams = single_donor_recipient.out_donor_bam
        Array[File] out_donor_bais = single_donor_recipient.out_donor_bai
        Array[File] out_recipient_bams = single_donor_recipient.out_recipient_bam
        Array[File] out_recipient_bais = single_donor_recipient.out_recipient_bai

        Array[Array[File]] out_bucket_donor_bams = single_donor_recipient.out_bucket_donor_bams
        Array[Array[File]] out_bucket_donor_bais = single_donor_recipient.out_bucket_donor_bais
        Array[Array[File]] out_bucket_recipient_bams = single_donor_recipient.out_bucket_recipient_bams
        Array[Array[File]] out_bucket_recipient_bais = single_donor_recipient.out_bucket_recipient_bais

        Array[Array[File]] out_crossed_filtered_bams = single_donor_recipient.out_crossed_filtered_bams
#        Array[Array[File]] out_crossed_filtered_fastas = single_donor_recipient.out_crossed_filtered_fastas
        Array[Array[File]] out_crossed_filtered_bais = single_donor_recipient.out_crossed_filtered_bais
        Array[Array[File]] out_crossed_filtered_beds = single_donor_recipient.out_crossed_filtered_beds

#        Array[File] out_donor_MMd_MUr = single_donor_recipient.out_donor_MMd_MUr
#        Array[File] out_donor_MUd_UMr = single_donor_recipient.out_donor_MUd_UMr
#        Array[File] out_donor_UMd_MUr = single_donor_recipient.out_donor_UMd_MUr
#        Array[File] out_recipient_MMd_MUr = single_donor_recipient.out_recipient_MMd_MUr
#        Array[File] out_recipient_MUd_UMr = single_donor_recipient.out_recipient_MUd_UMr
#        Array[File] out_recipient_UMd_MUr = single_donor_recipient.out_recipient_UMd_MUr

        Array[File] out_multiqc_before_and_after_trim_reports = single_donor_recipient.out_multiqc_before_and_after_trim_report
        Array[File] out_multiqc_all_donor_metrics_reports = single_donor_recipient.out_multiqc_donor_crossed_filtered_multiqc_report
        Array[File] out_multiqc_all_recipient_metrics_reports = single_donor_recipient.out_multiqc_recipient_crossed_filtered_multiqc_report
    }

}



