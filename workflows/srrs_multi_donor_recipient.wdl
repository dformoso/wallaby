version 1.0

import "paired_donor_recipient.wdl"
import "tasks/align.wdl" as align
import "tasks/download.wdl" as download
import "tasks/trimmomatic.wdl" as trimmomatic
import "tasks/quality.wdl" as quality
import "tasks/tools.wdl" as tools
import "tasks/structs/compute.wdl"

workflow multi_donor_recipient {

    input {
        File donor_ref_genome
        File donor_ref_genome_fai
        File donor_ref_genome_gff  # .gff must exist, but it can be an empty file
        File recipient_ref_genome
        File srr_list
    }

    String donor_name = basename(donor_ref_genome, ".fa")
    String recipient_name = basename(recipient_ref_genome, ".fa")

    # Compute resources
    Compute server = read_json("../inputs/sizes.json")

    # Indexing Donor and Recipient Reference Genomes
    call align.index as donor_index { 
        input: 
            fasta = donor_ref_genome,
            resources = server.size["local_instance"]
    }

    call align.index as recipient_index { 
        input: 
            fasta = recipient_ref_genome,
            resources = server.size["local_instance"]
    }

    # Download all FASTQ files from the given SRR numbers
    Array[String] srrs = read_lines(srr_list)
    
    scatter (srr_name in srrs) {

        # Download SSRs
        call download.srr as downloaded_srr {
            input:
                srr = srr_name,
                sample = "false",
                sampling_factor = 1,
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
        call quality.multi_qc as srr_multiqc_trim {
            input:
                quality_files = flatten(
                    [
                    srr_fastqc_before_trim.files,
                    srr_fastqc_after_trim.files
                    ]),
                report_name = "${srr_name}_multiqc_report.html",
                include = "../inputs/*",            
                resources = server.size["local_instance"]
        }

        # Launch Donor to Recipient main workflow 
        call paired_donor_recipient.main as donor_recipient {
            input:
                donor_name = donor_name,
                donor_index = donor_index.index_object,
                donor_ref_genome = donor_ref_genome,
                recipient_name = recipient_name,
                recipient_index = recipient_index.index_object,
                recipient_ref_genome = recipient_ref_genome,
                srr_name = srr_name,
                fastq_1 = srr_trim_adapters.fastq_1_paired,
                fastq_2 = srr_trim_adapters.fastq_2_paired
        }
    }

    # Merge all putative insertion tables into one .csv file
    call tools.merge_csvs as overlap_loci_table {
        input:
            csvs = select_all(donor_recipient.overlap_loci_csv),
            resources = server.size["local_instance"]
    }

    # Adding files neede for the downstream visualization tools
    call tools.summary_and_inputs as summary_and_inputs {
        input:
            donor_name = donor_name,
            donor_ref_genome = donor_ref_genome,
            donor_ref_genome_fai = donor_ref_genome_fai,
            donor_ref_genome_gff = donor_ref_genome_gff,
            recipient_name = recipient_name,
            resources = server.size["local_instance"]
    }

    output {
        File out_summary = summary_and_inputs.summary
        File out_donor_ref_genome = summary_and_inputs.out_donor_ref_genome
        File out_donor_ref_genome_fai = summary_and_inputs.out_donor_ref_genome_fai
        File out_donor_ref_genome_gff = summary_and_inputs.out_donor_ref_genome_gff

        Array[Array[File]] out_filtered_bams = donor_recipient.filtered_bams
        Array[Array[File]] out_filtered_bais = donor_recipient.filtered_bais
        Array[Array[File]] out_filtered_beds = donor_recipient.filtered_beds
        Array[Array[File]] out_overlap_loci_bams_and_bais = donor_recipient.overlap_loci_bams_and_bais
        
        File out_overlap_loci_table = overlap_loci_table.out

        Array[File?] out_multiqc_trim_html = select_all(srr_multiqc_trim.html)
        Array[File?] out_multiqc_trim_zip = select_all(srr_multiqc_trim.zip)
        Array[File?] out_multiqc_donor_html = select_all(donor_recipient.multiqc_donor_html)
        Array[File?] out_multiqc_donor_zip = select_all(donor_recipient.multiqc_donor_zip)
        Array[File?] out_multiqc_recipient_html = select_all(donor_recipient.multiqc_recipient_html)
        Array[File?] out_multiqc_recipient_zip = select_all(donor_recipient.multiqc_recipient_zip)
    }
}