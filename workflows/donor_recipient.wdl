version development

import "subworkflows/align.wdl" as align
import "subworkflows/bucketize.wdl" as bucketize
import "subworkflows/cross.wdl" as cross
import "subworkflows/complex.wdl" as complex
import "subworkflows/blast.wdl" as blast
import "subworkflows/locate.wdl" as locate
import "subworkflows/metrics.wdl" as metrics

workflow donor_recipient {

    input {
        File donor_ref_genome
        File recipient_ref_genome
        File reads_fastq_1
        File reads_fastq_2
        Directory blastdb
    }

    # Donor Reference Genome
    ## Align 
    call align.main as align_donor {
        input:
            ref_genome = donor_ref_genome,
            fastq_1 = reads_fastq_1,
            fastq_2 = reads_fastq_2,
            base_filename = "reads-to-donor"
    }

    ## Bucketize
    call bucketize.main as bucketize_donor {
        input:
            index = align_donor.index,
            bam = align_donor.bam,
            base_filename = "reads-to-donor"
    }

    # Recipient Reference Genome
    ## Align 
    call align.main as align_recipient {
        input:
            ref_genome = recipient_ref_genome,
            fastq_1 = reads_fastq_1,
            fastq_2 = reads_fastq_2,
            base_filename = "reads-to-recipient"
    }

    ## Bucketize
    call bucketize.main as bucketize_recipient {
        input:
            index = align_recipient.index,
            bam = align_recipient.bam,
            base_filename = "reads-to-recipient"
    }

    # Cross BAM files buckets
    call cross.main as cross {
        input:
            donor_MM = bucketize_donor.MM,
            donor_MU = bucketize_donor.MU,
            donor_UM = bucketize_donor.UM,
            donor_UU = bucketize_donor.UU,
            recipient_MM = bucketize_recipient.MM,
            recipient_MU = bucketize_recipient.MU,
            recipient_UM = bucketize_recipient.UM,
            recipient_UU = bucketize_recipient.UU
    }

    # Filter out low complexity sequences
    call complex.main as complex {
        input:
            bam_files = cross.bams
    }

    # Locate sequences in reference genome
    call locate.main as locate_donor {
        input:
            bams = complex.bams,
            ref_genome = donor_ref_genome
    }

    call locate.main as locate_recipient {
        input:
            bams = complex.bams,
            ref_genome = recipient_ref_genome
    }

    # Blastn search over all crossings of interest
    call blast.main as blaster {
        input:
            fastas = complex.fastas,
            blastdb = blastdb
    }

    # Calculate metrics
    call metrics.main as metrics_donor {
        input:
            bams = bucketize_donor.bams
    }

    call metrics.main as metrics_recipient {
        input:
            bams = bucketize_recipient.bams
    }

    call metrics.main as metrics_cross {
        input:
            bams = cross.bams
    }
    
    call metrics.main as metrics_complex {
        input:
            bams = complex.bams
    }

    output {
        Array[File] out_donor_bams = bucketize_donor.bams
        Array[File] out_recipient_bams = bucketize_recipient.bams
        Array[File] out_cross_bams = cross.bams
        Array[File] out_complex_bams = complex.bams
        Array[File] out_fastas = complex.fastas
        Array[File] out_blastns = blaster.blastns

        Array[File] donor_stats = metrics_donor.stats
        Array[File] donor_flagstats = metrics_donor.flagstats
        Array[File] donor_count = metrics_donor.count

        Array[File] recipient_stats = metrics_recipient.stats
        Array[File] recipient_flagstats = metrics_recipient.flagstats
        Array[File] recipient_count = metrics_recipient.count

        Array[File] cross_stats = metrics_cross.stats
        Array[File] cross_flagstats = metrics_cross.flagstats
        Array[File] cross_count = metrics_cross.count

        Array[File] complex_stats = metrics_complex.stats
        Array[File] complex_flagstats = metrics_complex.flagstats
        Array[File] complex_count = metrics_complex.count
    }

}