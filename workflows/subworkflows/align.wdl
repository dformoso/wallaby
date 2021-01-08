version 1.0

import "../tasks/align.wdl" as align
import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        String aligner_type 
        File fastq_1
        File fastq_2
        Array[File] index_object
        String base_filename
        Resources resources
    }

    # Aligning FASTQ files to Reference Genomes
    call align.align as aligner { 
        input: 
            index_object = index_object, 
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            aligner_type = aligner_type,
            out_file = "~{base_filename}.sam",
            resources = resources
    }

    # Creating a BAM indexed file
    call samtools.view as samtools_sam_to_bam { 
        input: 
            file = aligner.out, 
            out_file = "~{base_filename}.bam",
            resources = resources
    }
    
    call samtools.sort as samtools_sorter { 
        input: 
            file = samtools_sam_to_bam.out, 
            out_file = "~{base_filename}.bam",
            resources = resources
    }

    call samtools.index as samtools_indexer { 
        input: 
            file = samtools_sorter.out, 
            out_file = "~{base_filename}.bai",
            resources = resources
    }

    output {
        File bam = samtools_sorter.out
        File bai = samtools_indexer.out
    }

}