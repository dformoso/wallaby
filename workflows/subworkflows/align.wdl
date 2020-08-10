version development

import "../tasks/bwa.wdl" as bwa
import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/bwa.wdl"
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        File ref_genome
        File fastq_1
        File fastq_2
        String base_filename
        Resources resources
    }

    # Aligning FASTQ files to Reference Genomes
    call bwa.index as bwa_indexer { 
        input: 
            fasta = ref_genome,
            resources = resources

    }

    call bwa.align as bwa_aligner { 
        input: 
            bwa_index = bwa_indexer.index, 
            fastq_1 = fastq_1, 
            fastq_2 = fastq_2, 
            out_file = "~{base_filename}.sam",
            resources = resources
    }

    # Creating a BAM indexed file
    call samtools.view as samtools_sam_to_bam { 
        input: 
            file = bwa_aligner.out, 
            out_file = "~{base_filename}.bam",
            resources = resources
    }
    
    call samtools.sort as samtools_sorter { 
        input: 
            file = samtools_sam_to_bam.out, 
            out_file = "~{base_filename}.sorted.bam",
            resources = resources
    }

    call samtools.index as samtools_indexer { 
        input: 
            file = samtools_sorter.out, 
            out_file = "~{base_filename}.sorted.bam.bai",
            resources = resources
    }

    output {
        BWAIndex index = bwa_indexer.index
        File bam = samtools_sorter.out
    }

}
