version development

import "structs/bwa.wdl"
import "structs/compute.wdl"

task index {
    
    input {
        File fasta 
        String basename_fasta = basename(fasta) 
        Resources resources
    }

    command <<<
        ln ~{fasta} .
        bwa index ~{basename_fasta} 
    >>>

    output {
        BWAIndex index = {
            "fasta": "~{basename_fasta}",
            "amb": "~{basename_fasta}.amb",
            "ann": "~{basename_fasta}.ann",
            "bwt": "~{basename_fasta}.bwt",
            "pac": "~{basename_fasta}.pac",
            "sa": "~{basename_fasta}.sa"
        }
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bwa:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task align {

    input {
        File fastq_1
        File fastq_2
        BWAIndex bwa_index
        String out_file = "reads-to-ref-genome.sam"
        Resources resources
    }

    command <<<
        bwa mem \
            ~{"-t " + resources.cpu} \
            ~{bwa_index.fasta} \
            ~{fastq_1} \
            ~{fastq_2} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bwa:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
