version development

import "structs/structures.wdl"

task index {
    
    input {
        File fasta 
        String basename_fasta = basename(fasta) 
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
        cpu: "4"
        memory: "24GB"
        docker: "dformoso/bwa:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

task align {

    input {
        File fastq_1
        File fastq_2
        BWAIndex bwa_index
        String out_file = "reads-to-ref-genome.sam"
        Int bwa_threads = 24
    }

    command <<<
        bwa mem \
            ~{"-t " + bwa_threads} \
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
        cpu: bwa_threads
        memory: "16GB"
        docker: "dformoso/bwa:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}
