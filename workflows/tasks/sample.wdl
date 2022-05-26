version 1.0

import "structs/compute.wdl"

task fastq {
    
    input {
        File fastq_1
        File fastq_2
        Float sampling = 1
        Resources resources
    }

    command <<<
        seqtk sample -s 11 ~{fastq_1} ~{sampling} > "~{basename(fastq_1)}"
        seqtk sample -s 11 ~{fastq_2} ~{sampling} > "~{basename(fastq_2)}"
    >>>

    output {
        File out_1 = "~{basename(fastq_1)}"
        File out_2 = "~{basename(fastq_2)}"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/seqtk:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}