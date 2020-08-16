version development

import "structs/compute.wdl"

task qc {
    input {
        File fastq_1
        File fastq_2
        Resources resources
    }

    command <<<
        set -e
        mkdir temp
        fastqc \
            ~{"-threads " + resources.cpu} \
            ~{fastq_1} \
            ~{fastq_2} \
            -o `pwd`
    >>>

    output {
        File fastq_1_zip = "~{basename(fastq_1, ".fastq")}_fastqc.zip"
        File fastq_2_zip = "~{basename(fastq_2, ".fastq")}_fastqc.zip"
        File fastq_1_html = "~{basename(fastq_1, ".fastq")}_fastqc.html"
        File fastq_2_html = "~{basename(fastq_2, ".fastq")}_fastqc.html"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/fastqc:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
