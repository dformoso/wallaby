version 1.0

import "structs/compute.wdl"

task download_srr {
    input {
        String srr
        Resources resources
    }

    command <<<
        fasterq-dump --split-files \
            ~{"--threads " + resources.cpu} \
            ~{srr}
    >>>

    output {
        File out_1 = "~{srr}_1.fastq"
        File out_2 = "~{srr}_2.fastq"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/sratoolkit:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
