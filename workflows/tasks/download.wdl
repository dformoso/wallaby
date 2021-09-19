version 1.0

import "structs/compute.wdl"

task srr {
    input {
        String srr
        Resources resources
    }

    command <<<
        prefetch ~{srr}
        fasterq-dump -c 4000MB -b 4000MB -m 4000MB --split-files --threads 8  ~{srr}
        rm -rf ~{srr}/
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
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
