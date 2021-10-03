version 1.0

import "structs/compute.wdl"

task srr {
    input {
        String srr
        Float sampling_rate = 1
        Resources resources
    }

    command <<<
        prefetch ~{srr} --max-size 100G
        fasterq-dump -c 4000MB -b 4000MB -m 4000MB --split-files --threads 8  ~{srr}
        seqtk sample -s 11 "~{srr}_1.fastq" ~{sampling_rate} > "~{srr}_1_sampled.fastq"
        seqtk sample -s 11 "~{srr}_2.fastq" ~{sampling_rate} > "~{srr}_2_sampled.fastq"
        rm -rf ~{srr}/
        rm -rf "~{srr}_1.fastq"
        rm -rf "~{srr}_2.fastq"
        mv "~{srr}_1_sampled.fastq" "~{srr}_1.fastq"
        mv "~{srr}_2_sampled.fastq" "~{srr}_2.fastq"
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
