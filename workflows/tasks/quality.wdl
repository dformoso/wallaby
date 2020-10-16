version 1.0

import "structs/compute.wdl"

task fast_qc {
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
        Array[File] files = [fastq_1_zip, fastq_2_zip, fastq_1_html, fastq_2_html]
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/quality:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task multi_qc {
    input {
        Array[File] quality_files
        String report_name = "multiqc_report.html"
        Resources resources
    }

    command <<<
        set -e
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        for quality_file in ~{sep=' ' quality_files}; do ls ${quality_file} .; done;
        multiqc -n ~{report_name} --fullnames ../inputs
    >>>

    output {
        File out = "~{report_name}"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/quality:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}