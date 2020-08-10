version development

import "structs/compute.wdl"

task view {
    input {
        File file
        String out_file
        Boolean include_header = false
        Boolean output_count = false
        Boolean output_file = true
        String include = ""
        String exclude = ""
        Resources resources
    }

    command <<<
        samtools view \
            ~{"-f " + include} \
            ~{"-F " + exclude} \
            ~{true="-h" false="" include_header} \
            ~{true="-b" false="" output_file} \
            ~{true="-c" false="" output_count} \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}


task sort {
    input {
        File file
        String out_file
        Resources resources
    }

    command <<<
        samtools sort \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}


task index {
    input {
        File file
        String out_file
        Resources resources
    }

    command <<<
        samtools index \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}


task merge {
    input {
        Array[File] files
        String out_file
        Resources resources
    }

    command <<<
        samtools merge \
            ~{"-@ " + resources.cpu} \
            ~{out_file} \
            ~{sep=" " files}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task mpileup {
    input {
        File fasta
        File bam
        String out_file
        Resources resources
    }

    command <<<
        samtools mpileup \
            ~{"-f " + fasta} \
            ~{bam} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task stats {
    input {
        File file
        String out_file
        Resources resources
    }

    command <<<
        set -e
        if [ -s "~{file}" ]
        then
        samtools stats \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            > ~{out_file}
        fi
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task flagstats {
    input {
        File file
        String out_file
        Resources resources
    }

    command <<<
        samtools flagstat \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task extract_qnames {
    input {
        File file
        String out_file
        Resources resources
    }

    command <<<
        samtools view \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            | cut -f1 | sort | uniq \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task count {
    input {
        File file
        String out_file
        Boolean include_header = false
        Boolean output_count = true
        Boolean output_file = false
        String include = ""
        String exclude = ""
        Resources resources
    }

    command <<<
        set -e
        if [ -s "~{file}" ]
        then
        samtools view \
            ~{"-f " + include} \
            ~{"-F " + exclude} \
            ~{true="-h" false="" include_header} \
            ~{true="-b" false="" output_file} \
            ~{true="-c" false="" output_count} \
            ~{"-@ " + resources.cpu} \
            ~{file} \
            > ~{out_file}
        fi
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task bam_to_fastas {
    input {
        File file
        String out_file_1
        String out_file_2
        Resources resources
    }

    command <<<
        samtools fasta \
            ~{"-1 " + out_file_1} \
            ~{"-2 " + out_file_2} \
            ~{"-0 /dev/null"} \
            ~{"-s /dev/null"} \
            ~{"-n " + file}
    >>>

    output {
        File out_1 = out_file_1
        File out_2 = out_file_2
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}