version 1.0

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
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}


task index {
    input {
        File file
        String out_file = "~{basename(file)}.bai"
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
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task merge_if_exists {
    input {
        File? bam1
        File? bam2
        String out_file
        Resources resources
    }

    command <<<
        if [ -f "~{bam1}" ] ||  [ -f "~{bam2}" ] ; then
            samtools merge \
                ~{"-@ " + resources.cpu} \
                ~{out_file} \
                ~{bam1} ~{bam2}
        fi
    >>>

    output {
        File? out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
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
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task bam_to_fasta {
    input {
        File? file
        String out_file
        Resources resources
    }

    command <<<
        if test -f ${file}; then
            samtools fasta \
            ~{"-@ " + resources.cpu} \
            ~{"-n " + file} \
            > ~{out_file}
        fi
    >>>

    output {
        File? out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task bam_to_bed {
    input {
        File file
        String out_file = "~{basename(file, ".bam")}.bed"
        Resources resources
    }

    command <<<
        bedtools bamtobed -i \
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
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}