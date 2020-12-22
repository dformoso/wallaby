version 1.0

import "structs/compute.wdl"

task filter_reads {
    input {
        File file
        File txt
        String out_file
        String filter = "includeReadList"
        Resources resources
    }

    command <<<
        set -e
        if [ -s "~{txt}" ]
        then
            if [ -s "~{file}" ] 
            then 
            java -jar /usr/picard/picard.jar FilterSamReads \
                ~{"INPUT=" + file} \
                ~{"OUTPUT=" + out_file} \
                ~{"READ_LIST_FILE=" + txt} \
                ~{"FILTER=" + filter}
            fi
        fi
    >>>

    output {
        File? out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/picard-tools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task filter_valid_reads {
    input {
        File file
        File txt
        String out_file
        Resources resources
        String filter = "includeReadList"
    }

    command {
        set -e
        if [ -s "~{txt}" ]
        then
            if [ -s "~{file}" ] 
            then 
            java -jar /usr/picard/picard.jar FilterSamReads \
                ~{"INPUT=" + file} \
                ~{"OUTPUT=" + out_file} \
                ~{"READ_LIST_FILE=" + txt} \
                ~{"FILTER=" + filter}
            fi
        fi
    }

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/picard-tools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}


task sam_to_fastq {
    input {
        File file
        String out_fastq_1
        String out_fastq_2
        Resources resources
        String validation_stringency = "SILENT"
    }

    command {
        set -e
        if [ -s "~{file}" ] 
        then 
        java -Xmx40G -jar /usr/picard/picard.jar SamToFastq \
            ~{"INPUT=" + file} \
            ~{"FASTQ=" + out_fastq_1} \
            ~{"SECOND_END_FASTQ=" + out_fastq_2} \
            ~{"VALIDATION_STRINGENCY=" + validation_stringency} 
        fi
    }

    output {
        File fastq_1 = out_fastq_1
        File fastq_2 = out_fastq_2
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/picard-tools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}