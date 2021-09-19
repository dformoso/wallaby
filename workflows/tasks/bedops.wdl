version 1.0

import "structs/compute.wdl"

task gtf2bed {
    input {
        File gtf
        Resources resources
    }
    
    String out_file = "~{basename(gtf, ".gtf")}.bed"

    command <<<
        gtf2bed \
            < ~{gtf} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bedops:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task sortbed {
    input {
        File bed
        Resources resources
    }
        
    String out_file = "~{basename(bed, ".bed")}_sorted.bed"

    command <<<
        sort-bed \
            ~{bed} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bedops:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task bedmap {
    input {
        File bed
        File tax_bed
        Resources resources
    }
        
    String out_file = "~{basename(bed, "_sorted.bed")}_tax.bed"

    command <<<
        bedmap --echo --echo-map-id-uniq --delim '\t' \
            ~{bed} \
            ~{tax_bed} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bedops:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}