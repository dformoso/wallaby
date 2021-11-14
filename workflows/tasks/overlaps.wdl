version 1.0

import "structs/compute.wdl"

task putative_insertions {
    
    input {
        Array[File] bams
        Array[File] bais
        Array[File] beds

        String srr_name

        String donor_name
        File donor_ref_genome
        String recipient_name
        File recipient_ref_genome

        Int min_num_crossings = 1
        Int min_num_reads = 5

        Resources resources
    }

    command <<<
        # Link files to current directory
        for bam in ~{sep=' ' bams}
        do
            ln $bam .
        done
        for bai in ~{sep=' ' bais}
        do
            ln $bai .
        done
        for bed in ~{sep=' ' beds}
        do
            ln $bed .
        done

        # Run R Script to create .csv file with overlaps
        Rscript /tmp/analyze_srr.r \
            ~{srr_name} \
            ~{donor_name} \
            ~{recipient_name} \
            ~{donor_ref_genome} \
            ~{recipient_ref_genome} \
            "." \
            ~{min_num_crossings} \
            ~{min_num_reads}
        
        # Run Python Script to user .csv files to filter files using overlaps
        python3 /tmp/filter_srr.py \
            "~{srr_name}_putative_insertion_table.csv" \
            ~{srr_name} \
            ~{donor_name} \
            ~{recipient_name}

        ls *id*.bam
    >>>

    output {
        Array[File] out = read_lines(stdout())
        File csv = "~{srr_name}_putative_insertion_table.csv"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/overlaps:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
