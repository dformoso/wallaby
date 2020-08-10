version development

import "structs/compute.wdl"

task n {
    
    input {
        Array[File] fastas
        Directory blastdb
        Resources resources
        String evalue
    }

    command <<<
        export BLASTDB=~{blastdb} 
        for fasta in ~{sep="  " fastas}
            do
                if 
                    [[ $fasta =~ "donor_MMd_MUr" ]] ||
                    [[ $fasta =~ "donor_MMd_UMr" ]] ||
                    [[ $fasta =~ "donor_MUd_MUr" ]] ||
                    [[ $fasta =~ "donor_MUd_UMr" ]] ||
                    [[ $fasta =~ "donor_UMd_MUr" ]] ||
                    [[ $fasta =~ "donor_UMd_UMr" ]] ||
                    [[ $fasta =~ "recipient_MMd_MUr" ]] ||
                    [[ $fasta =~ "recipient_MMd_UMr" ]] ||
                    [[ $fasta =~ "recipient_MUd_MUr" ]] ||
                    [[ $fasta =~ "recipient_MUd_UMr" ]] ||
                    [[ $fasta =~ "recipient_UMd_MUr" ]] ||
                    [[ $fasta =~ "recipient_UMd_UMr" ]]
                then
                    blastn \
                        -query $fasta \
                        -db nt \
                        -num_threads ~{resources.cpu} \
                        -evalue ~{evalue} \
                        -outfmt '6 seqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus' \
                        -out "`basename ${fasta}`.blastn"
                fi
            done 
    >>>

    output {
        Array[File] out = glob("*.blastn")
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "ncbi/blast:2.10.1"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}