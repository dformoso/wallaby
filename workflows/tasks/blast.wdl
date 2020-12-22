version 1.0

import "structs/compute.wdl"

task n {
    
    input {
        Array[File] fastas
        File blastdb
        Int evalue
        Resources resources
    }

    command <<<
        tar xfv ~{blastdb}
        export BLASTDB=`pwd`
        for fasta in ~{sep="  " fastas}
            do
                if 
                    [[ $fasta =~ "MMd_MUr" ]] ||
                    [[ $fasta =~ "MMd_UMr" ]] ||
                    [[ $fasta =~ "MUd_MUr" ]] ||
                    [[ $fasta =~ "MUd_UMr" ]] ||
                    [[ $fasta =~ "UMd_MUr" ]] ||
                    [[ $fasta =~ "UMd_UMr" ]]
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
        rm nt.*
        rm taxdb.*
    >>>

    output {
        File donor_MMd_MUr = "reads-to-donor_MMd_MUr.fasta.blastn"
        File donor_MUd_UMr = "reads-to-donor_MUd_UMr.fasta.blastn"
        File donor_UMd_MUr = "reads-to-donor_UMd_MUr.fasta.blastn"
        File recipient_MMd_MUr = "reads-to-recipient_MMd_MUr.fasta.blastn"
        File recipient_MUd_UMr = "reads-to-recipient_MUd_UMr.fasta.blastn"
        File recipient_UMd_MUr = "reads-to-recipient_UMd_MUr.fasta.blastn"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "ncbi/blast:2.10.1"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}