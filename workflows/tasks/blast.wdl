version development

task n {
    
    input {
        Array[File] fastas
        Directory blastdb
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
                        -num_threads 24 \
                        -evalue 1 \
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
        cpu: "24"
        memory: "128GB"
        docker: "ncbi/blast:2.10.1"
        disks: "local-disk 1000GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}