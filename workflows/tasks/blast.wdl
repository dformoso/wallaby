version development

task n {
    
    input {
        File fasta 
        Directory blastdb
        String out_file = "~{basename(fasta)}.blast"
    }

    command <<<
        export BLASTDB=~{blastdb} 
        blastn \
            -query ~{fasta} \
            -db nt \
            -num_threads 24 \
            -evalue 1 \
            -outfmt '6 seqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus' \
            -out ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "128GB"
        docker: "ncbi/blast:latest"
        disks: "local-disk 1000GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}