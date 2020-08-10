#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

cd /media/formoso/NVMe1/research/genomics/software/

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    tar xzf ncbi-blast-2.10.1+-x64-linux.tar.gz

PATH="/media/formoso/NVMe1/research/genomics/software/ncbi-blast-2.10.1+/bin:${PATH}"
export BLASTDB=/media/formoso/NVMe1/research/genomics/projects/wallaby/data/blast/blastdb


# or permanently add the following to /etc/bash.bashrc
sudo nano /etc/bash.bashrc
PATH=$PATH:/media/formoso/NVMe1/research/genomics/software/ncbi-blast-2.10.1+/bin
export BLASTDB=/media/formoso/NVMe1/research/genomics/projects/wallaby/data/blast/blastdb



blastn \
    -query /media/formoso/NVMe1/research/genomics/projects/wallaby/data/ref_genomes/pseudomonas_aeruginosa_pao1.fasta \
    -db nt \
    -num_threads 24 \
    -evalue 1 \
    -outfmt '6 seqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus'