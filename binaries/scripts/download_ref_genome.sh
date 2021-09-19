#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

sudo apt-get install -y samtools

rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ../../data/ref_genomes/human/
gunzip ../../data/ref_genomes/human/hg38.fa.gz
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz ../../data/ref_genomes/human/
gunzip ../../data/ref_genomes/human/hg38.ncbiRefSeq.gtf.gz

samtools ../../data/ref_genomes/human/faidx hg38.fa 