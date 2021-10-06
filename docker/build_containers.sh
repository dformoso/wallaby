#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

sudo apt install gnupg2 pass -y -qq
sudo docker login 
sudo docker login docker.io

# Build docker containers locally
sudo docker build --rm -t dformoso/albacore:latest albacore
sudo docker build --rm -t dformoso/nanopolish:latest nanopolish
sudo docker build --rm -t dformoso/nanopore-util:latest nanopore-util
sudo docker build --rm -t dformoso/guppy-gpu:latest guppy-gpu
sudo docker build --rm -t dformoso/deepbinner:latest deepbinner
sudo docker build --rm -t dformoso/minimap2:latest minimap2
sudo docker build --rm -t dformoso/samtools:latest samtools
sudo docker build --rm -t dformoso/porechop:latest porechop
sudo docker build --rm -t dformoso/trimmomatic:latest trimmomatic
sudo docker build --rm -t dformoso/sratoolkit:latest sratoolkit
sudo docker build --rm -t dformoso/bwa:latest bwa
sudo docker build --rm -t dformoso/bwa_samtools:latest bwa_samtools
sudo docker build --rm -t dformoso/star:latest star
sudo docker build --rm -t dformoso/aligners:latest aligners
sudo docker build --rm -t dformoso/picard-tools:latest picard-tools
sudo docker build --rm -t dformoso/prinseq:latest prinseq
sudo docker build --rm -t dformoso/quality:latest quality
sudo docker build --rm -t dformoso/jupyterlab:latest jupyterlab
sudo docker build --rm -t dformoso/bedops:latest bedops
sudo docker build --rm -t dformoso/rstudio:latest rstudio
sudo docker build --rm -t dformoso/seqtk:latest seqtk

# Push docker containers to docker hub
sudo docker push dformoso/albacore:latest
sudo docker push dformoso/nanopolish:latest
sudo docker push dformoso/nanopore-util:latest
sudo docker push dformoso/guppy-gpu:latest
sudo docker push dformoso/deepbinner:latest
sudo docker push dformoso/minimap2:latest
sudo docker push dformoso/samtools:latest
sudo docker push dformoso/porechop:latest
sudo docker push dformoso/sratoolkit:latest
sudo docker push dformoso/trimmomatic:latest
sudo docker push dformoso/bwa:latest
sudo docker push dformoso/bwa_samtools:latest
sudo docker push dformoso/star:latest
sudo docker push dformoso/aligners:latest
sudo docker push dformoso/picard-tools:latest
sudo docker push dformoso/prinseq:latest
sudo docker push dformoso/quality:latest
sudo docker push dformoso/jupyterlab:latest 
sudo docker push dformoso/bedops:latest 
sudo docker push dformoso/rstudio:latest 
sudo docker push dformoso/seqtk:latest 

