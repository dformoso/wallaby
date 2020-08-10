# Instructions on setting up NCBI BLAST locally

cd /media/formoso/NVMe1/research/genomics/software/

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
    tar xzf ncbi-blast-2.10.1+-x64-linux.tar.gz

# or permanently add the following to /etc/bash.bashrc
sudo nano /etc/bash.bashrc
PATH=$PATH:/media/formoso/NVMe1/research/genomics/software/ncbi-blast-2.10.1+/bin
export BLASTDB=/media/formoso/NVMe1/research/genomics/projects/wallaby/data/blast/blastdb


