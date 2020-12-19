docker run -it --rm -v "/root/research/genomics/projects/wallaby/data:/tmp/" dformoso/sratoolkit

mkdir SRR12091993srr
mkdir SRR12091994
mkdir SRR12091995

mkdir SRR12091996
mkdir SRR12091997
mkdir SRR12091998


#HPV
fasterq-dump --split-files --threads 8 SRR5090597
fasterq-dump --split-files --threads 8 SRR5090598
fasterq-dump --split-files --threads 8 SRR5090599
fasterq-dump --split-files --threads 8 SRR5090600


# SARS-COV-2
fasterq-dump --split-files --threads 8 SRR12091993 
fasterq-dump --split-files --threads 8 SRR12091994
fasterq-dump --split-files --threads 8 SRR12091995

fasterq-dump --split-files --threads 8 SRR12091996
fasterq-dump --split-files --threads 8 SRR12091997
fasterq-dump --split-files --threads 8 SRR12091998



docker run -it --rm -v "/root/research/genomics/projects/wallaby/workflows/cromwell-final-outputs:/tmp/" dformoso/samtools

