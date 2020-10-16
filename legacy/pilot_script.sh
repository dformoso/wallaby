# VARIABLES
PIPELINE_FOLDER='/media/formoso/NVMe/research/genomics/pipelines/lgtsearch'
SUBJECT='SRR5377828'
DONOR='pseudomonas_aeruginosa_pao1'
RECIPIENT='GRCh38_p13'
THREADS=24

# Download data from SRA
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/sratoolkit \
    fasterq-dump --split-files $SUBJECT

# Pseudomonas Aeruginosa
# Use BWA to align reads
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data:/tmp/" dformoso/bwa \
    bwa index ref_genomes/${DONOR}.fasta

sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data:/tmp/" dformoso/bwa \
    /bin/bash -c \
    " \
    bwa mem -t ${THREADS} \
        ref_genomes/${DONOR}.fasta \
        ${SUBJECT}/${SUBJECT}_1.fastq \
        ${SUBJECT}/${SUBJECT}_2.fastq > \
        alignments/${SUBJECT}-${DONOR}.sam 
    "

# SAM to BAM conversion and indexing
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -b -S alignments/${SUBJECT}-${DONOR}.sam > alignments/${SUBJECT}-${DONOR}.bam && \
    samtools sort -@ ${THREADS} alignments/${SUBJECT}-${DONOR}.bam -o alignments/${SUBJECT}-${DONOR}-sorted.bam && \
    samtools index -@ ${THREADS} alignments/${SUBJECT}-${DONOR}-sorted.bam 
    "


# GRCh38_p13
# Use BWA to align reads
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data:/tmp/" dformoso/bwa \
    bwa index ref_genomes/${RECIPIENT}.fasta

sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data:/tmp/" dformoso/bwa \
    /bin/bash -c \
    " \
    bwa mem -t ${THREADS} \
        ref_genomes/${RECIPIENT}.fasta \
        ${SUBJECT}/${SUBJECT}_1.fastq \
        ${SUBJECT}/${SUBJECT}_2.fastq > \
        alignments/${SUBJECT}-${RECIPIENT}.sam 
    "
    
# SAM to BAM conversion and indexing
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -b -S alignments/${SUBJECT}-${RECIPIENT}.sam > alignments/${SUBJECT}-${RECIPIENT}.bam && \
    samtools sort -@ ${THREADS} alignments/${SUBJECT}-${RECIPIENT}.bam -o alignments/${SUBJECT}-${RECIPIENT}-sorted.bam && \
    samtools index -@ ${THREADS} alignments/${SUBJECT}-${RECIPIENT}-sorted.bam 
    "

# Bucketise 
# Recipient -> Donor (Pseudomonas)
# R1 Mapped - R2 Mapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    samtools view -@ ${THREADS} -f 3 -F 2048 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_MM.bam

# R1 Mapped - R2 Unmapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -f 73 -F 4022 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_MU_R1.bam && \
    samtools view -@ ${THREADS} -f 133 -F 3962 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_MU_R2.bam && \
    samtools merge -f \
        buckets/${SUBJECT}-${DONOR}_MU.bam \
        -@ ${THREADS} \
        buckets/${SUBJECT}-${DONOR}_MU_R1.bam \
        buckets/${SUBJECT}-${DONOR}_MU_R2.bam 
    "
    
# R1 Unmapped - R2 Mapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -f 69 -F 4026 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_UM_R1.bam && \
    samtools view -@ ${THREADS} -f 137 -F 3958 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_UM_R2.bam && \
    samtools merge -f \
        buckets/${SUBJECT}-${DONOR}_UM.bam \
        -@ ${THREADS} \
        buckets/${SUBJECT}-${DONOR}_UM_R1.bam \
        buckets/${SUBJECT}-${DONOR}_UM_R2.bam 
    "

# R1 Unmapped - R2 Unmapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -f 77 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_UU_R1.bam && \
    samtools view -@ ${THREADS} -f 141 \
        alignments/${SUBJECT}-${DONOR}-sorted.bam \
        -o buckets/${SUBJECT}-${DONOR}_UU_R2.bam && \
    samtools merge -f \
        buckets/${SUBJECT}-${DONOR}_UU.bam \
        -@ ${THREADS} \
        buckets/${SUBJECT}-${DONOR}_UU_R1.bam \
        buckets/${SUBJECT}-${DONOR}_UU_R2.bam 
    "

# Subject -> Recipient (Human)
# R1 Mapped - R2 Mapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    samtools view -@ ${THREADS} -f 3 -F 2048 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_MM.bam

# R1 Mapped - R2 Unmapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -f 73 -F 4022 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_MU_R1.bam && \
    samtools view -@ ${THREADS} -f 133 -F 3962 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_MU_R2.bam && \
    samtools merge -f \
        buckets/${SUBJECT}-${RECIPIENT}_MU.bam \
        -@ ${THREADS} \
        buckets/${SUBJECT}-${RECIPIENT}_MU_R1.bam \
        buckets/${SUBJECT}-${RECIPIENT}_MU_R2.bam 
    "

# R1 Unmapped - R2 Mapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -f 69 -F 4026 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_UM_R1.bam && \
    samtools view -@ ${THREADS} -f 137 -F 3958 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_UM_R2.bam && \
    samtools merge -f \
        buckets/${SUBJECT}-${RECIPIENT}_UM.bam \
        -@ ${THREADS} \
        buckets/${SUBJECT}-${RECIPIENT}_UM_R1.bam \
        buckets/${SUBJECT}-${RECIPIENT}_UM_R2.bam 
    "

# R1 Unmapped - R2 Unmapped
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -f 77 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_UU_R1.bam && \
    samtools view -@ ${THREADS} -f 141 \
        alignments/${SUBJECT}-${RECIPIENT}-sorted.bam \
        -o buckets/${SUBJECT}-${RECIPIENT}_UU_R2.bam && \
    samtools merge -f \
        buckets/${SUBJECT}-${RECIPIENT}_UU.bam \
        -@ ${THREADS} \
        buckets/${SUBJECT}-${RECIPIENT}_UU_R1.bam \
        buckets/${SUBJECT}-${RECIPIENT}_UU_R2.bam 
    "

# Count of reads per bucket
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${DONOR}_MM.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${DONOR}_MU.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${DONOR}_UM.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${DONOR}_UU.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${RECIPIENT}_MM.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${RECIPIENT}_MU.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${RECIPIENT}_UM.bam ; \
    samtools view -@ ${THREADS} -c buckets/${SUBJECT}-${RECIPIENT}_UU.bam
    "

# Grouping between buckets
# Create files containing a unique, ordered, list of QNAMES (the identifiers for the reads, or Read IDs)
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${DONOR}_MU.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${DONOR}_MU_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${DONOR}_UM.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${DONOR}_UM_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${DONOR}_MM.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${DONOR}_MM_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${DONOR}_UU.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${DONOR}_UU_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${RECIPIENT}_MU.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${RECIPIENT}_MU_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${RECIPIENT}_UM.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${RECIPIENT}_UM_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${RECIPIENT}_MM.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${RECIPIENT}_MM_QNAMEs.txt ; \
    samtools view -@ ${THREADS} buckets/${SUBJECT}-${RECIPIENT}_UU.bam | cut -f1 | sort | uniq > bucket_pairs/${SUBJECT}-${RECIPIENT}_UU_QNAMEs.txt 
    "

# Create files containing the inner join between two groups of files, using QNAMEs
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    comm -12 \
        bucket_pairs/${SUBJECT}-${DONOR}_MU_QNAMEs.txt \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_UM_QNAMEs.txt \
        > bucket_pairs/donor_MU_UM_recipient.txt && \
    comm -12 \
        bucket_pairs/${SUBJECT}-${DONOR}_UM_QNAMEs.txt \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_MU_QNAMEs.txt \
        > bucket_pairs/donor_UM_MU_recipient.txt && \
    comm -12 \
        bucket_pairs/${SUBJECT}-${DONOR}_MU_QNAMEs.txt \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_MM_QNAMEs.txt \
        > bucket_pairs/donor_MU_MM_recipient.txt && \
    comm -12 \
        bucket_pairs/${SUBJECT}-${DONOR}_UM_QNAMEs.txt \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_MM_QNAMEs.txt \
        > bucket_pairs/donor_UM_MM_recipient.txt && \
    comm -12 \
        bucket_pairs/${SUBJECT}-${DONOR}_MM_QNAMEs.txt \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_UU_QNAMEs.txt \
        > bucket_pairs/donor_MM_UU_recipient.txt && \
    comm -12 \
        bucket_pairs/${SUBJECT}-${DONOR}_UU_QNAMEs.txt \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_MM_QNAMEs.txt \
        > bucket_pairs/donor_UU_MM_recipient.txt 
    "

# Count of reads per bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    cat bucket_pairs/donor_MU_UM_recipient.txt | wc -l && \
    cat bucket_pairs/donor_UM_MU_recipient.txt | wc -l && \
    cat bucket_pairs/donor_MU_MM_recipient.txt | wc -l && \
    cat bucket_pairs/donor_UM_MM_recipient.txt | wc -l && \
    cat bucket_pairs/donor_MM_UU_recipient.txt | wc -l && \
    cat bucket_pairs/donor_UU_MM_recipient.txt | wc -l 
    "

# Filtering the original BAM files using the inner join text files containing the common QNAMEs
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/picard-tools \
    /bin/bash -c \
    " \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${DONOR}_MU.bam \
        O=bucket_pairs/${SUBJECT}-${DONOR}_donor_MU_UM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_MU_UM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${RECIPIENT}_UM.bam \
        O=bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MU_UM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_MU_UM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${DONOR}_UM.bam \
        O=bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_MU_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_UM_MU_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${RECIPIENT}_MU.bam \
        O=bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_UM_MU_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_UM_MU_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${DONOR}_MU.bam \
        O=bucket_pairs/${SUBJECT}-${DONOR}_donor_MU_MM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_MU_MM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${RECIPIENT}_MM.bam \
        O=bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MU_MM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_MU_MM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${DONOR}_UM.bam \
        O=bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_MM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_UM_MM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${RECIPIENT}_MM.bam \
        O=bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_UM_MM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_UM_MM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${DONOR}_MM.bam \
        O=bucket_pairs/${SUBJECT}-${DONOR}_donor_MM_UU_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_MM_UU_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${RECIPIENT}_UU.bam \
        O=bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MM_UU_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_MM_UU_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${DONOR}_UU.bam \
        O=bucket_pairs/${SUBJECT}-${DONOR}_donor_UU_MM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_UU_MM_recipient.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=buckets/${SUBJECT}-${RECIPIENT}_MM.bam \
        O=bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_UU_MM_recipient.bam \
        READ_LIST_FILE=bucket_pairs/donor_UU_MM_recipient.txt \
        FILTER= includeReadList 
    "

# Count of reads per inner join filtered bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_MU_UM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MU_UM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_MU_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_UM_MU_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_MU_MM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MU_MM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_MM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MU_MM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_MM_UU_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MM_UU_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_UU_MM_recipient.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_UU_MM_recipient.bam | wc -l 
    "

# Merge _donor_MU_UM_recipient and _donor_UM_MU_recipient donor and recipient .bam files
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools merge -f \
        bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_donor_MU.bam \
        -@ ${THREADS} \
        bucket_pairs/${SUBJECT}-${DONOR}_donor_MU_UM_recipient.bam \
        bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_MU_recipient.bam && \
    samtools merge -f  \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU.bam \
        -@ ${THREADS} \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_MU_UM_recipient.bam \
        bucket_pairs/${SUBJECT}-${RECIPIENT}_donor_UM_MU_recipient.bam 
    "

# Count of reads per inner join filtered bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_donor_MU.bam | wc -l && \
    samtools view -@ ${THREADS} bucket_pairs/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU.bam | wc -l 
    "

# Verify mate-pair information between mates and fix if needed.
# This tool ensures that all mate-pair information is in sync between each read and its mate pair
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/picard-tools \
    /bin/bash -c \
    " \
    java -jar /usr/picard/picard.jar FixMateInformation \
        I=bucket_pairs/${SUBJECT}-${DONOR}_donor_UM_donor_MU.bam \
        O=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU.bam \
        SO=coordinate \
        ASSUME_SORTED=0 \
        VALIDATION_STRINGENCY=SILENT && \
    java -jar /usr/picard/picard.jar FixMateInformation \
        I=bucket_pairs/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU.bam \
        O=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU.bam \
        SO=coordinate \
        ASSUME_SORTED=0 \
        VALIDATION_STRINGENCY=SILENT
    "

# Count of reads per inner join filtered bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU.bam | wc -l && \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU.bam | wc -l 
    "

# Dedup _donor_MU_UM_recipient and _donor_UM_MU_recipient donor and recipient .bam files
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/picard-tools \
    /bin/bash -c \
    " \
    java -jar /usr/picard/picard.jar MarkDuplicates \
        I=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU.bam \
        O=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup.bam \
        METRICS_FILE=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 && \
    java -jar /usr/picard/picard.jar MarkDuplicates \
        I=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU.bam \
        O=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup.bam \
        METRICS_FILE=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
    "

# Count of reads per inner join filtered bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup.bam | wc -l && \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup.bam | wc -l 
    "

# Remove low complexity sequences
# First convert BAM files to FASTQ files
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/picard-tools \
    /bin/bash -c \
    " \
    java -jar /usr/picard/picard.jar SamToFastq \
        INPUT=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup.bam \
        FASTQ=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_1.fastq \
        SECOND_END_FASTQ=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_2.fastq \
        VALIDATION_STRINGENCY=SILENT  && \
    java -jar /usr/picard/picard.jar SamToFastq \
        INPUT=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup.bam \
        FASTQ=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_1.fastq \
        SECOND_END_FASTQ=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_2.fastq \
        VALIDATION_STRINGENCY=SILENT 
    "

# Filter FASTQ files to match given parameters
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/prinseq \
    /bin/bash -c \
    " \
    prinseq \
    --fastq=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_1.fastq \
    --out_good=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_1_good \
    --out_bad=null \
    -lc_method=dust \
    -lc_threshold=7 && \
    prinseq \
    --fastq=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_2.fastq \
    --out_good=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_2_good \
    --out_bad=null \
    -lc_method=dust \
    -lc_threshold=7 && \
    prinseq \
    --fastq=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_1.fastq \
    --out_good=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_1_good \
    --out_bad=null \
    -lc_method=dust \
    -lc_threshold=7 && \
    prinseq \
    --fastq=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_2.fastq \
    --out_good=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_2_good \
    --out_bad=null \
    -lc_method=dust \
    -lc_threshold=7 
    "

# Extract the SEQ_IDs from the FASTQ file into a .txt file
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    cat \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_1_good.fastq \
    | grep @ | sed 's/^@//' | sed 's/\/.$//' > \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_1_good.txt && \
    cat \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_2_good.fastq \
    | grep @ | sed 's/^@//' | sed 's/\/.$//' > \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_2_good.txt && \
    cat \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_1_good.fastq \
    | grep @ | sed 's/^@//' | sed 's/\/.$//' > \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_1_good.txt && \
    cat \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_2_good.fastq \
    | grep @ | sed 's/^@//' | sed 's/\/.$//' > \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_2_good.txt 
    "

# Concatenate the files for R1s and R2s
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    cat \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_1_good.txt \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_2_good.txt | sort --parallel=${THREADS} > \
    dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_good.txt && \
    cat \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_1_good.txt \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_2_good.txt | sort --parallel=${THREADS} > \
    dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_good.txt 
    "

# Filter the dedup BAM so that they only have our 'good' QNAMEs
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/picard-tools \
    /bin/bash -c \
    " \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup.bam \
        O=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex.bam \
        READ_LIST_FILE=dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_good.txt \
        FILTER= includeReadList  && \
    java -jar /usr/picard/picard.jar FilterSamReads \
        I=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup.bam \
        O=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam \
        READ_LIST_FILE=dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_good.txt \
        FILTER= includeReadList 
    "

# Count of reads per inner join filtered bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex.bam | wc -l && \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam | wc -l 
    "

# Filter so to include mapped reads only - first: index
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools sort -@ ${THREADS} \
        dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex.bam \
        -o dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex-sorted.bam && \
    samtools sort -@ ${THREADS} \
        dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam \
        -o dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex-sorted.bam && \
    samtools index -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex-sorted.bam && \
    samtools index -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex-sorted.bam 
    "

# Filter so to include mapped reads only - second: filter
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} -F 4 dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex.bam \
        -o dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_final.bam && \
    samtools view -@ ${THREADS} -F 4 dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam \
        -o dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_final.bam 
    "

# Count of reads per inner join filtered bucket grouping
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_final.bam | wc -l && \
    samtools view -@ ${THREADS} dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_final.bam | wc -l 
    "

# Extract list of reads from BAM file
# Output paired reads to separate files, discarding singletons, supplementary and secondary reads
sudo docker run -it --rm -v "${PIPELINE_FOLDER}/data/:/tmp/" dformoso/samtools \
    /bin/bash -c \
    " \
        samtools fasta \
        -1 dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex_1.fasta \
        -2 dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex_2.fasta \
        -0 /dev/null \
        -s /dev/null \
        -n dedup_and_low_complexity_filtering/${SUBJECT}-${DONOR}_donor_UM_donor_MU_dedup_lowcomplex.bam && \
        samtools fasta \
        -1 dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam_1.fasta \
        -2 dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam_2.fasta \
        -0 /dev/null \
        -s /dev/null \
        -n dedup_and_low_complexity_filtering/${SUBJECT}-${RECIPIENT}_recipient_UM_recipient_MU_dedup_lowcomplex.bam
    "

# Display BLAST databases on GCP
sudo docker run --rm ncbi/blast update_blastdb.pl --showall pretty --source gcp

# Download nt (nucleotide collection) database
sudo docker run --rm \
    -v ${PIPELINE_FOLDER}/blast/blastdb:/blast/blastdb:rw \
    -w /blast/blastdb \
    ncbi/blast \
    update_blastdb.pl --source gcp nt

# Download taxdb (taxonomy) database
sudo docker run --rm \
    -v ${PIPELINE_FOLDER}/blast/blastdb:/blast/blastdb:rw \
    -w /blast/blastdb \
    ncbi/blast \
    update_blastdb.pl --source gcp taxdb


# 2m34.414s
time sudo docker run --rm \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/blast/blastdb:/blast/blastdb:ro \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/blast/blastdb_custom:/blast/blastdb_custom:ro \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/blast/queries:/blast/queries:ro \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/blast/results:/blast/results:rw \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/data/:/tmp/ \
    ncbi/blast \
    blastn \
        -query /tmp/query1.fa \
        -db nt \
        -num_threads 24 \
        -evalue 1 \
        -outfmt '6 seqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus' \
        -out /blast/results/blastn.query_donor.out


"""
            qseqid = Query Seq-id
               qgi = Query GI
              qacc = Query accesion
           qaccver = Query accesion.version
              qlen = Query sequence length
            sseqid = Subject Seq-id
         sallseqid = All subject Seq-id(s), separated by a ';'
               sgi = Subject GI
            sallgi = All subject GIs
              sacc = Subject accession
           saccver = Subject accession.version
           sallacc = All subject accessions
              slen = Subject sequence length
            qstart = Start of alignment in query
              qend = End of alignment in query
            sstart = Start of alignment in subject
              send = End of alignment in subject
              qseq = Aligned part of query sequence
              sseq = Aligned part of subject sequence
            evalue = Expect value
          bitscore = Bit score
             score = Raw score
            length = Alignment length
            pident = Percentage of identical matches
            nident = Number of identical matches
          mismatch = Number of mismatches
          positive = Number of positive-scoring matches
           gapopen = Number of gap openings
              gaps = Total number of gaps
              ppos = Percentage of positive-scoring matches
            frames = Query and subject frames separated by a '/'
            qframe = Query frame
            sframe = Subject frame
              btop = Blast traceback operations (BTOP)
            staxid = Subject Taxonomy ID
          ssciname = Subject Scientific Name
          scomname = Subject Common Name
        sblastname = Subject Blast Name
         sskingdom = Subject Super Kingdom
           staxids = unique Subject Taxonomy ID(s), separated by a ';'
         sscinames = unique Subject Scientific Name(s), separated by a ';'
         scomnames = unique Subject Common Name(s), separated by a ';'
        sblastnames = unique Subject Blast Name(s), separated by a ';'
        sskingdoms = unique Subject Super Kingdom(s), separated by a ';'
            stitle = Subject Title
        salltitles = All Subject Title(s), separated by a '<>'
           sstrand = Subject Strand
             qcovs = Query Coverage Per Subject
           qcovhsp = Query Coverage Per HSP
            qcovus = Query Coverage Per Unique Subject (blastn only)
"""

# Create two separate BLAST files using BLASTalias. One for Donor, and one for Recipient

# Default filter min overlap is 20/50 (minimum length to filter out overlapping reads) - NOT IMPLEMENTED - Verbatum from LGTSeek
# Determining how to handle equal or better hits - NOT IMPLEMENTED - Verbatum from LGTSeek
# If our hit is equally as good as our best, then allow for more than one best hit - NOT IMPLEMENTED - Verbatum from LGTSeek
# Take a blast -m8 report and calculate the LCA - NOT IMPLEMENTED - Verbatum from LGTSeek

# Determine LGT from best Hit - NOT IMPLEMENTED - Verbatum from LGTSeek
# Tab-delimited output file with the following fields:
# 1) Read name
# 2) Best Euk hit bitscore
# 3) Best Bacteria hit bitscore
# 4) Euk bit / Bac bit
# 5) Euk bit - Bac bit (also known as h-score)

# Find (and filter) all the Donor (Pseudo) reads in the Eukaryots BLAST file
# Find (and filter) all the recipient (Human) reads in the Bacterial BLAST file
# Can have multiple m8 entries for 'best_hit'
# Also save best BLAST 

# Create summary info file containing each Donor hit
# Create summary info file containing each Recipient hit

# Filter BAM file by read ids
# Run mpileup






###############################

# BWA Align Parameters
# $;MISMATCH_PENALTY$; = 3
# $;MAX_GAP_OPENS$; = 1
# $;MAX_GAP_EXTNS$; = -1
# $;GAP_OPEN_PENALTY$; = 11
# $;GAP_EXTN_PENALTY$; = 4
# $;MAX_THREADS$; = 1
# $;MAX_ALIGN$; = 1
# $;SOFTCLIP_MIN$; = 24
# $;KEEP_MM$; = 1


# BLAST Parameters
# ;-evalue  Expectation value (E) default = 10.0
# $;EXPECT$; = 1e-5
# ;-max_target_seqs Maximum number of aligned sequences to keep default=150
# $;MAX_TARGET_SEQS$; = 150
# ;; Example options to put here:  -gilist/-negative_gilist, -task
# $;OTHER_OPTS$; =
# $;COMPRESS_RAW_OUTPUT$; = 0
# $;COMPRESS_BSML_OUTPUT$; = 0
# ;; Filter the HSPs for use in the %identity/similarity/coverage for each seq-pair-alignment in the output bsml. Setting this to 0 could show very low # scores on the seq-pair-alignment despite having a very high scoring HSP.
# $;FILTER_HSPS_FOR_STATS$;=1



sudo docker run -it --rm \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/blast/queries:/blast/queries:ro \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/blast/results:/blast/results:rw \
    -v /media/formoso/NVMe/research/genomics/pipelines/lgtsearch/data/:/tmp/ \
    ncbi/blast:2.10.1

update_blastdb.pl --source gcp nt
update_blastdb.pl --source gcp taxdb

time blastn \
        -query /tmp/query1.fa \
        -db nt \
        -num_threads 24 \
        -evalue 1 \
        -outfmt '6 seqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus' \
        -out /blast/results/blastn.query_donor.out

