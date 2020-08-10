version development

import "../tasks/samtools.wdl" as samtools
import "../tasks/tools.wdl" as tools
import "../tasks/picard.wdl" as picard
import "../tasks/prinseq.wdl" as prinseq
import "../tasks/structs/compute.wdl"

# Removes low complexity sequences
workflow main {

    input {
        Array[File] bam_files
        Resources resources
        String lc_method
        String lc_threshold
    }

    scatter (bam in bam_files) {

        # First convert BAM file to FASTQ files
        call picard.sam_to_fastq as sam_to_fastq { 
            input: 
                file = bam, 
                out_fastq_1 = "~{basename(bam)}.1.fastq", 
                out_fastq_2 = "~{basename(bam)}.2.fastq",
                resources = resources
        }

        # Filter FASTQ files to match given parameters
        call prinseq.matching as prinseq_fastq1 { 
            input: 
                file = sam_to_fastq.fastq_1, 
                out_file = "~{basename(sam_to_fastq.fastq_1)}",
                lc_method = lc_method,
                lc_threshold = lc_threshold,
                resources = resources
        }
        
        call prinseq.matching as prinseq_fastq2 { 
            input: 
                file = sam_to_fastq.fastq_2, 
                out_file = "~{basename(sam_to_fastq.fastq_2)}",
                lc_method = lc_method,
                lc_threshold = lc_threshold,
                resources = resources
        }
        
        # Extract the SEQ_IDs from the FASTQ files into .txt files
        call tools.seq_ids_from_fastq as seqids_fastq1 { 
            input: 
                file = prinseq_fastq1.out, 
                out_file = "~{basename(prinseq_fastq1.out)}.txt",
                resources = resources
        }
       
        call tools.seq_ids_from_fastq as seqids_fastq2 { 
            input: 
                file = prinseq_fastq2.out, 
                out_file = "~{basename(prinseq_fastq2.out)}.txt",
                resources = resources
        }
        
        # Concatanate the text files for fastq1 and fastq2 into one
        call tools.concat_text as concat_fastqs { 
            input: 
                file_1 = seqids_fastq1.out, 
                file_2 = seqids_fastq2.out, 
                out_file = "~{basename(bam)}.merged.txt",
                resources = resources
        }
        
        # Filter the BAM file so that they only have our selected SEQ_IDs
        call picard.filter_valid_reads as complex { 
            input: 
                file = bam, 
                txt = concat_fastqs.out, 
                out_file = "~{basename(bam)}_complex.bam",
                resources = resources
        }
        
        # Extract list of reads from BAM file
        # Output paired reads to separate files, 
        # discarding singletons, supplementary and secondary reads
        call samtools.bam_to_fastas as bam_to_fastas {
            input :
                file = complex.out,
                out_file_1 = "~{basename(bam)}.1.fasta",
                out_file_2 = "~{basename(bam)}.2.fasta",
                resources = resources
        }

    }

    output {
        Array[File] bams = complex.out
        Array[File] fastas = flatten([bam_to_fastas.out_1, bam_to_fastas.out_2])
    }

}