version 1.0

import "structs/compute.wdl"

task trim {
    input {
        File fastq_1
        File fastq_2
        String out_fastq_1_paired = "~{basename(fastq_1, ".fastq")}_paired_trimmed.fastq"
        String out_fastq_1_unpaired = "~{basename(fastq_1, ".fastq")}_unpaired_trimmed.fastq"
        String out_fastq_2_paired = "~{basename(fastq_2, ".fastq")}_paired_trimmed.fastq"
        String out_fastq_2_unpaired = "~{basename(fastq_2, ".fastq")}_unpaired_trimmed.fastq"
        String out_trimlog = "out_trimlog.txt"
        String? adapter = "TruSeq2-PE.fa"
        Int? seed_mismatches = 2
        Int? paired_clip_threshold = 30
        Int? unpaired_clip_threshold = 10
        Int? leading = 3
        Int? trailing = 3
        Int? sliding_window_quality = 15
        Int? sliding_window_length = 4
        Int? min_length = 36
        Boolean? is_phred33 = true
        Resources resources
    }

    command <<<
        set -e
        cat << EOF > all_adapters.fa
        >Illumina_Universal_Adapter
        AGATCGGAAGAG
        >Illumina Small RNA 3 Adapter
        TGGAATTCTCGG
        >Illumina Small RNA 5 Adapter
        GATCGTCGGACT
        >Nextera Transposase Sequence
        CTGTCTCTTATA
        >SOLID Small RNA Adapter
        CGCCTTGGCCGT
        >Nextera - PrefixNX/1
        AGATGTGTATAAGAGACAG
        >PrefixNX/2
        AGATGTGTATAAGAGACAG
        >Trans1
        TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
        >Trans1_rc
        CTGTCTCTTATACACATCTGACGCTGCCGACGA
        >Trans2
        GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
        >Trans2_rc
        CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
        >TruSeq2-PE.fa - PrefixPE/1
        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >PrefixPE/2
        CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
        >PCR_Primer1
        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >PCR_Primer1_rc
        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >PCR_Primer2
        CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
        >PCR_Primer2_rc
        AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
        >FlowCell1
        TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
        >FlowCell2
        TTTTTTTTTTCAAGCAGAAGACGGCATACGA
        >TruSeq2-SE.fa - TruSeq2_SE
        AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
        >TruSeq2_PE_f
        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        >TruSeq2_PE_r
        AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
        >TruSeq3-PE-2.fa - PrefixPE/1
        TACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >PrefixPE/2
        GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
        >PE1
        TACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >PE1_rc
        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
        >PE2
        GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
        >PE2_rc
        AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        >TruSeq3-PE.fa - PrefixPE/1
        TACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >PrefixPE/2
        GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
        >TruSeq3-SE.fa - TruSeq3_IndexedAdapter
        AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        >TruSeq3_UniversalAdapter
        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
        >Reverse_adapter
        AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Universal_Adapter
        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
        >pcr_dimer
        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
        >PCR_Primers
        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT
        >TruSeq_Adapter_Index_1_6
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_2
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_3
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_4
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_5
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_6
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_7
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_8
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_9
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_10
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_11
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_12
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_13
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_14
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_15
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_16
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_18_7
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_19
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_20
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_21
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_22
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_23
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_25
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG
        >TruSeq_Adapter_Index_27
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG
        >I5_Nextera_Transposase_1
        CTGTCTCTTATACACATCTGACGCTGCCGACGA
        >I7_Nextera_Transposase_1
        CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
        >I5_Nextera_Transposase_2
        CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC
        >I7_Nextera_Transposase_2
        CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]501
        GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502
        GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503
        GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]504
        GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505
        GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506
        GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507
        GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508
        GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517
        GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701
        CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702
        CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703
        CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704
        CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705
        CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706
        CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707
        CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N708
        CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N709
        CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710
        CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711
        CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712
        CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG
        >I5_Primer_Nextera_XT_Index_Kit_v2_S502
        GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S503
        GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S505
        GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S506
        GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S507
        GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S508
        GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S510
        GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S511
        GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S513
        GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S515
        GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S516
        GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S517
        GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S518
        GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S520
        GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S521
        GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I5_Primer_Nextera_XT_Index_Kit_v2_S522
        GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I7_Primer_Nextera_XT_Index_Kit_v2_N701
        CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N702
        CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N703
        CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N704
        CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N705
        CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N706
        CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N707
        CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N710
        CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N711
        CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N712
        CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N714
        CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N715
        CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N716
        CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N718
        CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N719
        CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N720
        CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N721
        CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N722
        CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N723
        CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N724
        CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N726
        CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N727
        CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N728
        CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG
        >I7_Primer_Nextera_XT_Index_Kit_v2_N729
        CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG
        >I5_Adapter_Nextera
        CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT
        >I7_Adapter_Nextera_No_Barcode
        CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG
        >Nextera_LMP_Read1_External_Adapter
        GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        >Nextera_LMP_Read2_External_Adapter
        GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        >RNA_Adapter_(RA5)_part_#_15013205
        GATCGTCGGACTGTAGAACTCTGAAC
        >RNA_Adapter_(RA3)_part_#_15013207
        CCTTGGCACCCGAGAATTCCA
        >Stop_Oligo_(STP)_8
        CCACGGGAACGTGGTGGAATTC
        >RNA_RT_Primer_(RTP)_part_#_15013981
        TGGAATTCTCGGGTGCCAAGGC
        >RNA_PCR_Primer_(RP1)_part_#_15013198
        TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT
        >RNA_PCR_Primer_Index_1_(RPI1)_2,9
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_2_(RPI2)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_3_(RPI3)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_4_(RPI4)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_5_(RPI5)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_6_(RPI6)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_7_(RPI7)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_8_(RPI8)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_9_(RPI9)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_10_(RPI10)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_11_(RPI11)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_12_(RPI12)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_13_(RPI13)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_14_(RPI14)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_15_(RPI15)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_16_(RPI16)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_17_(RPI17)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_18_(RPI18)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_19_(RPI19)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_20_(RPI20)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_21_(RPI21)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_22_(RPI22)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_23_(RPI23)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_24_(RPI24)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_25_(RPI25)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_26_(RPI26)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_27_(RPI27)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_28_(RPI28)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_29_(RPI29)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_30_(RPI30)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_31_(RPI31)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_32_(RPI32)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_33_(RPI33)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_34_(RPI34)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_35_(RPI35)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_36_(RPI36)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_37_(RPI37)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_38_(RPI38)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_39_(RPI39)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_40_(RPI40)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_41_(RPI41)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_42_(RPI42)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_43_(RPI43)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_44_(RPI44)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_45_(RPI45)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_46_(RPI46)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_47_(RPI47)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG
        >RNA_PCR_Primer_Index_48_(RPI48)
        TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG
        >PhiX_read1_adapter
        AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA
        >PhiX_read2_adapter
        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA
        >BLOuterF
        AAATCTCTAGCAGTGGCGCCCGAACAG
        >BlouterR
        TGAGGGATCTCTAGTTACCAGAGTC
        >275F
        ACAGGGACCTGAAAGCGAAAG
        >280R
        CTAGTTACCAGAGTCACACAACAGACG
        EOF
        trimmomatic \
            PE \
            ~{true="-phred33" false="-phred64" is_phred33} \
            ~{"-threads " + resources.cpu} \
            ~{"-trimlog " + out_trimlog} \
            ~{fastq_1} \
            ~{fastq_2} \
            ~{out_fastq_1_paired} \
            ~{out_fastq_1_unpaired} \
            ~{out_fastq_2_paired} \
            ~{out_fastq_2_unpaired} \
            TOPHRED33 \
            ILLUMINACLIP:~{adapter}:~{seed_mismatches}:~{paired_clip_threshold}:~{unpaired_clip_threshold} \
            ~{"LEADING:" + leading} \
            ~{"TRAILING:" + trailing} \
            SLIDINGWINDOW:~{sliding_window_length}:~{sliding_window_quality} \
            ~{"MINLEN:" + min_length}
    >>>

    output {
        File fastq_1_paired = out_fastq_1_paired 
        File fastq_1_unpaired = out_fastq_1_unpaired 
        File fastq_2_paired = out_fastq_2_paired 
        File fastq_2_unpaired = out_fastq_2_unpaired 
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/trimmomatic:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
