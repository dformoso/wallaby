import sys
import pandas
import subprocess

insertion_table_name = sys.argv[1]
srr = sys.argv[2]
donor_name = sys.argv[3]
recipient_name = sys.argv[4]

insertion_table = pandas.read_csv(insertion_table_name)

crossings = ["1Md_2Mr", "1Md_1Mr", "2Md_1Mr"]
padding_left = 0
padding_right = 0

recipient_1Md_2Mr = srr + "-to-" + recipient_name + "_" + "1Md_2Mr" + "_filtered.bam"
recipient_1Md_1Mr = srr + "-to-" + recipient_name + "_" + "1Md_1Mr" + "_filtered.bam"
recipient_2Md_1Mr = srr + "-to-" + recipient_name + "_" + "2Md_1Mr" + "_filtered.bam"

donor_1Md_2Mr = srr + "-to-" + donor_name + "_" + "1Md_2Mr" + "_filtered.bam"
donor_1Md_1Mr = srr + "-to-" + donor_name + "_" + "1Md_1Mr" + "_filtered.bam"
donor_2Md_1Mr = srr + "-to-" + donor_name + "_" + "2Md_1Mr" + "_filtered.bam"

for srr_id in insertion_table["id"]:
    insertion_table_srr = insertion_table[insertion_table["srr"] == srr]
    insertion_table_srr_id = insertion_table_srr[insertion_table_srr["id"] == srr_id]
    chromosome = str(insertion_table_srr_id['chr'].values[0])
    start = str(insertion_table_srr_id['start'].values[0] - padding_left)
    stop = str(insertion_table_srr_id['stop'].values[0] + padding_right)
    locus = str(chromosome + ":" + start + "-" + stop)
    
    subprocess.getoutput("/tmp/samtools-1.14/samtools view " + recipient_1Md_2Mr + " " + locus + " | cut -f1 | sort | uniq > reads1.txt")
    subprocess.getoutput("/tmp/samtools-1.14/samtools view " + recipient_1Md_1Mr + " " + locus + " | cut -f1 | sort | uniq > reads2.txt")
    subprocess.getoutput("/tmp/samtools-1.14/samtools view " + recipient_2Md_1Mr + " " + locus + " | cut -f1 | sort | uniq > reads3.txt")
    subprocess.getoutput("cat reads1.txt reads2.txt reads3.txt > reads.txt")
    subprocess.getoutput("rm reads1.txt reads2.txt reads3.txt")

    subprocess.getoutput("/tmp/samtools-1.14/samtools view -N reads.txt -o " + srr + "-to-" + donor_name + "_1Md_2Mr_filtered_id" + str(srr_id) + ".bam " + donor_1Md_2Mr)
    subprocess.getoutput("/tmp/samtools-1.14/samtools view -N reads.txt -o " + srr + "-to-" + donor_name + "_1Md_1Mr_filtered_id" + str(srr_id) + ".bam " + donor_1Md_1Mr)
    subprocess.getoutput("/tmp/samtools-1.14/samtools view -N reads.txt -o " + srr + "-to-" + donor_name + "_2Md_1Mr_filtered_id" + str(srr_id) + ".bam " + donor_2Md_1Mr)
    subprocess.getoutput("/tmp/samtools-1.14/samtools view -N reads.txt -o " + srr + "-to-" + recipient_name + "_1Md_2Mr_filtered_id" + str(srr_id) + ".bam " + recipient_1Md_2Mr)
    subprocess.getoutput("/tmp/samtools-1.14/samtools view -N reads.txt -o " + srr + "-to-" + recipient_name + "_1Md_1Mr_filtered_id" + str(srr_id) + ".bam " + recipient_1Md_1Mr)
    subprocess.getoutput("/tmp/samtools-1.14/samtools view -N reads.txt -o " + srr + "-to-" + recipient_name + "_2Md_1Mr_filtered_id" + str(srr_id) + ".bam " + recipient_2Md_1Mr)