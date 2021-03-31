#!/bin/env python

"""
variant-calling from bulk RNA-Seq pipeline.
Authors:  Fahimeh Palizban (palizbanfahimeh@gmail.com)
Description:
This program implements a workflow pipeline for next generation
sequencing variant detection using the integration of mapping to both reference genome and transcriptome.
It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

"""

import sys
import argparse
import os
import os.path
import csv
from operator import itemgetter

prog = "RNA_Var.py"

version = """%prog
RNA_var is free for non-commercial use without warranty.
============================================================================
"""

usage = """Usage: %prog[-h] -1 Read1.fastq -2 Read2.fastq -o outputDir -index_dir reference_genome -transcriptome
                  reference_transcriptome [-t # threads, default: 4] [-g Using gzip input files, default: False] [-c Minimum coverage, default: 12] -picard picard.jar -gatk gatk.jar"""


def main(): 
    print("Staaaaart")
    parser = argparse.ArgumentParser(description='RNA_var: An All-In-One Workflow for variant Detection, Quantification, and Visualization')

    parser.add_argument('-t', '--n_thread', required=False, default='4',
                        type=str, help='Number of threads')
  
    parser.add_argument('-picard',  required=True,
                        type=str, help='Path to picard.jar file')
  
    parser.add_argument('-gatk',  required=True,
                        type=str, help='Path to gatk.jar file')

    parser.add_argument('-1', '--fq1', required=True, metavar='read1.fastq', type=str,
                        help='Path to Read 1 of the paired-end RNA-seq')

    parser.add_argument('-2', '--fq2', required=True, metavar='read2.fastq', type=str,
                        help='Path to Read 2 of the paired-end RNA-seq')

    parser.add_argument('-o', '--out', required=True, metavar='outputDir', type=str,
                        help='Path to the output directory')

    parser.add_argument('-index_genome', '--index_dir', required=True, metavar='reference_genome', type=str,
                        help='Path to directory containing index files for human genome')

    parser.add_argument('-transcriptome', '--transcriptome', required=True, metavar=' reference_transcriptome', type=str,
                        help='Path to directory containing index files and fasta for reference transcriptome')

    parser.add_argument('-c', '--min_coverage', required=False, metavar='min_coverage', default='12', type=str,
                        help='Minimum read coverage')

    parser.add_argument('-g', '--gzip', required=False,
                        metavar='True/False', default='False', type=str,
                        help="""Toggle "-g True" for using gunzipped FASTQ input""")

    parser.add_argument('-s', '--consensus', required=False,
                        metavar='True/False', default='False', type=str,
                        help="""Toggle "-s True" to generate consensus fasta file""")

    args = parser.parse_args()

    fq1 = os.path.abspath(args.fq1)
    fq2 = os.path.abspath(args.fq2)

    try:
        os.path.isfile(fq1)
    except IOError:
        print('Error: Unable to locate Read 1 FASTQ file. Please check your path and try again.')
        sys.exit()

    try:
        os.path.isfile(fq2)
    except IOError:
        print('Error: Unable to locate Read 2 FASTQ file. Please check your path and try again.')
        sys.exit()

    index_dir = os.path.abspath(args.index_dir)
    transcriptome = os.path.abspath(args.transcriptome)
    picard = os.path.abspath(args.picard)
    gatk =  os.path.abspath(args.gatk)

    try:
        os.path.isfile(index_dir)
    except IOError:
        print('Error: Unable to locate human genome reference files. Please check your path and try again.')
        sys.exit()

    try:
        os.path.isfile(transcriptome)
    except IOError:
        print('Error: Unable to locate viral genome reference files. Please check your path and try again.')
        sys.exit()

    out = os.path.abspath(args.out)
    n_thread = args.n_thread
    gzip = args.gzip
    min_coverage = args.min_coverage
    consensus = args.consensus

    ##os.system('ulimit -n 2048')

    print("Aligning to human reference genome using STAR")

    def alignment():
        #cmd1='hisat2 -x '+index_dir+' -1 '+fq1+' -2'+fq2+' -S '+out+'/accepted_hits.sam -p '+n_thread+' --quiet'
        if gzip == "False":
            cmd1 = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/accepted_hits. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate'
        elif gzip =="True":
            cmd1 = 'STAR --runThreadN ' + n_thread + ' --genomeDir ' + index_dir + ' --readFilesIn ' + fq1 + ' ' + fq2 + ' --outFileNamePrefix ' + out + '/accepted_hits. --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'
        else:
            print("Error: Invalid input for -g/--gzip. Input should be 'True' if using gunzipped files")
            sys.exit()
        print('Running ', cmd1)
        ##os.system(cmd1)

    alignment()

    print("Aligning to reference transcriptome using BWA")

    def transcriptome_alignment():
        cmd2 = 'bwa mem ' + fq1 + ' ' + fq2 + ' ' + transcriptome + '/*.fasta' + ' -t ' + n_thread + ' > ' + out +'trans.sam'
        print('Running ', cmd2)
        ##os.system(cmd2)

    transcriptome_alignment()

    print("Converting SAM to BAM")

    def sam_to_bam():
        cmd3 = 'samtools view -Sb -h ' + out + '/trans.sam > ' + out + '_trans.bam'
        print('Running ', cmd3)
        ##os.system(cmd3)

    sam_to_bam()

    print("Sorting BAM")

    def sort():
        cmd4 = 'samtools sort -@ ' + n_thread + ' ' + out + '/trans.bam -o ' + out + '_trans_Coord_sorted.bam'
        print('Running ', cmd4)
        ##os.system(cmd4)

        #cmd5 = 'samtools sort -n -@ ' + n_thread + ' ' + out + '/unmapped_aln.bam -o' + out + '/unmapped_aln_sorted.bam'
        #print('Running ', cmd5)
        ###os.system(cmd5)

    sort()

    def addreadgroup():
        cmd5 = 'java -jar' + picard +  ' AddOrReplaceReadGroups ' + 'I=' + out + '/accepted_hits.bam ' + 'O=' + out + '_RG.bam' + ' RGID=4 RGLB=twist RGPL=illumina RGPU=unit1 RGSM=' + out
        print('Running ', cmd5)
        ##os.system(cmd5)        
        cmd6 = 'java -jar' + picard +  ' AddOrReplaceReadGroups ' + 'I='+ out + 'trans_Coord_sorted.bam' + 'O=' + out + '_Trans_RG.bam' + 'RGID=4 RGLB=twist RGPL=illumina RGPU=unit1 RGSM=' + out
        print('Running ', cmd6)
	##os.system(cmd6)


    addreadgroup()

    def markduplicate():
        cmd7 = 'java -jar' + picard + ' ' + 'MarkDuplicates ' +  'I=' + out + '_RG.bam' +   'O=' + out + 'MarkDup_RG.Bam' + 'M=' + out + 'marked_dup_metrics.tx'
        cmd8 = 'java -jar' + picard + ' ' + 'MarkDuplicates ' +  'I=' + out + '_Trans_RG.bam' +   'O='  + out + 'MarkDup_Trans_RG.bam '+ 'M=' + out +  'trans_marked_dup_metrics.tx'

    print("Indexing BAM")

    def index():
        cmd6 = 'samtools index ' + out + '/unmapped_aln_Coord_sorted.bam'
        print('Running ', cmd6)
        ##os.system(cmd6)

    index()

    print("Processing Viral Counts using StringTie")

    def stringtie():
        cmd7 = 'stringtie -A -l ' + out + '/accepted_hits.bam -G ' + index_dir+ '/*.gtf -o ' + out + '/stringtie/stringtie.gtf -p ' + n_thread
        print('Running ', cmd7)
        ##os.system(cmd7)

    stringtie()

    def Indexmarkedbam():
        cmd8 = 'java -jar' + picard + ' ' + 'BuildBamIndex' + 'I=' + out + 'MarkDup_RG.Bam' 
        print('Running ', cmd8)
        ##os.system(cmd8)
        cmd9 = 'java -jar' + picard + ' ' + 'BuildBamIndex' + 'I=' + out + 'MarkDup_Trans_RG.bam '
        print('Running ', cmd9)
        ##os.system(cmd9)

    Indexmarkedbam()
   
    def splitNcigar():
       cmd10 = 'java -jar' + gatk + '' + 'SplitNCigarReads -R ' + index_dir+ '/*.fasta -I ' + out + 'MarkDup_RG.Bam -O ' + out + '_splitNcigar.Bam '


    def variantCalling():
        cmd10  = 'samtools mpileup -f ' + transcriptome + '/*.fasta  -g ' + out + 'MarkDup_Trans_RG.bam ' + '|bcftools call -m > ' +out + '_trans.vcf' 
         print('Running ', cmd11)
        ##os.system(cmd11)   
        cmd12 = 'java -jar ' + gatk + ' HaplotypeCaller -R ' + index_dir + '/*.fasta -I ' + out + '_splitNcigar.Bam -O ' + '_gen.vcf
        print('Running ', cmd12)
        ##os.system(cmd12) 
    
    def vcfmerg():
        cmd13 = 'bcftools merg ' + out + '/*.vcf -Oz -o ' + out + '_Merged.vcf' 
        print('Running ', cmd13)
        ##os.system(cmd13) 

   def ROH():
       cmd14 = 'vcftools --vcf ' + '' + out + '_Merged.vcf --plink --out' + out + 'plink'
       print('Running ', cmd14)
        ##os.system(cmd14) 
       cmd15 = 'plink --file ' + out + 'plink'+ ' –homozyg '

#if __name__ == '__main__':
main()
