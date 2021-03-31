# RNA_SeqVariantCalling-pipeline
A novel pipeline that increases the yield of variant calling from RNA-Seq by concurrent use of genome and transcriptome references in parallel. By using the integrated workflow, variant calling and extracting differentially expressed genes will be identified at the same time.
In the following part all the insilico analysis related to .... is provided.

_It is important to note that in the current project we did not devlope a software package and all the scripts have been implemented to facilitate the automated use of the pipeline._

## Dependencies
Before running the pipeline, there are several tools and packages that need to be obtained and installed:

- Python
- FastQC
- Trimmomatic
- STAR
- BWA
- Picard
- GATK
- SAMtools
- VCFtools
- PLINK
- StringTie



## Installation
All the instructions are based on running the pipeline on Ubuntu.

###### FastQC
$ sudo apt install fastqc
###### Trimmomatic
The binary jar file can be obtained from : http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
###### STAR
$ git clone https://github.com/alexdobin/STAR.git
Compile
cd STAR/source
make STAR
