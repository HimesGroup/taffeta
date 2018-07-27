#! /usr/bin/python

###
# Reference files
###

# hg38 reference files
hg38_fa = "/project/bhimeslab/Reference/hg38/genome.ERCC.fa"
hg38_gtf = "/project/bhimeslab/Reference/hg38/genes.gtf"
hg38_ref = "/project/bhimeslab/Reference/hg38/refFlat.txt"
hg38_ERCC_gtf = "/project/bhimeslab/Reference/hg38/genes.ERCC.gtf"
hg38_star_index_dir = "/project/bhimeslab/Reference/hg38/STAR_index"
hg38_rRNA_gtf = "/project/bhimeslab/Reference/hg38/rRNA_hg38.gtf"

# hg19 reference files
hg19_fa = "/project/bhimeslab/Reference/hg19/genome.ERCC.fa"
hg19_gtf = "/project/bhimeslab/Reference/hg19/genes.gtf"
hg19_ref = "/project/bhimeslab/Reference/hg19/refFlat.txt"
hg19_ERCC_gtf = "/project/bhimeslab/Reference/hg19/genes.ERCC.gtf"
hg19_star_index_dir = "/project/bhimeslab/Reference/hg19/STAR_index"

# mm38 reference files
mm38_fa = "/project/bhimeslab/Reference/mm38/mm.GRCm38.genome.ERCC92.fa"
mm38_gtf = "/project/bhimeslab/Reference/mm38/mm.GRCm38.75.genes.gtf"
mm38_ref = "/project/bhimeslab/Reference/mm38/refFlat.txt"
mm38_ERCC_gtf = "/project/bhimeslab/Reference/mm38/mm.GRCm38.75.genes.ERCC.gtf"
mm38_star_index_dir = "/project/bhimeslab/Reference/mm38/STAR_index"

# mm10 reference files
mm10_fa = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genome.ERCC92.fa"
mm10_gtf = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genes.gtf"
mm10_ref = "/project/bhimeslab/Reference/mm38/refFlat_UCSC.txt"
mm10_ERCC_gtf = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genes.ERCC.gtf"
mm10_star_index_dir = "/project/bhimeslab/Reference/mm38/STAR_index_mm10"

# rn6 reference files
rn6_fa = "/project/bhimeslab/Reference/rn6/genome.ERCC.fa"
rn6_gtf = "/project/bhimeslab/Reference/rn6/genes.gtf"
rn6_ref = "/project/bhimeslab/Reference/rn6/refFlat.txt"
rn6_ERCC_gtf = "/project/bhimeslab/Reference/rn6/genes.ERCC.gtf"
rn6_star_index_dir = "/project/bhimeslab/Reference/rn6/STAR_index"

# susScr3 reference files
susScr3_fa = "/project/bhimeslab/Reference/susScr3/Ensembl_susScr3_genome.ERCC92.fa"
susScr3_gtf = "/project/bhimeslab/Reference/susScr3/genes.gtf"
susScr3_ref = "/project/bhimeslab/Reference/susScr3/refFlat.txt"
susScr3_ERCC_gtf = "/project/bhimeslab/Reference/susScr3/Ensembl_susScr3_genes.ERCC.gtf"
susScr3_star_index_dir = "/project/bhimeslab/Reference/susScr3/STAR_index"

# ERCC
ERCC_only = "/project/bhimeslab/Reference/ERCC92.gtf"


###
# Informatics Tools
###

# path and directory
trimmomatic="/opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar"
picard_dir="/opt/software/picard/picard-tools-1.96/"

# version
trimmomatic_version="0.32"
fastqc_version="0.11.7"
star_version="2.5.2b"
samtools_version="1.8"
bamtools_version="2.3.0"
picard_version="1.96"

###
# Other information
###

# author name and contact info (shown in .html report)
author="Mengyuan Kan (mengykan@upenn.edu)"
