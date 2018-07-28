taffeta
=======

Reproducible analysis and validation of RNA-Seq data

Authors: Maya Shumyatcher, Mengyuan Kan, Blanca Himes.

## Overview

The goal of taffeta is to preform reproducible analysis and validation of RNA-Seq data, as a part of [RAVED pipeline](https://github.com/HimesGroup/raved):
  * Download SRA .fastq data
  * Perform preliminary QC
  * Align reads to a reference genome
  * Perform QC on aligned files
  * Create a report that can be used to verify that sequencing was successful and/or identify sample outliers
  * Perform differential expression of reads aligned to transcripts according to a given reference genome
  * Create a report that summarizes the differential expression results

Generate LSF scripts with download commands that can be processed in parallel. Specify --fastqc option to run FastQC for .fastq files downloaded.

## Dependencies

* For alignment and quantification: STAR, HTSeq, kallisto (older versions used bowtie2, tophat, cufflinks, cummerbund, whose options are no longer available).
* For DE analysis: R packages DESeq2 and sleuth
* For QC: fastqc, trimmomatic, samtools, bamtools, picardtools.
* Annotation files should be available: reference genome fasta, gtf, refFlat, and index files. We use ERCC spike-ins, so our reference files include ERCC mix transcripts. 
* For adapter trimming, we provide Ilumina TruSeq single index and Illumina unique dual (UD) index adpter and primer sequences. Users can tailor this file by adding sequences from other protocols.
* The Python scripts make use of modules that include subprocess, os, argparse, sys.
* To create reports, R and various libraries should be available, such as DT, gplots, ggplot2, reshape2, rmarkdown, RColorBrewer, plyr, dplyr, lattice, ginefilter, biomaRt. Note that pandoc version 1.12.3 or higher is required to generate HTML report from current RMD files.

## Workflow

### Download data from GEO/SRA

Run script rnaseq_sra_download.py to download .fastq files from SRA. Users can provide a phenotype file with a SRA_ID column. Download samples corresponding to the SRA_ID. If a phenotype file is not provided, use phenotype information from GEO. SRA_ID is retrieved from the field relation.1. Corresponding ftp addresses are obtained from SRA sqilte database.

> rnaseq_sra_download.py --geo_id <i>GEO_Accession</i> --path_start <i>output_path</i> --project_name <i>output_prefix</i> --template_dir <i>templete_file_directory</i> --fastqc

Output files: 1. <i>project_name</i>\_SRAdownload_RnaSeqReport.html; 2. <i>project_name</i>\_sraFile.info 3. <i>GEO_Accession</i>\_withoutQC.txt 4. FastQC results


### User-tailored phentoype file

The sample info file used in the following steps should be provided by users.

Required columns: 'Sample' column containing sample ID, 'Status' column containing variables of comparison state. 'R1' and/or 'R2' columns containing full paths of .fastq files.

Other columns: 'Treatment', 'Disease', 'Donor' (donor or cell line ID if in vitro treatment is used), 'Tissue', 'ERCC\_Mix' if ERCC spike-in sample is used, 'protocol' designating sample preparation kit information.

'Index' column contains index sequence for each sample. If provided, trim raw .fastq files based on corresponding adapter sequences

Most GEO phenotype data do not have index information. However, FastQC is able to detect them as "Overrepresented sequences". Users can tailor Index column based on FastQC results. We provide a file with most updated adapter and primer sequences for FastQC detection.

### Alignment, quantification and QC

Run rnaseq_align_and_qc.py to perform: 1) adapter trimming, 2) fastqc, 3) alignment, 4) quantification, 5) QC metrics 

> rnaseq_align_and_qc.py --project_name <i>output_prefix</i> --sample_in <i>sample_info_file.txt</i> --aligner star --path_start <i>output_path</i> --template_dir <i>templete_file_directory</i> --fastqc

The "--aligner star" option indicates that star should be used as the aligner (default is tophat). The "--discovery no" option refers to using --no-novel-juncs and --transcriptome-only options with tophat; while this option is not relevant if using star as the aligner, it must be specified in order for the script to run. Following execution of this script, various output files will be written for each sample in directories structured as:

<i>path_start</i>/<i>sample_name</i>/<i>sample_name</i>_R1_Trimmed.fastqc <br>
<i>path_start</i>/<i>sample_name</i>/<i>sample_name</i>_R2_Trimmed.fastqc <br>
<i>path_start</i>/<i>sample_name</i>/<i>sample_name</i>_R1_fastqc <br>
<i>path_start</i>/<i>sample_name</i>/<i>sample_name</i>_R2_fastqc <br>
<i>path_start</i>/<i>sample_name</i>/<i>sample_name</i>_ReadCount <br>
<i>path_start</i>/<i>sample_name</i>/<i>aligner</i>_out <br>
<i>path_start</i>/<i>sample_name</i>/<i>quantification_tool</i>_out <br>
 ...

Run rnaseq_align_and_qc_report.py. Create an HTML report of QC and alignment summary statistics for RNA-seq samples.

> python rnaseq_align_and_qc_report.py --aligner star <i>project_name</i> <i>sample_info_file.txt</i>
	
This script uses the many output files created in step 1), converts these sample-specific files into matrices that include data for all samples, and then creates an Rmd document (main template is rnaseq_align_and_qc_report_Rmd_template.txt) that is converted into an html report using pandoc and R package rmarkdown. The report and accompanying files are contained in:

> <i>path_start</i>/<i>project_name</i>_Alignment_QC_Report/

The report can be opened with the file:

> <i>path_start</i>/<i>project_name</i>_Alignment_QC_Report/<i>project_name</i>_QC_RnaSeqReport.html

### Differential expression analysis - DESeq2

Run rnaseq_de_report.py to perform DE analysis and create an HTML report of differential expression summary statistics.

> rnaseq_de_report.py --project_name <i>output_prefix</i> --sample_in <i>sample_info_file_withQC.txt</i> --comp <i>sample_comp_file.txt</i> --de_package deseq2 --ref_genome <i>reference genome</i> --path_start <i>output_path</i> --template_dir <i>templete_file_directory</i>

The "--de_package deseq2" option is needed in order to specify that differential expression analysis should be conducted with DESeq2 using the HTSeq output file created after running rnaseq_align_and_qc.py. A merged transcriptome can be created using these files (option --merge_transcriptme yes), but the default is to use the reference genome gtf file. Comparisons of interest must be specified in a tab-delimited text file with one comparison per line, where comparison conditions are separated by "_vs_", as in "case_vs_control." A phenotype file must be provided; this file must at minimum indicate a phenotype for each sample; it may also contain other information, such as batch or cell line, if these are of interest, as per experiment design. 

The differential expression analysis accommodates a paired and unpaired option; if "paired" is selected, additionally specify the condition to pair by - e.g. paired:Cell_Line - note that the condition (in this case Cell_Line) must be a column name in the phenotype file. If there are any samples without a pair in any given comparison, the script will automatically drop these samples from that comparison, which will be noted in the report.

This script creates an Rmd document (main template is rnaseq_de_report_Rmd_template.txt) that uses the DESeq2 R package to load and process the HTSeq output file 2). The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_DE_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_DE_Report/<i>project_name</i>_DE_RnaSeqReport.html
	
### Differential expression analysis - Kallisto/Sleuth

Run kalliso_analysis.py to perform read pseudoalignment and quantification

> kallisto_analysis.py --path_start <i>project_directory_location</i> <i>sample_info_file.txt</i>

Run sleuth_analysis.py to perform differential expression analysis

> sleuth_analysis.py --path_start <i>project_directory_location</i>  --comp <i>project_comparison_file</i> <i>project_name</i> <i>sample_info_file.txt</i>

Run rnaseq_de_report.py to create an HTML report of differential expression summary statistics and plots for top differentially expressed transcripts according to all specified pairwise conditions

> python rnaseq_de_report.py <i>project_name</i> <i>sample_info_file.txt</i> --path_start <i>project_directory_location</i> --de sleuth --comp <i>project_comparison_file</i>

	
### acknowledgements
This set of scripts was initially developed to analyze RNA-Seq and DGE data at the [Partners Personalized Medicine PPM](http://pcpgm.partners.org/). Barbara Klanderman is the molecular biologist who led the establishment of PPM RNA-seq lab protocols and played an essential role in determining what components of the reports would be most helpful to PPM wet lab staff. 

