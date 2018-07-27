taffeta
=======

Reproducible analysis and validation of RNA-Seq data

Authors: Maya Shumyatcher, Mengyuan Kan, Blanca Himes.

The goal of taffeta is to preform reproducible analysis and validation of RNA-Seq data, as a part of [RAVED pipeline](https://github.com/HimesGroup/raved):
  * Download SRA .fastq data
  * Perform preliminary QC
  * Align reads to a reference genome
  * Perform QC on aligned files
  * Create a report that can be used to verify that sequencing was successful and/or identify sample outliers
  * Perform differential expression of reads aligned to transcripts according to a given reference genome
  * Create a report that summarizes the differential expression results

Several freely available software packages are used to perform most of these steps (see below).

### dependencies
* STAR, HTSeq, kallisto (older versions used bowtie2, tophat, cufflinks, cummerbund, whose options are no longer available).
* Programs that should be installed for QC: fastqc, trimmomatic, samtools, bamtools, picardtools.
* Annotation files should be available: reference genome fasta, gtf, refFlat, and index files. We use ERCC spike-ins, so our reference files include ERCC mix transcripts. 
* For adapter trimming, we include Ilumina TruSeq single index and Illumina unique dual (UD) index adpter and primer sequences.
* The Python scripts make use of modules that include subprocess, os, argparse, sys.
* To create reports, R and various libraries should be available, including DT, gplots, ggplot2, reshape2, rmarkdown, RColorBrewer, plyr, dplyr, lattice, ginefilter, biomaRt. Additionally, pandoc version 1.12.3 or higher should be available. If following the gene-based workflow, R package DESeq2 should be available. If following the workflow for transcript-based results, R package sleuth should be available. 

## Workflow

### Download data from GEO/SRA

Run script rnaseq_sra_download.py to download .fastq files from SRA. Users can provide a phenotype file with a SRA_ID column. Download samples corresponding to the SRA_ID. If a phenotype file is not provided, use phenotype information from GEO. SRA_ID is retrieved from the field relation.1. Corresponding ftp addresses are obtained from SRA sqilte database. Generate LSF scripts with download commands that can be processed in parallel. Specify --fastqc option to run FastQC for .fastq files downloaded.

> rnaseq_sra_download.py --geo_id <i>GEO_Accession</i> --path_start <i>output_path</i> --project_name <i>output_prefix</i> --template_dir <i>templete_file_directory</i> --fastqc

Output files: 1. <i>project_name</i>\_SRAdownload_RnaSeqReport.html; 2. <i>project_name</i>\_sraFile.info 3. <i>GEO_Accession</i>\_withoutQC.txt 4. FastQC results


### User-tailored phentoype file

The sample info file used in the following steps should be provided by users.
Required columns: 'Sample' column with sample ID, 'Status' column with comparison states.
Other columns: 'Treatment', 'Disease', 'Donor' (donor or cell line ID if in vitro treatment is used), 'Tissue', 'ERCC\_Mix' if ERCC spike-in sample is used, 'protocol' designating sample preparation kit information.
'Index' column has index sequence if used for each sample. For samples downloaded from SRA, the index information may not be provided in GEO. We suggest users reviewing "Overrepresented sequences" category from fastqc results and add corresponding index sequence to each sample.
Assign full paths for Read 1 and Read 2 (if available) .fastq files in R1 and R2 column. If multiple .fastq files obtained from sequencer, combine them into one file for Read 1 and Read 2 respectively.

### QC and Alignment




Steps 4-5 are optional.

If downloading data from GEO, start at "Step 0," else start at "Step 1." Additionally, if downloading data from GEO, the "Run" column in the sample info file should say "ncbi."

0) Download RNA-Seq data from GEO for all runs within a project using get_sras.py:

> get_sras.py --path_start /project/bhimeslab <i>sample_info_file.txt</i> <i>GEO_Accession</i>

1) Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project using rnaseq_align_and_qc.py:

> python rnaseq_align_and_qc.py --discovery no --aligner star <i>sample_info_file.txt</i>

The "--aligner star" option indicates that star should be used as the aligner (default is tophat). The "--discovery no" option refers to using --no-novel-juncs and --transcriptome-only options with tophat; while this option is not relevant if using star as the aligner, it must be specified in order for the script to run. Following execution of this script, various output files will be written for each sample in directories structured as:
> 
 <i>batch_num</i>/<i>sample_name</i>/tophat_out <br>
 <i>batch_num</i>/<i>sample_name</i>/cufflinks_out <br>
 <i>batch_num</i>/<i>sample_name</i>/cufflinks_out_ERCC <br>
 <i>batch_num</i>/<i>sample_name</i>/<i>sample_name</i>_R1_Trimmed.fastqc <br>
 <i>batch_num</i>/<i>sample_name</i>/<i>sample_name</i>_R2_Trimmed.fastqc <br>
 <i>batch_num</i>/<i>sample_name</i>/<i>sample_name</i>_R1_fastqc <br>
 <i>batch_num</i>/<i>sample_name</i>/<i>sample_name</i>_R2_fastqc <br>
 <i>batch_num</i>/<i>sample_name</i>/<i>sample_name</i>_ReadCount <br>
 ...

Note that the adapter trimming step is skipped for datasets downloaded from GEO, since Illumina barcodes for each sample are not provided in GEO.

2) Create an HTML report of QC and alignment summary statistics for RNA-seq samples associated with a project using rnaseq_align_and_qc_report.py:

> module load pandoc/2.0.6 <br>
> python rnaseq_align_and_qc_report.py --aligner star <i>project_name</i> <i>sample_info_file.txt</i>
	
This script uses the many output files created in step 1), converts these sample-specific files into matrices that include data for all samples, and then creates an Rmd document (main template is rnaseq_align_and_qc_report_Rmd_template.txt) that is converted into an html report using pandoc and R package rmarkdown. The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_Alignment_QC_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_Alignment_QC_Report/<i>project_name</i>_QC_RnaSeqReport.html

Note that pandoc version 1.12.3 or higher is required in order to use the rmarkdown package from the command line.

3) Perform differential expression analysis and create an HTML report of differential expression summary statistics and plots for top differentially expressed genes according to all specified pairwise conditions for RNA-seq samples associated with a project using rnaseq_de_report.py:

> module load pandoc/2.0.6 <br>
> python rnaseq_de_report.py <i>project_name</i> <i>sample_info_file.txt</i> --de_package deseq2 --comp <i>sample_comp_file.txt</i> --pheno <i>sample_pheno_file.txt</i> --design <i>[unpaired or paired:Condition]</i>

The "--de_package deseq2" option is needed in order to specify that differential expression analysis should be conducted with DESeq2 using the HTSeq output file created after running rnaseq_align_and_qc.py. A merged transcriptome can be created using these files (option --merge_transcriptme yes), but the default is to use the reference genome gtf file. Comparisons of interest must be specified in a tab-delimited text file with one comparison per line, where comparison conditions are separated by "_vs_", as in "case_vs_control." A phenotype file must be provided; this file must at minimum indicate a phenotype for each sample; it may also contain other information, such as batch or cell line, if these are of interest, as per experiment design. 

The differential expression analysis accommodates a paired and unpaired option; if "paired" is selected, additionally specify the condition to pair by - e.g. paired:Cell_Line - note that the condition (in this case Cell_Line) must be a column name in the phenotype file. If there are any samples without a pair in any given comparison, the script will automatically drop these samples from that comparison, which will be noted in the report.

This script creates an Rmd document (main template is rnaseq_de_report_Rmd_template.txt) that uses the DESeq2 R package to load and process the HTSeq output file 2). The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_DE_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_DE_Report/<i>project_name</i>_DE_RnaSeqReport.html
	
4) Create bigwig files of raw reads to be visualized on, e.g., the UCSC genome browser using:

> python make_bigwig_files.py <i>project_name</i> <i>sample_info_file.txt</i>
	
Requires the existence of a chromosome size file, which can be made using fetchChromSizes. Note that bigwig files 

5) Create a report of differentially expressed results for a given set of genes of interest. 

> python rnaseq_gene_subset_de_report.py <i>project_name</i> <i>sample_info_file.txt</i> <i>gene_list_file.txt</i>

where the <i>gene_list_file.txt</i> contains "gene_id" names matching those of the cuffdiff output file, one per line.

#### for transcript-based results

1) Write and execute an lsf job to perform read pseudoalignment and quantification for RNA-seq samples associated with a project using kallisto_analysis.py:

> kallisto_analysis.py --path_start <i>project_directory_location</i> <i>sample_info_file.txt</i>

2) Write and execute an lsf job to perform differential expression analysis for RNA-seq samples associated with a project using sleuth_analysis.py:

> sleuth_analysis.py --path_start <i>project_directory_location</i>  --comp <i>project_comparison_file</i> <i>project_name</i> <i>sample_info_file.txt</i>

3) Create an HTML report of differential expression summary statistics and plots for top differentially expressed transcripts according to all specified pairwise conditions for RNA-seq samples associated with a project using rnaseq_de_report.py:

> module load pandoc/2.0.6 <br>
> python rnaseq_de_report.py <i>project_name</i> <i>sample_info_file.txt</i> --path_start <i>project_directory_location</i> --de sleuth --comp <i>project_comparison_file</i>

	
### acknowledgements
This set of scripts was initially developed to analyze RNA-Seq and DGE data at the [Partners Personalized Medicine PPM](http://pcpgm.partners.org/). Barbara Klanderman is the molecular biologist who led the establishment of PPM RNA-seq lab protocols and played an essential role in determining what components of the reports would be most helpful to PPM wet lab staff. 

