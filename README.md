taffeta
=======

RNAseq/DGE analysis pipeline based primarily on Tuxedo Tools. 

This set of scripts was developed by Blanca Himes to analyze RNA-Seq and DGE data at the Partners HealthCare, Inc. Center for Personalized Genetic Medicine ([PCPGM](http://pcpgm.partners.org/)). The goal is to take fastq files (sequenced with Illumina HiSeq or MiSeq) associated with a "project" and:
  * Perform preliminary QC 
  * Align reads to a reference genome
  * Perform QC on aligned files
  * Create a report that can be used to verify that sequencing was successful and/or identify sample outliers
  * Perform differential expression of reads aligned to transcripts according to a given reference genome
  * Create a report that summarizes the differential expression results

Several freely available software packages are used to perform most of these steps (see below). The pipeline is currently implemented in the Partners HealthCare, Inc. [High Performance Computing (HPC)](http://rc.partners.org/hpc) environment, where jobs are dispatched using LSF. 

### dependencies
* The Tuxedo suite of tools and various other programs should be installed: bowtie or bowtie2, tophat, cufflinks, cummerbund, fastqc, trimmomatic, samtools, bamtools, picardtools. 
* Annotation files should be available: reference genome fasta, gtf, refFlat, and index files. We use ERCC spike-ins, so our reference files include ERCC mix transcripts. 
* For adapter trimming, we include Ilumina and Nextflex sequences as used in the PCPGM lab.
* The Python scripts make use of modules that include subprocess, os, argparse, sys.
* To create reports, R and various libraries should be available, including cummeRbund, R2HTML, RColorBrewer, xtable.

### input files
Before running the pipeline, characteristics of a set of fastq files for samples that are part of a project are described in a tab-delimited txt file containing the following fields:
```
	customer_ID		| ID given to sample by customer
	gigpad_ID		| A sample ID that may differ from customer_ID and that is associated with library prior to pooling
	lane			| Lane of sequencer where this sample's pool was run, needed to locate raw fastq file
	index			| Six digit sequence of the library index used for this sample, needed to locate raw fastq file and perform adapter trimming
	ercc_mix		| Mix of ERCC spike used for library construction (options: "1", "2", "-")
	file_directory	| Top-level directory where sample's fastq files were written to following Casava filters
	batch			| gigpad batch number associated with sample, needed to locate raw fastq file and used as directory name for pipeline output files
	label			| Biological condition associated with the sample, provided by customer
	ref_genome		| Rerence genome associated with sample. (options: "hg19", "Zv9", "mm10")
	library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE",
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end")
```

PCPGM uses a rigid naming structure as is apparent in rnaseq_align_and_qc.py. The file naming and directory structure are obviously only applicable to local use of the scripts. They are being included for the sake of transparency and may someday be replaced with a more generalizable workflow. The fastq files that are associated with the project are read from where they are saved after sequencing/Casava filters are applied and then local copies are created using <i>sample_name</i>_R1.fastq (and <i>sample_name</i>_R2.fastq for paired reads). 

### workflow
1) Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project using rnaseq_align_and_qc.py:

> python rnaseq_align_and_qc.py --discovery no <i>sample_info_file.txt</i>

The "--discovery no" option refers to using --no-novel-juncs and --transcriptome-only options with tophat. Following execution of this script, various output files will be written for each sample in directories structured as:
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

2) Create an HTML report of QC and alignment summary statistics for RNA-seq samples associated with a project using rnaseq_align_and_qc_report.py:

> python rnaseq_align_and_qc_report.py <i>project_name</i> <i>sample_info_file.txt</i>
	
This script uses the many output files created in step 1), converts these sample-specific files into matrices that include data for all samples, and then creates an Rnw document (main template is rnaseq_align_and_qc_report_Rnw_template.txt) that is converted into an html report using R/Sweave. The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_Alignment_QC_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_Alignment_QC_Report/<i>project_name</i>_QC_RnaSeqReport.html

3) Write and execute an lsf job to perform differential expression analysis for RNA-seq samples associated with a project using rnaseq_de.py:

> python rnaseq_de.py <i>project_name</i> <i>sample_info_file.txt</i>

Differential expression analysis is conducted with cuffdiff using cufflinks output files created after running rnaseq_align_and_qc.py. A merged transcriptome can be created using these files (option --merge_transcriptme yes), but the default is to use the reference genome gtf file. If a particular order of conditions among samples is desired, it can be provided as a comma-separated list (option --conditions cond1,cond2,cond3,...). Otherwise, all condition types according to <i>sample_info_file.txt</i> sorted in alphabetical order are used.

4) Create an HTML report of differential expression summary statistics and plots for top differentially expressed genes according to all possible pairwise conditions for RNA-seq samples associated with a project using rnaseq_de_report.py:

> python rnaseq_de_report.py <i>project_name</i> <i>sample_info_file.txt</i>

This script creates an Rnw document (main template is rnaseq_de_report_Rnw_template.txt) that uses the cummeRbund R package to load and process the cuffdiff output files created in step 3). The report and accompanying files are contained in:

> <i>project_name</i>/<i>project_name</i>_DE_Report/

The report can be opened with the file:

> <i>project_name</i>/<i>project_name</i>_DE_Report/<i>project_name</i>_DE_RnaSeqReport.html
	
5) Create an [IGV](http://www.broadinstitute.org/igv/) batch script that can be used to create PDF snapshots of raw reads for differentially expressed genes (or any list of genes read from a txt file) using:

> python make_igv_script_file.py <i>project_name</i> <i>sample_info_file.txt</i>
	
The resulting output file:

> <i>project_name</i>/<i>project_name</i>_DE_Report/IGV_Plots/igv_snapshots_of_top_genes_script.txt
	
can be loaded from within IGV using File->Run Batch Script. Snapshots are saved to:

> <i>project_name</i>/<i>project_name</i>_DE_Report/IGV_Plots/
	
### acknowledgements
Barbara Klanderman is the molecular biologist who led the establishment of PCPGM RNA-seq lab protocols and played an essential role in determining what components of the reports would be most helpful to PCPGM wet lab staff. Thank you to Ken Auerbach and Jonathan Jackson of the Enterprise Research Infrastructure & Services (ERIS) group at Partners Healthcare for their in-depth support with installing and testing the programs whose output taffeta requires. Thank you to Rory Kirchner (@roryk) and Benjamin Harshfield for github-101 help and inspiration.
