#!/usr/bin/python
import sys
import subprocess
import os
import re
import fnmatch

from rnaseq_align_and_qc import *
from rnaseq_de import *


def make_de_rnw_html(rnw_template, project_name, path_start, gtf, ref_genome):
	"""
	Creates Rnw report. The top of report is below and the rest concatenated from a separate text document (rnw_template).
	Runs Sweave to create html document with R output (R2HTML library necessary for this)
	"""
	out_dir = path_start+project_name+"/"+project_name+"_DE_Report/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	outp = open(out_dir+project_name+"_DE_RnaSeqReport.Rnw", "w")
	outp.write("<HTML>\n")
	outp.write("<Head>\n")
	outp.write("<Title>\n")
	outp.write("RNA-seq Differential Expression Report for "+project_name+"\n")
	outp.write("</Title>\n")
	outp.write("</Head>\n")
	outp.write("<Body>\n\n")
	outp.write("<<input.data,echo=FALSE,results=hide>>=\n")
	outp.write("library(cummeRbund)\n")
	outp.write("curr.batch=\""+project_name+"\"\n")
	outp.write("path.start=\""+path_start+project_name+"\"\n")
	#cummerbund database can includes a reference genome and gtf file although these options are not necessary for report
	outp.write("cuff_data=readCufflinks(dir=path.start, genome = \""+ref_genome+"\", gtfFile = \""+gtf+"\", rebuild=TRUE) \n")
	outp.write("conditions=samples(genes(cuff_data)) \n")
	#reports currently focus on selecting top DE results based on genes and plots are made for genes and isoforms. 
	#it is possible to make reports that include TSS and CDS results, but these are not routinely given to customers
	outp.write("features=c(\"genes\") \n")
	outp.write("all.features=c(\"genes\", \"isoforms\") \n")
	#outp.write("all_features=c(\"genes\", \"isoforms\", \"TSS\", \"CDS\") \n")
	outp.write("ref.genome=\""+ref_genome+"\"\n")
	outp.write("@\n")
	outp.write(rnw_template)
	outp.close()
	subprocess.call("cd "+out_dir+"; echo \"library(R2HTML); Sweave('"+project_name+"_DE_RnaSeqReport.Rnw', driver=RweaveHTML)\" | R --no-save --no-restore", shell=True)



def read_htseq_file(fin):
	"""
	Read in output file created by htseq 
	Return as a matrix
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	#Removing final 5 rows of htseq output that contains summary, e.g., __no_feature __ambiguous etc
	c = [x for x in c if not re.match('__*', x)]
	htseq_out = map(lambda x: x.split('\t'), c)
	return htseq_out



def make_htseq_count_matrix(project_name, path_start, sample_info_file):
	"""
	Read in individual output files created by htseq out_dir+"/htseq_out/"+curr_sample+"_counts.txt according to sample_info_file
	Return as a single matrix for a whole batch
	"""
	#Make sure output directory for batch is in place
	out_dir = path_start+project_name+"/"+project_name+"_DE_Report/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	#Open output file
	outp1 = open(out_dir+project_name+"_htseq_matrix.txt", "w")
	name1 = ["Gene"]

	#Get sample info. Fields: [curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	sample_names = map(lambda x: x[0], runs)
	path_names = map(lambda x: x[3]+x[4], runs)
	
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = path_names[i]+"/"+curr_name+"/htseq_out/"+curr_name+"_counts.txt"
		#print curr_path
		if not os.path.exists(curr_path):
			print "Missing htseq output file ", curr_path
			break
		curr_htseq = read_htseq_file(curr_path)
		if len(curr_htseq) == 0:
			print "Current sample will not be included ", curr_name
			print "Empty htseq file ", curr_path
			print "Rerun with complete sample set"
			break
		name1.append(curr_name)
		if i == 0:
			a = curr_htseq
		else:
			#print len(a), len(curr_htseq)
			for j in range(len(a)):
				if a[j][0] == curr_htseq[j][0]:
					a[j] = a[j]+[curr_htseq[j][1]]
				else:
					print "Error, gene names don't match:", curr_name, a[j], curr_htseq[j]
	outp1.write("\t".join(name1)+"\n")
	outp1.write("\n".join(map("\t".join, a) ))
	outp1.close()
	print "Created read count matrix file "+out_dir+project_name+"_htseq_matrix.txt"


def make_deseq2_html(rmd_template, project_name, path_start, pheno_file, ref_genome, comp_file):
	"""
	Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
	"""
	out_dir = path_start+"/"+project_name+"/"+project_name+"_DE_Report/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	outp = open(out_dir+'/'+project_name+"_DESeq2_Report.Rmd", "w")

	#write .Rmd file
	outp.write("---\ntitle: \"DESeq2 results - based on HTSeq Counts from STAR-Aligned Sample Reads\"\n")
	outp.write("output: html_document\n---\n\n")
	outp.write("## DESeq2 results - based on HTSeq Counts from STAR-Aligned Sample Reads for "+project_name+"\n")
	outp.write("```{r, echo=FALSE, message=FALSE}\n")
	outp.write("library(gplots)\nlibrary(reshape2)\nlibrary(RColorBrewer)\nlibrary(plyr)\nlibrary(lattice)\nlibrary(genefilter)\nlibrary(ggplot2)\nlibrary(DESeq2)\nlibrary(biomaRt)\noptions(width = 1000)\n")
	outp.write("\n")
	outp.write("curr.batch=\""+project_name+"\"\n")
	outp.write("path.start='"+path_start+"'\n")
	if ref_genome=="hg38":
		outp.write("housekeeping_genes <- c('ACTB','GAPDH','B2M','RPL19','GABARAP')\n")
	elif ref_genome=="mm38" or ref_genome=="rn6":
		outp.write("housekeeping_genes <- c('Actb','Gapdh','B2m','Rpl19','Gabarap')\n")
	else:
		print("Housekeeping genes only available for hg38, mm38, rn6.")

	#Read in phenofile; phenofile should contain columns"sample", "condition" and "individual" 
	outp.write("coldata <- read.table(paste0(path.start,'"+pheno_file+"'), sep='\\t', header=TRUE)\n")
	outp.write("coldata <- subset(coldata, select=c(Sample,Label))\n")

	#load counts from HTSeq output 
	outp.write("countdata <- read.table(paste0(path.start, curr.batch, '/', curr.batch, '_DE_Report/', curr.batch, '_htseq_matrix.txt'), sep='\\t', header=TRUE)\n")
	outp.write("countdata$Gene <- sapply(strsplit(as.character(countdata$Gene), '\\\.'), '[[', 1)\n")
	outp.write("row.names(countdata) <- countdata$Gene\n")
	
	if ref_genome=="hg38":
		outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='hsapiens_gene_ensembl')\n")
		outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_gene_id', 'hgnc_symbol'), values=countdata$Gene, mart=mart)\n")
		outp.write("countdata <- merge(countdata, genes, by.x='Gene', by.y='ensembl_gene_id')\n")
		outp.write("countdata <- rename(countdata, c('hgnc_symbol'='gene_symbol'))\n")
	elif ref_genome=="mm38":
		outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='mmusculus_gene_ensembl')\n")
		outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_gene_id', 'mgi_symbol'), values=countdata$Gene, mart=mart)\n")
		outp.write("countdata <- merge(countdata, genes, by.x='Gene', by.y='ensembl_gene_id')\n")
		outp.write("countdata <- rename(countdata, c('mgi_symbol'='gene_symbol'))\n")
	elif ref_genome=="rn6":
		outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='rnorvegicus_gene_ensembl')\n")
		outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_gene_id', 'rgd_symbol'), values=countdata$Gene, mart=mart)\n")
		outp.write("countdata <- merge(countdata, genes, by.x='Gene', by.y='ensembl_gene_id')\n")
		outp.write("countdata <- rename(countdata, c('rgd_symbol'='gene_symbol'))\n")
	else:
		print("This code can only append official gene symbols for hg38, mm38, rn6.")
	outp.write("\n```\n")

	#load text file containing all comparisons of interest
	comps_file = open(comp_file)
	comps = comps_file.readlines()

	#create and paste the portion of the report that is unique to each comparison
	for str in comps:
		cond1 = re.split('&|\n', str)[0].rstrip()
		cond2 = re.split('&|\n', str)[1].rstrip()
		outp.write("```{r, echo=FALSE}\n")
		outp.write("cond1 <- '"+cond1+"'\n")
		outp.write("cond2 <- '"+cond2+"'\n")
		outp.write("\n```\n")
		outp.write("### "+cond1+" vs. "+cond2+" comparison\n")
		outp.write("\n")
		outp.writelines(rmd_template)
		outp.write("\n")

	#Housekeeping gene expression barplots - once per report
	outp.write("### Housekeeping genes\n")
	outp.write("Counts have been normalized by sequencing depth, with pseudocount of 0.5 added to allow for log scale plotting, using DESeq2 function plotCounts().\n")
	outp.write("\n")
	outp.write("```{r, echo=FALSE, cache=FALSE, warning=FALSE, message=FALSE}\n")
	outp.write("for (i in 1:length(housekeeping_genes)) {\n")
	outp.write("  gene_symbol <- housekeeping_genes[i]\n")
	outp.write("  gene <- res_df[which(res_df$gene_symbol==gene_symbol),]$Gene\n")
	outp.write("  dfs <- ls()[grep('_norm_data_', ls())]\n")
	outp.write("  dfs <- dfs[grep(gene_symbol, dfs)]\n")
	outp.write("  curr_data <- data.frame()\n")	
	outp.write("  for (j in dfs) {\n")
	outp.write("  	norm_data <- get(j)\n")
	outp.write("  	norm_data$Sample <- rownames(norm_data)\n")
	outp.write("  	rownames(norm_data) <- NULL\n")
	outp.write("  	curr_data <- rbind(curr_data, norm_data)\n")
	outp.write("  }\n")
	outp.write("  	curr_data <- unique(curr_data)\n")
	outp.write("  	curr_data$Sample <- NULL\n")
	outp.write("      housekeeping_plot <- ggplot(curr_data, aes(x = curr_data$Label, y = count, fill=curr_data$Label)) +\n")
	outp.write("      geom_boxplot(outlier.colour=NA, lwd=0.2, color='grey18') +\n")
	outp.write("      stat_boxplot(geom ='errorbar', color='grey18') +\n")
	outp.write("      geom_jitter(size=2, width=0.2) +\n")
	outp.write("      guides(fill=FALSE) +\n")
	outp.write("      theme_bw() +\n")
	outp.write("      labs(title=paste0(gene_symbol,'_',gene)) + labs(x='condition') + labs(y='counts') +\n")
	outp.write("      theme(strip.text.x = element_text(size = 10),axis.text.x = element_text(angle = 90, hjust = 1, size=12),axis.text.y = element_text(size=9),title = element_text(size=12), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))\n")	
	outp.write("	print(housekeeping_plot)\n")
	outp.write("}\n```\n")
	outp.close()
	subprocess.call("cd "+out_dir+"; echo \"library(knitr); library(markdown); knit2html('"+project_name+"_DESeq2_Report.Rmd', force_v1 = TRUE)\" | R --no-save --no-restore", shell=True)

def make_sleuth_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file):
	"""
	Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
	"""
	out_dir = path_start+"/"+project_name+"/"+project_name+"_DE_Report/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	sleuth_dir = path_start+"/"+project_name+"/sleuth_out/"
	outp = open(out_dir+project_name+"_Sleuth_Report.Rmd", "w")
	outp.write("---\ntitle: \"Sleuth results - based on Kallisto TPM\"\n")
	outp.write("output: html_document \ntoc: true \ntoc_depth: 2 \n---\n")
	outp.write("\n")
	outp.write("## Sleuth results - based on Kallisto TPM for "+project_name+"\n")
	outp.write("```{r set-options, echo=FALSE, cache=FALSE, warning=FALSE, message=FALSE}\n")
	outp.write("library(VennDiagram)\nlibrary(reshape2)\nlibrary(sleuth)\nlibrary(biomaRt)\noptions(width = 2000)\n")
	outp.write("\n")
	outp.write("curr.batch=\""+project_name+"\"\n")
	outp.write("path.start='"+path_start+"'\n")	
	outp.write("kallisto_output <- data.frame()\n")
	if ref_genome=="hg38":
		outp.write("housekeeping_genes <- c('ACTB','GAPDH','B2M','RPL19','GABARAP')\n")
	elif ref_genome=="mm38" or ref_genome=="rn6":
		outp.write("housekeeping_genes <- c('Actb','Gapdh','B2m','Rpl19','Gabarap')\n")
	else:
		print("Housekeeping genes only available for hg38, mm38, rn6.")

	#load info sheet 
	outp.write("info_sheet <- read.table(sample_info_file, header = TRUE, stringsAsFactors=FALSE)\n") #outp.write("info_sheet <- read.table(paste0(path.start,'"+sample_info_file+"'), header = TRUE, stringsAsFactors=FALSE)\n")
	outp.write("info_sheet <- subset(info_sheet, select=c('Sample','Label'))\n")
	outp.write("colnames(info_sheet) <- c('run_accession', 'condition')\n")
	outp.write("info_sheet <- info_sheet[order(info_sheet$condition),]\n")
	outp.write("print(info_sheet, row.names=FALSE)\n")
	outp.write("\n```\n")

	#load text file containing all comparisons of interest
	comps_file = open(comp_file)
	comps = comps_file.readlines()

	#create and paste the portion of the report that is unique to each comparison
	for str in comps:
		str1 = re.split('&|\n', str)[0]
		str2 = re.split('&|\n', str)[1]
		mylist = [str1, str2]
		cond1 = mylist[0].rstrip()
		cond2 = mylist[1].rstrip()
		#sleuth needs cond1, cond2 to be in alphabetical order
		mylist.sort()
		outp.write("```{r, echo=FALSE}\n")
		outp.write("cond1 <- paste0('condition','"+cond1+"')\n")
		outp.write("cond2 <- paste0('condition','"+cond2+"')\n")
		outp.write("\n```\n")
		outp.write("### "+cond1+" vs. "+cond2+" comparison\n")
		outp.write("```{r, echo=FALSE}\n")
		outp.write("so <- readRDS(paste0('"+sleuth_dir+"', 'so_', '"+cond1+"', '_vs_', '"+cond2+"', '.rds'))\n")
		outp.write("k_curr <- kallisto_table(so)\n")
		
		if ref_genome=="hg38":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='hsapiens_gene_ensembl')\n")
			outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), values=k_curr$target_id, mart=mart)\n")
			outp.write("k_curr <- merge(k_curr, genes, by.x='target_id', by.y='ensembl_transcript_id')\n")
			outp.write("names(k_curr)[names(k_curr)=='hgnc_symbol'] <-'gene_symbol'\n")
		elif ref_genome=="mm38":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='mmusculus_gene_ensembl')\n")
			outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_transcript_id', 'ensembl_gene_id', 'mgi_symbol'), values=k_curr$target_id, mart=mart)\n")
			outp.write("k_curr <- merge(k_curr, genes, by.x='target_id', by.y='ensembl_transcript_id')\n")
			outp.write("names(k_curr)[names(k_curr)=='mgi_symbol'] <-'gene_symbol'\n")
		elif ref_genome=="rn6":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='rnorvegicus_gene_ensembl')\n")
			outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_transcript_id', 'ensembl_gene_id', 'rgd_symbol'), values=k_curr$target_id, mart=mart)\n")
			outp.write("k_curr <- merge(k_curr, genes, by.x='target_id', by.y='ensembl_transcript_id')\n")
			outp.write("names(k_curr)[names(k_curr)=='rgd_symbol'] <-'gene_symbol'\n")
		else:
			print("This code can only append official gene symbols for hg38, mm38, rn6.")
		outp.write("kallisto_output <- rbind(kallisto_output, k_curr)\n")
		outp.write("kallisto_output <- unique(kallisto_output)\n")
		outp.write("```\n")
		outp.write("\n")
		outp.writelines(rmd_template)
		outp.write("\n")

	#Housekeeping gene expression barplots - once per report
	outp.write("### Housekeeping genes\n")
	outp.write("```{r, echo=FALSE}\n")
	outp.write("house_ensg <- unique(k_curr$ensembl_gene_id[which(k_curr$gene_symbol %in% housekeeping_genes)])\n")
	outp.write("for (i in 1:length(house_ensg)) {\n")
  	outp.write("  curr_gene <- house_ensg[i]\n")
  	outp.write("  curr_data <- kallisto_output[which(kallisto_output$ensembl_gene_id==curr_gene),]\n")
 	outp.write("  curr_data <- subset(curr_data, select=c(target_id, tpm, condition, gene_symbol))\n")
	outp.write("  curr_data$condition <- factor(curr_data$condition)\n") #, levels=c('Control_Baseline', 'Asthma_Baseline'))\n")
  	outp.write("  gene_symbol <- unique(curr_data$gene_symbol)\n")
	outp.write("  print({\n")
    	outp.write("	ggplot(curr_data, aes(x = condition, y = tpm, fill=condition)) + \n")
        outp.write("	geom_boxplot(outlier.colour=NA, lwd=0.2, color='grey18') + \n")
       	outp.write("	stat_boxplot(geom ='errorbar', color='grey18') + \n")
       	outp.write("	geom_jitter(size=0.8, width=0.2) + \n")
      	outp.write("	facet_wrap(~target_id) + \n")
       	outp.write("	guides(fill=FALSE) + \n")
       	outp.write("	theme_bw() +  \n")
       	outp.write("	labs(title=gene_symbol) + \n")
       	outp.write("	labs(x='condition') + labs(y='TPM') + \n")
       	outp.write("	theme(text = element_text(size=9), \n")
        outp.write("	  strip.text.x = element_text(size = 10), \n")
        outp.write("	  axis.text.x = element_text(angle = 90, hjust = 1, size=12),\n")
        outp.write("	  axis.text.y = element_text(size=9),\n")
        outp.write("	  title = element_text(size=12),\n")
        outp.write("	  axis.title.x = element_text(size=12),\n")
        outp.write("	  axis.title.y = element_text(size=12))})\n")
	outp.write("}\n")
	outp.write("```\n")
	outp.close()
	subprocess.call("cd "+out_dir+"; echo \"library(knitr); library(markdown); knit2html('"+project_name+"_Sleuth_Report.Rmd', force_v1 = TRUE)\" | R --no-save --no-restore", shell=True)

def main(project_name, sample_info_file, de_package, pheno_file, path_start, comp_file):
	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"

	#Get sample info. Fields: [curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	ref_genome_list = map(lambda x: x[6], runs)
	condition_list = map(lambda x: x[5], runs)
	#Check whether all samples are of same reference genome
	if False in map(lambda y: y==ref_genome_list[0], ref_genome_list):
		print "Make sure all samples in project are of the same reference genome"
		sys.exit()
	else:
		ref_genome = ref_genome_list[0]
		ref_index, fa, gtf, ref, ERCC_gtf, mask_gtf, genome_dir = get_genome_ref_files(ref_genome, "NA")

	if de_package == "cummerbund":
		#Create the report
		if not os.path.exists("/project/bhimeslab/taffeta/rnaseq_de_report_Rnw_template.txt"):
			print "Cannot find rnaseq_de_report_Rnw_template.txt"
			sys.exit()
		rnw_in = open("/project/bhimeslab/taffeta/rnaseq_de_report_Rnw_template.txt", "r")
		rnw_template = rnw_in.read()
		make_de_rnw_html(rnw_template, project_name, path_start, gtf, ref_genome)	

	elif de_package == "deseq2":
		#Make htseq count matrix and save it to project directory
		make_htseq_count_matrix(project_name, path_start, sample_info_file)

		#Create the report
		if not os.path.exists("/project/bhimeslab/taffeta/rnaseq_deseq2_Rmd_template.txt"):
			print "Cannot find rnaseq_deseq2_Rmd_template.txt"
			sys.exit()
		rmd_in = open("/project/bhimeslab/taffeta/rnaseq_deseq2_Rmd_template.txt", "r")
		rmd_template = rmd_in.readlines()
		make_deseq2_html(rmd_template, project_name, path_start, pheno_file, ref_genome, comp_file)
	
	if de_package == "sleuth":
		#Create the report
		if not os.path.exists("/project/bhimeslab/taffeta/rnaseq_sleuth_Rmd_template.txt"):
			print "Cannot find rnaseq_sleuth_Rmd_template.txt"
			sys.exit()
		rmd_in = open("/project/bhimeslab/taffeta/rnaseq_sleuth_Rmd_template.txt", "r")
		rmd_template = rmd_in.readlines()
		make_sleuth_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file)	
		
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Create HTML report of differential expression results for RNA-seq samples associated with a project.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory created by rnaseq_de.py is located (default=./)")
	parser.add_argument("--de", default="deseq2", type=str, help="Should cummeRbund, DESeq2 or sleuth be used for diferential expression (DE) analysis?  If sleuth, request an interactive node with more memory using: bsub -Is -M 36000 bash"
		"(options: deseq2, cummerbund)")
	parser.add_argument("--pheno", help="A tab-delimited txt file containing sample PHENOTYPE information. Make sure sample IDs match and are in the same order as those in the sample information file.")
	parser.add_argument("--comp", help="A tab-delimited txt file containing sample comparisons to be made. One comparison per line, separate two conditions with & as in cond1&cond2.")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.de, args.pheno, args.path_start, args.comp)
