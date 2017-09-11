#!/usr/bin/python
import sys
import subprocess
import os
import argparse
import re
import fnmatch

from rnaseq_align_and_qc import *
from rnaseq_de import *


def main(project_name, pheno_file, comp_file, path_start):
	"""
	Dispatches an lsf job to locate tsv files that were output by kallisto and runs sleuth analysis for all comparisons specified
	"""
	out_dir = path_start+"/"+project_name+"/sleuth_out/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	info_sheet = path_start+"/"+pheno_file

	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"

	#Get sample info. Fields: [curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type]
	runs = get_sample_info(info_sheet)
	ref_genome_list = map(lambda x: x[6], runs)

	#Check whether all samples are of same reference genome
	if False in map(lambda y: y==ref_genome_list[0], ref_genome_list):
		print "Make sure all samples in project are of the same reference genome"
		sys.exit()
	else:
		ref_genome = ref_genome_list[0]
		ref_index, fa, gtf, ref, ERCC_gtf, mask_gtf, genome_dir = get_genome_ref_files(ref_genome, "NA")

	#load text file containing all comparisons of interest
	comps_file = open(comp_file)
	comps = comps_file.readlines()
	for str in comps:
		#note sleuth needs case and control to be in alphabetical order
		#i.e. the thing that comes first alphabetically will become the control
		#the first string in comps_file is the one you want to be 'case'; hence, add 'zz_' to it

		str1 = 'zz_' + re.split('_vs_', str)[0]
		str2 = re.split('_vs_', str)[1]
		mylist = [str1, str2]
		mylist.sort()
		ctrl = mylist[0].rstrip() #control will be the one that comes first alphabetically
		case = mylist[1].rstrip() #case wil be the one that starts with 'zz_' and hence is second in the sorted list

		job_name = project_name+"_sleuth_"+case+"_vs_"+ctrl

		#make an R script to run for each comparison specified
		r_file = out_dir+job_name+".R"

		outp = open(r_file, "w")
		outp.write("#!/usr/bin/Rscript\n")
		outp.write("library(sleuth)\n")
		outp.write("library(biomaRt)\n")
		outp.write("case <- '"+case+"'\n")
		outp.write("ctrl <- '"+ctrl+"'\n")
		outp.write("project <- '"+project_name+"'\n")
		outp.write("base_dir <- paste0('"+path_start+"','/','"+project_name+"')\n")
		outp.write("info_sheet <- read.table('"+info_sheet+"', header = TRUE, stringsAsFactors=FALSE)\n")
		outp.write("cond <- paste('condition', case, sep='')\n")

		#select rows from info sheet based on chosen conditions
		#select sample_ids (needed to find kallisto directories)
		#paste0('x_',cond) option added in case need to add x_ to condition names to force the comparison to go in a certain direction (i.e. so control condition is the intercept)
		outp.write("s2c <- info_sheet[which(info_sheet$Label %in% c(case, ctrl) | paste0('zz_',info_sheet$Label) %in% c(case, ctrl)), ]\n")
		outp.write("for (i in 1:nrow(s2c)) {\n")
		outp.write("  if (s2c$Label[i] %in% c(case, ctrl)) {s2c$Label[i] <- s2c$Label[i]} \n")
		outp.write("  else {s2c$Label[i] <- paste0('zz_', s2c$Label[i])}\n")
		outp.write("}\n")
		outp.write("sample_id <- s2c$Sample\n")

		#directions to kallisto directories
		outp.write("kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'kallisto_output'))\n")
		outp.write("s2c <- dplyr::select(s2c, sample = Sample, condition = Label)\n")
		outp.write("s2c <- dplyr::mutate(s2c, path = kal_dirs)\n")
		outp.write("s2c$sample <- paste(s2c$condition, '_', s2c$sample, sep='')\n")

		#append official gene symbols 
		if ref_genome=="hg38":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='hsapiens_gene_ensembl')\n")
			outp.write("t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id','external_gene_name'), mart = mart)\n")
			outp.write("t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)\n")
		elif ref_genome=="mm38" or ref_genome=="mm10":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='mmusculus_gene_ensembl')\n")
			outp.write("t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id','mgi_symbol'), mart = mart)\n")
			outp.write("t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = mgi_symbol)\n")
		elif ref_genome=="rn6":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='rnorvegicus_gene_ensembl')\n")
			outp.write("t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id','rgd_symbol'), mart = mart)\n")
			outp.write("t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = rgd_symbol)\n")
		else:
			print("This code can only append official gene symbols for hg38, mm38, mm10 or rn6.")
		
		#make into sleuth object
		outp.write("so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)\n")
		outp.write("so <- sleuth_fit(so)\n")
		outp.write("so <- sleuth_wt(so, which_beta = cond)\n")

		#save all kallisto output as one table - kallisto itself had saved each sample's output in a separate folder
		outp.write("kallisto_summary <- kallisto_table(so)\n")
		outp.write("kallisto_summary$target_id <- sapply(strsplit(as.character(kallisto_summary$target_id), '\\\.'), '[[', 1)\n")
		outp.write("kallisto_summary <- merge(kallisto_summary, t2g, by.x='target_id', by.y='target_id')\n")		
		outp.write("write.table(kallisto_summary, file=paste0(base_dir, '/sleuth_out/', project, '_', case, '_vs_', ctrl, '_kallisto_output.txt'), row.names=FALSE, quote=FALSE)\n")
		
		#save RDS file for sleuth report which can be made using rnaseq_de_report.py
		outp.write("saveRDS(so, file = paste(base_dir, '/sleuth_out/so_', case, '_vs_', ctrl, '.rds', sep=''))\n")
		
		#save down full sleuth results table
		outp.write("results_table <- sleuth_results(so, cond)\n")
		outp.write("results_table$target_id <- sapply(strsplit(as.character(results_table$target_id), '\\\.'), '[[', 1)\n") 
		outp.write("results_table <- results_table[,which(!colnames(results_table) %in% c('ens_gene', 'ext_gene'))]\n")
		outp.write("results_table <- merge(results_table, t2g, by.x='target_id', by.y='target_id')\n")
		outp.write("write.table(results_table, file=paste(base_dir, '/sleuth_out/', project, '_', case, '_vs_', ctrl, '_sleuth_output.txt', sep=''), row.names=FALSE, quote=FALSE)\n") 
		outp.close()

		#Make lsf file		
		outp = open(job_name+".lsf", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#BSUB -L /bin/bash\n")
		outp.write("#BSUB -J "+job_name+"\n")
		outp.write("#BSUB -q normal \n")
		outp.write("#BSUB -o "+job_name+"_%J.out\n")
		outp.write("#BSUB -e "+job_name+"_%J.screen\n")
		outp.write("#BSUB -M 36000\n")
		outp.write("#BSUB -n 12\n")
		outp.write("Rscript --no-save --no-restore "+r_file+"")
		outp.write("\n")
		outp.close()
	print("Executable files have been created.")
 
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform sleuth differential expression analysis based on kallisto files from kallisto_analysis.py.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
	parser.add_argument("--comp", help="A tab-delimited txt file containing sample comparisons to be made. One comparison per line, separate two conditions with & as in case&ctrl.")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.comp, args.path_start)
