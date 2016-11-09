#!/usr/bin/python
import sys
import subprocess
import os
import argparse
import re
import fnmatch

def main(project_name, pheno_file, comp_file, path_start):
	"""
	Dispatches an lsf job to locate tsv files that were output by kallisto and runs sleuth analysis for all comparisons specified
	"""
	out_dir = path_start+"/"+project_name+"/sleuth_out/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	info_sheet = path_start+"/"+pheno_file
	#load text file containing all comparisons of interest
	comps_file = open(comp_file)
	comps = comps_file.readlines()
	for str in comps:
		str1 = re.split('&|\n', str)[0]
		str2 = re.split('&|\n', str)[1]
		sleuth_conds = [str1, str2]
		#sleuth needs cond1, cond2 to be in alphabetical order
		sleuth_conds.sort()
		job_name = project_name+"_sleuth_"+sleuth_conds[0]+"_vs_"+sleuth_conds[1]

		#make an R script to run for each comparison specified
		r_file = out_dir+job_name+".r"

		outp = open(r_file, "w")
		outp.write("#!/usr/bin/Rscript\n")
		outp.write("library(sleuth)\n")
		outp.write("library(biomaRt)\n")
		outp.write("cond1 <- '"+sleuth_conds[0]+"'\n")
		outp.write("cond2 <- '"+sleuth_conds[1]+"'\n")
		outp.write("project <- '"+project_name+"'\n")
		outp.write("base_dir <- paste0('"+path_start+"','/','"+project_name+"')\n")
		outp.write("info_sheet <- read.table('"+info_sheet+"', header = TRUE, stringsAsFactors=FALSE)\n")
		outp.write("cond <- paste('condition', cond2, sep='')\n")

		#select rows from info sheet based on chosen conditions
		#select sample_ids (needed to find kallisto directories)
		outp.write("s2c <- union(info_sheet[which(info_sheet$Label==cond1), ], info_sheet[which(info_sheet$Label==cond2), ])\n")
		outp.write("sample_id <- s2c$Sample\n")

		#directions to kallisto directories
		outp.write("kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'kallisto_output'))\n")
		outp.write("s2c <- dplyr::select(s2c, sample = Sample, condition = Label)\n")
		outp.write("s2c <- dplyr::mutate(s2c, path = kal_dirs)\n")
		outp.write("s2c$sample <- paste(s2c$condition, '_', s2c$sample, sep='')\n")
		outp.write("so <- sleuth_prep(s2c, ~ condition)\n")
		outp.write("so <- sleuth_fit(so)\n")
		outp.write("so <- sleuth_wt(so, cond)\n")

		#append official gene symbols
		outp.write("mart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host='www.ensembl.org')\n")
		#outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='grch37.ensembl.org', path='/biomart/martservice' ,dataset='hsapiens_gene_ensembl')\n")
		outp.write("t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id','external_gene_name'), mart = mart)\n")
		outp.write("t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)\n")

		#make into sleuth object
		outp.write("so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)\n")
		outp.write("so <- sleuth_fit(so)\n")
		outp.write("so <- sleuth_wt(so, which_beta = cond)\n")

		#save all kallisto output as one table - kallisto itself had saved each sample's output in a separate folder
		outp.write("kallisto_summary <- kallisto_table(so)\n")
		outp.write("kallisto_summary <- merge(kallisto_summary, t2g, by.x='target_id', by.y='target_id')\n")		
		outp.write("write.table(kallisto_summary, file=paste0(base_dir, '/sleuth_out/', project, '_', cond1, '_vs_', cond2, '_kallisto_output.txt'), row.names=FALSE, quote=FALSE)\n")
		
		#save RDS file for sleuth report which can be made using rnaseq_de_report.py
		outp.write("saveRDS(so, file = paste(base_dir, '/sleuth_out/so_', cond1, '_vs_', cond2, '.rds', sep=''))\n")
		
		#save down full sleuth results table
		outp.write("results_table <- sleuth_results(so, cond)\n")
		outp.write("write.table(results_table, file=paste(base_dir, '/sleuth_out/', project, '_sleuth_output_', cond1, '_vs_', cond2, '.txt', sep=''), row.names=FALSE, quote=FALSE)\n") 
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
		outp.close()
 
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform sleuth differential expression analysis based on kallisto files from kallisto_analysis.py.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
	parser.add_argument("--comp", help="A tab-delimited txt file containing sample comparisons to be made. One comparison per line, separate two conditions with & as in cond1&cond2.")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.comp, args.path_start)
