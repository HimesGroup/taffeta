#!/usr/bin/python
import sys
import subprocess
import os
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
	

def main(project_name, sample_info_file, path_start):
	#Get list of sample info. Fields: [customer_id, gigpad_id, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	ref_genome_list = map(lambda x: x[8], runs)
	#Check whether all samples are of same reference genome
	if False in map(lambda y: y==ref_genome_list[0], ref_genome_list):
		print "Make sure all samples in project are of the same reference genome"
		sys.exit()
	else:
		ref_genome = ref_genome_list[0]
		ref_index, fa, gtf, ref, ERCC_gtf = get_genome_ref_files(ref_genome)

	#Create the report
	if not os.path.exists("rnaseq_de_report_Rnw_template.txt"):
		print "Cannot find rnaseq_de_report_Rnw_template.txt"
		sys.exit()
	rnw_in = open("rnaseq_de_report_Rnw_template.txt", "r")
	rnw_template = rnw_in.read()
	if path_start == "./":
		path_start = os.getcwd()+"/"
	make_de_rnw_html(rnw_template, project_name, path_start, gtf, ref_genome)	


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Create HTML report of differential expression results for RNA-seq samples associated with a project.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory created by rnaseq_de.py is located (default=./)")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to. Often a PCPGM batch, but it could correspond to a mixture of batches.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.path_start)


