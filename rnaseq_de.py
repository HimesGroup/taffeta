#!/usr/bin/python
import sys
import subprocess
import os
from rnaseq_align_and_qc import *


def get_bam_files_for_group(group, runs, path_start):
	"""
	Gets a list of all bam files for a condition group from the original list of samples derived from get_sample_info()
	Used to create file list for cuffdiff
	"""
	sample_bam_list = []
	for k in runs:
		curr_sample, gigpad, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type = k
		if label == group:
			sample_bam_list.append(path_start+batch+"/"+curr_sample+"/tophat_out/accepted_hits.bam")
	return ",".join(sample_bam_list)


def make_assembly(assemblies_name, runs, path_start):
	"""
	Make a txt file of all transcripts.gtf files output by cufflinks for a set of samples derived from get_sample_info()
	Used to create assembly file for cuffmerge
	"""
	outp = open(assemblies_name, "w")
	for k in runs:
		curr_sample, gigpad, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type = k
		outp.write(path_start+batch+"/"+curr_sample+"/cufflinks_out/transcripts.gtf\n")
	outp.close()


def main(project_name, sample_info_file, merge_transcriptome, bias_correction, conditions, path_start):
	"""
	Dispatches an lsf job to 
	1) Possibly create a merged transcriptome from samples
	2) Run cuffdiff for all samples from a project
	"""
	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"
	out_dir = path_start+project_name+"/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	#Get sample info. Fields: [customer_id, gigpad_id, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	sample_names = map(lambda x: x[0], runs)
	condition_list = map(lambda x: x[7], runs)
	ref_genome_list = map(lambda x: x[8], runs)
	library_type_list = map(lambda x: x[9], runs)
	#Check whether all samples are of same reference genome
	if False in map(lambda y: y==ref_genome_list[0], ref_genome_list):
		print "Make sure all samples in project are of the same reference genome"
		sys.exit()
	else:
		ref_genome = ref_genome_list[0]
		ref_index, fa, gtf, ref, ERCC_gtf = get_genome_ref_files(ref_genome)
	#Check whether all samples are of same library type
	if False in map(lambda y: y==library_type_list[0], library_type_list):
		print "Make sure all samples in project are of the same library type"
		sys.exit()
	else:
		library_type = library_type_list[0]
	#Get condition list if none is supplied
	if conditions == "":
		conditions = sorted(set(condition_list))
	else:
		conditions = conditions.split(',')
		for c in conditions:
			if c not in condition_list:
				print "A condition supplied does not match those in sample file: ", c
				sys.exit()

	job_name = project_name+"_DE"

	#Make lsf file		
	outp = open(job_name+".lsf", "w")
	outp.write("#!/bin/bash \n")
	outp.write("#BSUB -L /bin/bash\n")
	outp.write("#BSUB -J "+job_name+"\n")
	outp.write("#BSUB -q big-multi \n")
	outp.write("#BSUB -o "+job_name+"_%J.out\n")
	outp.write("#BSUB -e "+job_name+"_%J.screen\n")
	outp.write("#BSUB -R 'rusage[mem=24000]'\n")

	if merge_transcriptome == "yes":
		#Create merged transcriptome with cuffmerge
		assemblies_name = out_dir+project_name+"_assemblies.txt"
		make_assembly(assemblies_name, runs, path_start)
		outp.write("cuffmerge -o "+out_dir+"merged_asm -g "+gtf+" -s "+fa+" -p 12 "+assemblies_name+"\n")
		cuffdiff_gtf = out_dir+"merged_asm/merged.gtf"		
	elif merge_transcriptome == "no":
		cuffdiff_gtf = gtf
		print "Will use "+ref_genome+" as reference transcriptome with cuffdiff."
	else:
		print "Merge transcriptome options are 'yes' or 'no'."
		sys.exit()

	#Run cuffdiff on all sample files, with options for bias and DGE library type
	outp.write("cuffdiff -p 12 -o "+out_dir+" ")
	if library_type == "DGE":
		outp.write("--no-length-correction ")
	if bias_correction == "yes":
		outp.write("-b "+fa+" ")
	elif bias_correction != "no":
		print "Bias correction options are 'yes' or 'no'"
		sys.exit()
	outp.write("-L "+",".join(conditions)+" -u "+cuffdiff_gtf+" "+" ".join(map(lambda(x): get_bam_files_for_group(x, runs, path_start), conditions))+"\n")
	outp.close()
	
	subprocess.call("bsub < "+job_name+".lsf", shell=True)
	subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to run cuffdiff for RNA-seq samples associated with a project.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where PCPGM batch-level directories are located and report directory will be written (default=./)")
	parser.add_argument("--bias_correction", default="yes", type=str, help="Should cuffdiff be run with bias correction? (options: yes, no; default=yes)")
	parser.add_argument("--merge_transcriptome", default="no", type=str, help="Should a merged transcriptome for samples be created for use with cuffdiff? (options: yes, no; default=no)")
	parser.add_argument("--conditions", type=str, default="", help="You can supply an ordered list of comma-separated conditions (e.g. cond1,cond2). If none given, conditions will be determined from sample file.")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to. Often a PCPGM batch, but it could correspond to a mixture of batches.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.merge_transcriptome, args.bias_correction, args.conditions, args.path_start)

