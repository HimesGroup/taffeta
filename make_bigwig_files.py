#!/usr/bin/python
import sys
import subprocess
import os
import math
import string
from rnaseq_align_and_qc import *


def main(project_name, sample_info_file, path_start, chr_size_file):
	"""
	Creates bedgraph and bigwig files of all samples for a batch for visualization in UCSC genome browser
	"""
	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"
	out_dir = path_start+project_name+"/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	bedgraph_out = out_dir+"bedgraph/"
	if not os.path.exists(bedgraph_out):
		os.makedirs(bedgraph_out)
	
	bigwig_out = out_dir+"bigwig/"
	if not os.path.exists(bigwig_out):
		os.makedirs(bigwig_out)

	#Get sample info. Fields: [customer_id, gigpad_id, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	for k in runs:

		curr_sample, gigpad, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type = k
		curr_bam = path_start+batch+"/"+curr_sample+"/tophat_out/"+curr_sample+"_accepted_hits.sorted.bam"
		job_name = curr_sample
		
		#Make lsf file		
		outp = open(job_name+".lsf", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#BSUB -L /bin/bash\n")
		outp.write("#BSUB -J "+job_name+"\n")
		outp.write("#BSUB -q big-multi \n")
		outp.write("#BSUB -o "+job_name+"_%J.out\n")
		outp.write("#BSUB -e "+job_name+"_%J.screen\n")
		outp.write("#BSUB -R 'rusage[mem=24000]'\n")
		outp.write("#BSUB -n 12\n")

		#Create bedgraph files
		#outp.write("genomeCoverageBed -split -bg -ibam "+curr_bam+" -g "+chr_size_file+" -strand + > "+bedgraph_out+curr_sample+"_accepted_hits.POS.bedgraph\n")
		#outp.write("genomeCoverageBed -split -bg -ibam "+curr_bam+" -g "+chr_size_file+" -strand - > "+bedgraph_out+curr_sample+"_accepted_hits.REV.bedgraph\n")
		outp.write("genomeCoverageBed -split -bg -ibam "+curr_bam+" -g "+chr_size_file+" > "+bedgraph_out+curr_sample+"_accepted_hits.bedgraph\n")

		#Create bigwig giles
		#outp.write("./bedGraphToBigWig "+bedgraph_out+curr_sample+"_accepted_hits.POS.bedgraph "+chr_size_file+" "+bigwig_out+curr_sample+"_accepted_hits.POS.bigwig\n")
		#outp.write("./bedGraphToBigWig "+bedgraph_out+curr_sample+"_accepted_hits.REV.bedgraph "+chr_size_file+" "+bigwig_out+curr_sample+"_accepted_hits.REV.bigwig\n")
		outp.write("./bedGraphToBigWig "+bedgraph_out+curr_sample+"_accepted_hits.bedgraph "+chr_size_file+" "+bigwig_out+curr_sample+"_accepted_hits.bigwig\n")

		outp.close()
		subprocess.call("bsub < "+job_name+".lsf", shell=True)
		subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Creates bedgraph and bigwig files for RNA-seq samples associated with a project. Currently set for hg19 only.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory resides (default=./)")
	#Change this to include chr size files into main directories with ref files (see main too)
	parser.add_argument("--chr_size_file", default="/data/pcpgm/rnaseq/hg19_ERCC.chrom.sizes", type=str, help="A file created using UCSC browser script fetchChromSizes (e.g. fetchChromSizes hg19 > hg19.chrom.sizes). Default: /data/pcpgm/rnaseq/hg19.chrom.sizes")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to. Often a PCPGM batch, but it could correspond to a mixture of batches.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.path_start, args.chr_size_file)


