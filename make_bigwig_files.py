#!/usr/bin/python
import sys
import subprocess
import os
import math
import string
from rnaseq_align_and_qc import *


def make_track_file(sample_info_file, path_start, aligner):
	runs = get_sample_info(sample_info_file)
	outp = open(path_start+"ucsc_track_names.txt", "w")
	output_url = "https://dl.dropboxusercontent.com/u/178614/"
	
	#Get a dictionary of colors for conditions in info_sheet
	colors = ["0,0,400", "400,0,0", "0,0,100", "100,0,0", "200,200,0", "0,200,200"]
	conditions = list(set(map(lambda x: x[5], runs)))
	color_dic = {}
	for i in range(len(conditions)):
		color_dic[conditions[i]] = colors[i]

	#Write out each line of the track names file
	for k in runs:
		curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type, lane, run = k
		curr_bigwig = output_url+project+"/"+aligner+"_bigwig/"+curr_sample+"_accepted_hits.bigwig"

		outp.write("track type=bigWig name=\"")
		outp.write(label+"_"+curr_sample+"\" ")
		outp.write("color="+color_dic[label])
		#outp.write(" gridDefault=on maxHeightPixels=50 visibility=full autoScale=off viewLimits=0:13000 description=")
		outp.write(" gridDefault=on maxHeightPixels=50 visibility=full autoScale=on description=")
		outp.write(label+"_"+curr_sample+"\" bigDataUrl=")
		outp.write(curr_bigwig+"\n")
	outp.close()
	print "Created text file to load tracks in UCSC genome browser"





def main(project_name, aligner, sample_info_file, path_start):
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

	bedgraph_out = out_dir+aligner+"_bedgraph/"
	if not os.path.exists(bedgraph_out):
		os.makedirs(bedgraph_out)
	
	bigwig_out = out_dir+aligner+"_bigwig/"
	if not os.path.exists(bigwig_out):
		os.makedirs(bigwig_out)
		
	make_track_file(sample_info_file, bigwig_out, aligner)

	#Get sample info. Fields: [customer_id, gigpad_id, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	for k in runs:
		curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type, lane, run = k
		curr_bam = path_start+project+"/"+curr_sample+"/"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam"
		job_name = curr_sample
		
		if ref_genome == "hg38":
			chr_size_file = "/project/bhimeslab/Reference/hg38/hg38_ERCC.chrom.sizes"
		elif ref_genome == "mm38":
			chr_size_file = "/project/bhimeslab/Reference/mm38/mm38_ERCC.chrom.sizes"		
		elif ref_genome == "mm10":
			chr_size_file = "/project/bhimeslab/Reference/mm38/mm10_ERCC.chrom.sizes"	
		elif ref_genome == "rn6":
			chr_size_file = "/project/bhimeslab/Reference/rn6/rn6_ERCC.chrom.sizes"			
		elif ref_genome == "susScr3":
			chr_size_file = "/project/bhimeslab/Reference/susScr3/susScr3_ERCC.chrom.sizes"	
		else:
			print "No chromosome size file for this genome available:", ref_genome
		
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
		#Create bedgraph files
		#outp.write("genomeCoverageBed -split -bg -ibam "+curr_bam+" -g "+chr_size_file+" -strand + > "+bedgraph_out+curr_sample+"_accepted_hits.POS.bedgraph\n")
		#outp.write("genomeCoverageBed -split -bg -ibam "+curr_bam+" -g "+chr_size_file+" -strand - > "+bedgraph_out+curr_sample+"_accepted_hits.REV.bedgraph\n")
		outp.write("genomeCoverageBed -split -bg -ibam "+curr_bam+" -g "+chr_size_file+" > "+bedgraph_out+curr_sample+"_accepted_hits.bedgraph\n")
		
		#Sort bedgraph file
		outp.write("LC_COLLATE=C sort -k1,1 -k2,2n "+bedgraph_out+curr_sample+"_accepted_hits.bedgraph > "+bedgraph_out+curr_sample+"_accepted_hits.sorted.bedgraph\n")

		#Create bigwig giles
		#outp.write("./bedGraphToBigWig "+bedgraph_out+curr_sample+"_accepted_hits.POS.bedgraph "+chr_size_file+" "+bigwig_out+curr_sample+"_accepted_hits.POS.bigwig\n")
		#outp.write("./bedGraphToBigWig "+bedgraph_out+curr_sample+"_accepted_hits.REV.bedgraph "+chr_size_file+" "+bigwig_out+curr_sample+"_accepted_hits.REV.bigwig\n")
		outp.write("bedGraphToBigWig "+bedgraph_out+curr_sample+"_accepted_hits.sorted.bedgraph "+chr_size_file+" "+bigwig_out+curr_sample+"_accepted_hits.bigwig\n")

		#outp.write("rm "+bedgraph_out+curr_sample+"_accepted_hits.bedgraph\n")
		outp.close()
		#subprocess.call("bsub < "+job_name+".lsf", shell=True)
		#subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Creates bedgraph and bigwig files for RNA-seq samples associated with a project.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory resides (default=./)")
	#Change this to include chr size files into main directories with ref files (see main too)
	#parser.add_argument("--chr_size_file", default="/project/bhimeslab/Reference/hg38/hg38_ERCC.chrom.sizes", type=str, help="A file created using UCSC browser script fetchChromSizes (e.g. fetchChromSizes hg38 > hg38.chrom.sizes). Default: /project/bhimeslab/Reference/hg38/hg38_ERCC.chrom.sizes")
	parser.add_argument("--aligner", default="tophat", type=str, help="Should TopHat, STAR, or bowtie2 be used as aligner?"
		"(options: tophat, star, bowtie2)")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.aligner, args.samples_in, args.path_start)



