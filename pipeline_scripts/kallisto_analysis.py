#!/usr/bin/python
import subprocess
import os
import argparse


def get_sample_info(fin):
	"""
	Open tab-delimited txt file containing the following columns:
	v0: sample_ID		| ID given to sample 
	v1: index			| Six digit sequence of the index for this library
	v2: ercc_mix		| Mix of ERCC spike used for library construction (options: "1", "2", "-")
	v3: file_directory	| Directory where sample's fastq files reside
	v4: project			| Name for project associated with sample
	v5: label			| Biological condition associated with the sample, provided by customer
	v6: ref_genome		| Rerence genome associated with sample. (options: "hg38")
	v7: library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE",
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end")
	v8: lane			| Lane of sequencer (needed with UPenn NGSC files that are named by sample/lane/barcode)
	v9: run				| Run number on sequencer (needed with UPenn NGSC files that are named by sample/lane/barcode)
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	d = []
	for x in c:
		sample_id = x.split('\t')[0]
		index = x.split('\t')[1]
		ercc_mix = x.split('\t')[2]
		top_dir = x.split('\t')[3]
		if top_dir[-1] != "/":
			top_dir = top_dir+"/"
		project = x.split('\t')[4]
		label = x.split('\t')[5]
		ref_genome = x.split('\t')[6]
		library_type = x.split('\t')[7]
		lane = x.split('\t')[8]
		run = x.split('\t')[9]
		d.append([sample_id, index, ercc_mix, top_dir, project, label, ref_genome, library_type, lane, run])

	return d


def get_genome_ref_files(genome):
	"""
	Location of all reference files needed for a given genome.
	The ERCC gtf files were appended separately to each species own gtf file
	Current choice: "hg38"
	"""
	ERCC_index = "/project/bhimeslab/Reference/ERCC92.idx"
	if genome == "hg38":
		ref_index = "/project/bhimeslab/Reference/hg38/hg38_new.idx"
	elif genome == "mm38" or genome == "mm10":
		ref_index = "/project/bhimeslab/Reference/mm38/mm38.idx"
	else:
		print 'Unknown genome selected: ', genome
	return(ref_index, ERCC_index)


#Index sequencing primers
primer_start = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
primer_end = "ATCTCGTATGCCGTCTTCTGCTTG"

#Indexes used for PE, SE, or SPE RNA-seq reads
illumina_indexes = \
	{'ATCACG':['TruSeqAdapterIndex1', primer_start+'ATCACG'+primer_end], \
	 'CGATGT':['TruSeqAdapterIndex2', primer_start+'CGATGT'+primer_end], \
	 'TTAGGC':['TruSeqAdapterIndex3', primer_start+'TTAGGC'+primer_end], \
	 'TGACCA':['TruSeqAdapterIndex4', primer_start+'TGACCA'+primer_end], \
	 'ACAGTG':['TruSeqAdapterIndex5', primer_start+'ACAGTG'+primer_end], \
	 'GCCAAT':['TruSeqAdapterIndex6', primer_start+'GCCAAT'+primer_end], \
	 'CAGATC':['TruSeqAdapterIndex7', primer_start+'CAGATC'+primer_end], \
	 'ACTTGA':['TruSeqAdapterIndex8', primer_start+'ACTTGA'+primer_end], \
	 'GATCAG':['TruSeqAdapterIndex9', primer_start+'GATCAG'+primer_end], \
	 'TAGCTT':['TruSeqAdapterIndex10',primer_start+'TAGCTT'+primer_end], \
	 'GGCTAC':['TruSeqAdapterIndex11',primer_start+'GGCTAC'+primer_end], \
	 'CTTGTA':['TruSeqAdapterIndex12',primer_start+'CTTGTA'+primer_end], \
	 'AGTCAA':['TruSeqAdapterIndex13',primer_start+'AGTCAACA'+primer_end], \
	 'AGTTCC':['TruSeqAdapterIndex14',primer_start+'AGTTCCGT'+primer_end], \
	 'ATGTCA':['TruSeqAdapterIndex15',primer_start+'ATGTCAGA'+primer_end], \
	 'CCGTCC':['TruSeqAdapterIndex16',primer_start+'CCGTCCCG'+primer_end], \
	 'GTCCGC':['TruSeqAdapterIndex18',primer_start+'GTCCGCAC'+primer_end], \
	 'GTGAAA':['TruSeqAdapterIndex19',primer_start+'GTGAAACG'+primer_end], \
	 'GTGGCC':['TruSeqAdapterIndex20',primer_start+'GTGGCCTT'+primer_end], \
	 'GTTTCG':['TruSeqAdapterIndex21',primer_start+'GTTTCGGA'+primer_end], \
	 'CGTACG':['TruSeqAdapterIndex22',primer_start+'CGTACGTA'+primer_end], \
	 'GAGTGG':['TruSeqAdapterIndex23',primer_start+'GAGTGGAT'+primer_end], \
	 'ACTGAT':['TruSeqAdapterIndex25',primer_start+'ACTGATAT'+primer_end], \
	 'ATTCCT':['TruSeqAdapterIndex27',primer_start+'ATTCCTTT'+primer_end]}

#Indexes used for DGE reads
nextflex_indexes = \
	{'CGATGT':['NextFlexAdapter1', primer_start+'CGATGT'+primer_end], \
	 'TGACCA':['NextFlexAdapter2', primer_start+'TGACCA'+primer_end], \
	 'ACAGTG':['NextFlexAdapter3', primer_start+'ACAGTG'+primer_end], \
	 'GCCAAT':['NextFlexAdapter4', primer_start+'GCCAAT'+primer_end], \
	 'CAGATC':['NextFlexAdapter5', primer_start+'CAGATC'+primer_end], \
	 'CTTGTA':['NextFlexAdapter6', primer_start+'CTTGTA'+primer_end], \
	 'ATCACG':['NextFlexAdapter7', primer_start+'ATCACG'+primer_end], \
	 'TTAGGC':['NextFlexAdapter8', primer_start+'TTAGGC'+primer_end], \
	 'ACTTGA':['NextFlexAdapter9', primer_start+'ACTTGA'+primer_end], \
	 'GATCAG':['NextFlexAdapter10',primer_start+'GATCAG'+primer_end], \
	 'TAGCTT':['NextFlexAdapter11',primer_start+'TAGCTT'+primer_end], \
	 'GGCTAC':['NextFlexAdapter12',primer_start+'GGCTAC'+primer_end], \
	 'AGTCAA':['NextFlexAdapter13',primer_start+'AGTCAA'+primer_end], \
	 'AGTTCC':['NextFlexAdapter14',primer_start+'AGTTCC'+primer_end], \
	 'ATGTCA':['NextFlexAdapter15',primer_start+'ATGTCA'+primer_end], \
	 'CCGTCC':['NextFlexAdapter16',primer_start+'CCGTCC'+primer_end], \
	 'GTAGAG':['NextFlexAdapter17',primer_start+'GTAGAG'+primer_end], \
	 'GTCCGC':['NextFlexAdapter18',primer_start+'GTCCGC'+primer_end], \
	 'GTGAAA':['NextFlexAdapter19',primer_start+'GTGAAA'+primer_end], \
	 'GTGGCC':['NextFlexAdapter20',primer_start+'GTGGCC'+primer_end], \
	 'GTTTCG':['NextFlexAdapter21',primer_start+'GTTTCG'+primer_end], \
	 'CGTACG':['NextFlexAdapter22',primer_start+'CGTACG'+primer_end], \
	 'GAGTGG':['NextFlexAdapter23',primer_start+'GAGTGG'+primer_end], \
	 'GGTAGC':['NextFlexAdapter24',primer_start+'GGTAGC'+primer_end], \
	 'ACTGAT':['NextFlexAdapter25',primer_start+'ACTGAT'+primer_end], \
	 'ATGAGC':['NextFlexAdapter26',primer_start+'ATGAGC'+primer_end], \
	 'ATTCCT':['NextFlexAdapter27',primer_start+'ATTCCT'+primer_end], \
	 'CAAAAG':['NextFlexAdapter28',primer_start+'CAAAAG'+primer_end], \
	 'CAACTA':['NextFlexAdapter29',primer_start+'CAACTA'+primer_end], \
	 'CACCGG':['NextFlexAdapter30',primer_start+'CACCGG'+primer_end], \
	 'CACGAT':['NextFlexAdapter31',primer_start+'CACGAT'+primer_end], \
	 'CACTCA':['NextFlexAdapter32',primer_start+'CACTCA'+primer_end], \
	 'CAGGCG':['NextFlexAdapter33',primer_start+'CAGGCG'+primer_end], \
	 'CATGGC':['NextFlexAdapter34',primer_start+'CATGGC'+primer_end], \
	 'CATTTT':['NextFlexAdapter35',primer_start+'CATTTT'+primer_end], \
	 'CCAACA':['NextFlexAdapter36',primer_start+'CCAACA'+primer_end], \
	 'CGGAAT':['NextFlexAdapter37',primer_start+'CGGAAT'+primer_end], \
	 'CTAGCT':['NextFlexAdapter38',primer_start+'CTAGCT'+primer_end], \
	 'CTATAC':['NextFlexAdapter39',primer_start+'CTATAC'+primer_end], \
	 'CTCAGA':['NextFlexAdapter40',primer_start+'CTCAGA'+primer_end]}


def make_adapter_fa(index, index_dictionary, out_name, library_type):
	"""
	For each sample, a specific adapter fasta file is created based on its library type and index 
	This file is used by Trimmomatic to perform adapter trimming
	Index sequences were obtained from manufacturer files and are assigned as used at PCPGM
	"""	 
	outp = open(out_name, "w")
	if library_type in ["PE", "SPE"]:
		outp.write(">PrefixMultiplexingReadSequencingPrimer/1\nACACTCTTTCCCTACACGACGCTCTTCCGATCT\n")
		outp.write(">PrefixMultiplexingReadSequencingPrimer/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n")
	else:
		outp.write(">PrefixMultiplexingReadSequencingPrimer\nACACTCTTTCCCTACACGACGCTCTTCCGATCT\n")
	if index in index_dictionary:
		outp.write(">"+index_dictionary[index][0]+"\n"+index_dictionary[index][1]+"\n")
	outp.close()


def main(sample_info_file, standard_trim, path_start):
	"""
	Dispatches an lsf job to locate fastq files that were output by Casava (GIGPAD A1 routine) and then:
	1) Perform adapter trimming
	2) Run fastqc
	3) Get unique reads 
	4) Run tophat to align reads to reference genome
	5) Obtain various QC metrics on aligned files
	6) Run cufflinks to quantify ERCC spike ins (if applicable)
	7) Run cufflinks to quantify all other transcripts in sample
	
	The directory structure below is specific to the UPenn HPC cluster
	"""
	runs = get_sample_info(sample_info_file)
	#for curr_sample, k in runs.iteritems():
	for k in runs:
		#Get sample information
		curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type, lane, run = k

		#Get genome reference files
		ref_index, ERCC_index = get_genome_ref_files(ref_genome)

		#Set up project and sample output directories
		if path_start == "./":
			path_start = os.getcwd()
		if path_start[-1] != "/":
			path_start = path_start+"/"		
		project_dir = path_start+project+"/"
		if not os.path.exists(project_dir):
			os.makedirs(project_dir)
		out_dir = project_dir+curr_sample+"/"
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		
		job_name = curr_sample
				
		if lane != "-":
			#This directory structure and naming convention is (obviously!) unique to UPenn 
			#NGSC format: [RUN]_s_[LANE]_[1/2 corresponding to READ]_[BARCODE].fastq.gz
			if library_type in ("SE", "DGE"):
				R1 = project_dir+run+"_s_"+lane+"_1_"+index+".fastq.gz"
			else:
				R1 = project_dir+run+"_s_"+lane+"_1_"+index+".fastq.gz"
				R2 = project_dir+run+"_s_"+lane+"_2_"+index+".fastq.gz"		
		elif "SRR" in curr_sample:
			if library_type in ("SE", "DGE"):
				R1 = project_dir+curr_sample+".fastq.gz"
			else:
				R1 = project_dir+curr_sample+"_1.fastq.gz"
				R2 = project_dir+curr_sample+"_2.fastq.gz"			
		elif index == "-":
			#This directory structure and naming convention is (obviously!) unique to UPenn 
			if library_type in ("SE", "DGE"):
				R1 = project_dir+curr_sample+".fastq.gz"
			else:
				R1 = project_dir+curr_sample+"_1.fastq.gz"
				R2 = project_dir+curr_sample+"_2.fastq.gz"
		else:
			#This directory structure and naming convention is (obviously!) unique to UPenn 
			if library_type in ("SE", "DGE"):
				R1 = project_dir+curr_sample+"_"+index+".fastq.gz"
			else:
				R1 = project_dir+curr_sample+"_"+index+"_R1.fastq.gz"
				R2 = project_dir+curr_sample+"_"+index+"_R2.fastq.gz"
		
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
		
		#Check whether unaligned fastq files that were processed by Casava are present
		#Create directory for each sample
		if not os.path.isfile(R1):
			print "R1 file not found ", R1
		if library_type in ("PE", "SPE"):
			if not os.path.isfile(R2):		
				print "R2 file not found", R2

				
		outp.write("cd "+out_dir+"\n")
		R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
		R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"	
		
		#Perform adapter trimming with trimmomatic
		#May perform a standard trimming of bases from reads by amount given above if standard_trim variable is greater than 0. Most will use standard_trim=0 
		#Create fa file of adapters specific to file
		if not os.path.isfile(R1_trim):
			if library_type in ["PE", "SE", "SPE"]:
				make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
			elif library_type == "DGE":
				make_adapter_fa(index, nextflex_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
			if standard_trim == 0:
				R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
				R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"		
				if library_type in ["PE", "SPE"]:
					make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
				elif library_type in ["SE", "DGE"]:				
					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
			else:
				R1_trim = out_dir+curr_sample+"_R1_Trim"+str(standard_trim)+".fastq"
				R2_trim = out_dir+curr_sample+"_R2_Trim"+str(standard_trim)+".fastq"
				if library_type in ["PE", "SPE"]:
					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq HEADCROP:"+standard_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
				elif library_type in ["SE", "DGE"]:
					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" HEADCROP:"+standard_trim+"ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")					
		
		#Run fastqc on trimmed files.  
		#In some cases fastqc should be run on original files, but we have dropped this as a routine practice because the reports haven't changed much after trimming - adapter contamination has been minimal.
# 		if not os.path.isfile("*fastqc.zip"):						
# 			if library_type in ["PE", "SPE"]:
# 				outp.write("fastqc -o "+out_dir+" "+R1_trim+" "+R2_trim+"\n")
# 			elif library_type in ("SE", "DGE"):
# 				outp.write("fastqc -o "+out_dir+" "+R1_trim+"\n")
# 		
# 		#Get total number of reads, unique reads, % unique reads from trimmed file(s). 
# 		outp.write("cat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+curr_sample+"_ReadCount\n")
# 		if library_type in ["PE", "SPE"]:
# 			outp.write("cat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+curr_sample+"_ReadCount\n")
		
		#Run Kallisto
		if library_type in ["PE", "SPE"]:
			outp.write("kallisto quant -b 50 -i "+ref_index+" -o kallisto_output "+R1_trim+" "+R2_trim+"\n")
			outp.write("kallisto quant -b 50 -i "+ERCC_index+" -o kallisto_ercc "+R1_trim+" "+R2_trim+"\n")
		elif library_type in ("SE", "DGE"):
			print "Need to specify average fragment length"
			outp.write("kallisto quant --single -l 200 -s 20 -b 50 -i "+ref_index+" -o kallisto_output "+R1_trim+"\n")
			outp.write("kallisto quant --single -l 200 -s 20 -b 50 -i "+ERCC_index+" -o kallisto_ercc "+R1_trim+"\n")
		outp.close()
	
		#subprocess.call("bsub < "+job_name+".lsf", shell=True)
		#subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project.")
	parser.add_argument("--standard_trim", default=0, type=int, help="Number of bases to be trimmed from leftmost end of all reads (default=0)")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
	parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.samples_in, args.standard_trim, args.path_start)

