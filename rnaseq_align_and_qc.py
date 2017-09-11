#!/usr/bin/python
import subprocess
import os
import argparse
import glob


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
	v7: library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE", "SSE"
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end", "stranded single-end")
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


def get_genome_ref_files(genome, mask):
	"""
	Location of all reference files needed for a given genome.
	The ERCC gtf files were appended separately to each species own gtf file
	Current choice: "hg38"
	"""
	mask_gtf = "NA"
	if genome == "hg38":
		ref_index = "/project/bhimeslab/Reference/hg38/genome.ERCC"
		fa = "/project/bhimeslab/Reference/hg38/genome.ERCC.fa"
		gtf = "/project/bhimeslab/Reference/hg38/genes.gtf"
		ref = "/project/bhimeslab/Reference/hg38/refFlat.txt"
		ERCC_gtf = "/project/bhimeslab/Reference/hg38/genes.ERCC.gtf"
		genome_dir = "/project/bhimeslab/Reference/hg38"
		if mask == "rRNA":
			mask_gtf = "/project/bhimeslab/Reference/hg38/rRNA_hg38.gtf"
		elif mask == "globin":
			mask_gtf = "/project/bhimeslab/Reference/hg38/globin_hg38.gtf"
		elif mask == "rRNA_globin":
		 	mask_gtf = "/project/bhimeslab/Reference/hg38/rRNA_globin_hg38.gtf"
		else:
			print 'Unknown mask selected: ', mask
	elif genome == "hg19":
		ref_index = "/project/bhimeslab/Reference/hg19/genome.ERCC"
		fa = "/project/bhimeslab/Reference/hg19/genome.ERCC.fa"
		gtf = "/project/bhimeslab/Reference/hg19/genes.gtf"
		ref = "/project/bhimeslab/Reference/hg19/refFlat.txt"
		ERCC_gtf = "/project/bhimeslab/Reference/hg19/genes.ERCC.gtf"
		genome_dir = "/project/bhimeslab/Reference/hg19"
		print 'No mask files available. Do not have STAR indexes for hg19.'
	elif genome == "mm38":
		ref_index = "/project/bhimeslab/Reference/mm38/mm38_ERCC"
		fa = "/project/bhimeslab/Reference/mm38/mm.GRCm38.genome.ERCC92.fa"
		gtf = "/project/bhimeslab/Reference/mm38/mm.GRCm38.75.genes.gtf"
		ref = "/project/bhimeslab/Reference/mm38/refFlat.txt"
		ERCC_gtf = "/project/bhimeslab/Reference/mm38/mm.GRCm38.75.genes.ERCC.gtf"
		genome_dir = "/project/bhimeslab/Reference/mm38"
		print 'No mask files available. Do not have STAR indexes for mm38.'
	elif genome == "mm10":
		ref_index = "/project/bhimeslab/Reference/mm38/UCSC_mm10_ERCC"
		fa = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genome.ERCC92.fa"
		gtf = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genes.gtf"
		ref = "/project/bhimeslab/Reference/mm38/refFlat_UCSC.txt"
		ERCC_gtf = "/project/bhimeslab/Reference/mm38/UCSC_mm10_genes.ERCC.gtf"
		genome_dir = "/project/bhimeslab/Reference/mm38"
		print 'No mask files available.'
	elif genome == "rn6":
		ref_index = "/project/bhimeslab/Reference/rn6/genome.ERCC"
		fa = "/project/bhimeslab/Reference/rn6/genome.ERCC.fa"
		gtf = "/project/bhimeslab/Reference/rn6/genes.gtf"
		ref = "/project/bhimeslab/Reference/rn6/refFlat.txt"
		ERCC_gtf = "/project/bhimeslab/Reference/rn6/genes.ERCC.gtf"
		genome_dir = "/project/bhimeslab/Reference/rn6"
		print 'No mask files available. Do not have STAR indexes for rn6.'
	elif genome == "susScr3":
		ref_index = "/project/bhimeslab/Reference/susScr3/genome.ERCC"
		fa = "/project/bhimeslab/Reference/susScr3/Ensembl_susScr3_genome.ERCC92.fa"
		gtf = "/project/bhimeslab/Reference/susScr3/genes.gtf"
		ref = "/project/bhimeslab/Reference/susScr3/refFlat.txt"
		ERCC_gtf = "/project/bhimeslab/Reference/susScr3/Ensembl_susScr3_genes.ERCC.gtf"
		genome_dir = "/project/bhimeslab/Reference/susScr3"
		print 'No mask files available. Do not have STAR indexes for susScr3.'
	else:
		print 'Unknown genome selected: ', genome
	return(ref_index, fa, gtf, ref, ERCC_gtf, mask_gtf, genome_dir)

		
#ERCC gtf file
ERCC_only = "/project/bhimeslab/Reference/ERCC92.gtf"	
#hg38 rRNA file
rRNA_gtf = "/project/bhimeslab/Reference/hg38/rRNA_hg38.gtf"

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
	else:
		print "Assuming the index provided works with standard Illumina primer sequences"
		outp.write(">"+index+"\n"+primer_start+index+primer_end+"\n")
	outp.close()


def main(sample_info_file, discovery, standard_trim, mask, aligner, path_start):
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
		ref_index, fa, gtf, ref, ERCC_gtf, mask_gtf, genome_dir = get_genome_ref_files(ref_genome, mask)

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
			#[RUN]_s_[LANE]_[1/2 corresponding to READ]_[BARCODE].fastq.gz
			if library_type in ("SE", "DGE"):
				R1 = project_dir+run+"_s_"+lane+"_1_"+index+".fastq.gz"
			#for mir34a specifically
			elif library_type in ("SSE"): 
				R1=project_dir+curr_sample+"_"+index+".R1.fastq.gz" #changed from '_R1' to '.R1' for hbe_enterolactone
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
				#For rat_stretch, the index is not part of file name
				#R1 = project_dir+curr_sample+".fastq.gz"
				R1 = project_dir+curr_sample+"_"+index+".R1.fastq.gz"
			else:
				#For pig_lung, the index is not part of file name
				R1 = project_dir+curr_sample+"_1.fastq.gz"
				R2 = project_dir+curr_sample+"_2.fastq.gz"
				#R1 = project_dir+curr_sample+"_"+index+"_1.fastq.gz"
				#R2 = project_dir+curr_sample+"_"+index+"_2.fastq.gz"
		
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
		
		#In case need to unzip at some point:
		#outp.write("zcat "+R1+".gz > "+R1[:-3]+"\n")
		#outp.write("zcat "+R2+".gz > "+R2[:-3]+"\n")
				
		outp.write("cd "+out_dir+"\n")
 		R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
 		R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"	
 		#outp.write("zcat "+R1+" > "+ R1_trim+"\n")
 		#R1_trim = R1
 		#R2_trim = R2
 		#print "R1: ", R1_trim	


		#Perform adapter trimming with trimmomatic
		#May perform a standard trimming of bases from reads by amount given above if standard_trim variable is greater than 0. Most will use standard_trim=0 
		#Create fa file of adapters specific to file
 		if not os.path.isfile(R1_trim):
 			if library_type in ["PE", "SE", "SPE", "SSE"]:
 				make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
 			elif library_type == "DGE":
 				make_adapter_fa(index, nextflex_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
 			if standard_trim == 0:
 				R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
 				R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"		
 				if library_type in ["PE", "SPE"]:
 					make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
 					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
 				elif library_type in ["SE", "DGE", "SSE"]:				
 					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
 			else:
 				R1_trim = out_dir+curr_sample+"_R1_Trim"+str(standard_trim)+".fastq"
 				R2_trim = out_dir+curr_sample+"_R2_Trim"+str(standard_trim)+".fastq"
 				if library_type in ["PE", "SPE"]:
 					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq HEADCROP:"+standard_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")
 				elif library_type in ["SE", "DGE", "SSE"]:
 					outp.write("java -Xmx1024m  -classpath /opt/software/Trimmomatic/0.32/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" HEADCROP:"+standard_trim+"ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n")					

		#Run fastqc on trimmed files. 
		if not glob.glob(out_dir+curr_sample+'*fastqc.zip'):
			#In some cases fastqc should be run on original files, but we have dropped this as a routine practice because the reports haven't changed much after trimming - adapter contamination has been minimal.
			if library_type in ["PE", "SPE"]:
				outp.write("fastqc -o "+out_dir+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type in ("SE", "DGE", "SSE"):
				outp.write("fastqc -o "+out_dir+" "+R1_trim+"\n")

		#Get total number of reads, unique reads, % unique reads from trimmed file(s). 
		outp.write("cat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+curr_sample+"_ReadCount\n")
		if library_type in ["PE", "SPE"]:
			outp.write("cat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+curr_sample+"_ReadCount\n")
		
 		if aligner == "tophat":
 			#Run TopHat with options specific to library type
			if discovery == "no":
				if library_type == "PE":
					outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
				elif library_type == "SPE":
					outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
				elif library_type == "SSE":
					outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
				elif library_type in ["DGE", "SE"]:
					outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
		
			elif discovery == "yes":
				if library_type == "PE":
					outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
				elif library_type == "SPE":
					outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
				elif library_type in ["DGE", "SE"]:
					outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
			outp.write("cd "+out_dir+"/tophat_out/\n")
				
		elif aligner == "star":
			outp.write("mkdir star_out\n")
			outp.write("cd star_out\n")
			if library_type == "SPE":
				outp.write("/project/bhimeslab/STAR-2.5.2b/bin/Linux_x86_64/STAR --genomeDir "+genome_dir+" --runThreadN 12 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim+" "+R2_trim+"\n")
			elif library_type == "PE":
				outp.write("/project/bhimeslab/STAR-2.5.2b/bin/Linux_x86_64/STAR --genomeDir "+genome_dir+" --runThreadN 12 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim+" "+R2_trim+"\n")			
			elif library_type in ["DGE", "SE"]:
				outp.write("/project/bhimeslab/STAR-2.5.2b/bin/Linux_x86_64/STAR --genomeDir "+genome_dir+" --runThreadN 12 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim+"\n")
			elif library_type in ["SSE"]:
				outp.write("/project/bhimeslab/STAR-2.5.2b/bin/Linux_x86_64/STAR --genomeDir "+genome_dir+" --runThreadN 12 --sjdbOverhang 100 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim+"\n")	
			outp.write("mv Aligned.sortedByCoord.out.bam accepted_hits.bam\n")
			outp.write("mkdir "+out_dir+"/htseq_out/\n")
			outp.write("samtools view accepted_hits.bam | htseq-count -r pos - "+gtf+" > "+out_dir+"/htseq_out/"+curr_sample+"_counts.txt\n")
		
		#Get samtools mapping stats
		#Create sorted bam file:
		outp.write("samtools sort accepted_hits.bam "+curr_sample+"_accepted_hits.sorted\n")
		#Create indexed bam file:
		outp.write("samtools index "+curr_sample+"_accepted_hits.sorted.bam\n")
		#Write out index stats of where reads align to by chr:
		outp.write("samtools idxstats "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.stats\n")
		#Write out bamtools summary stats:
		outp.write("bamtools stats -in "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.bamstats\n")
		#Run CollectRnaSeqMetrics
		if library_type == "SPE":
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND VALIDATION_STRINGENCY=LENIENT INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		elif library_type == "SSE":
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND VALIDATION_STRINGENCY=LENIENT INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	 
		else:
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=NONE VALIDATION_STRINGENCY=LENIENT INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		
		#Get number of reads spanning junctions by getting "N"s in CIGAR field of bam file
		#Be sure that Junction Spanning Reads are Added First then Unmapped Reads for proper ordering of fields in report
		outp.write("echo \"Junction Spanning Reads: \" $(bamtools filter -in "+curr_sample+"_accepted_hits.sorted.bam -script /project/bhimeslab/taffeta/cigarN.script | bamtools count ) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")
		
		if aligner == "tophat":
			#Get number of unmapped reads
			outp.write("echo Unmapped Reads: $(samtools view -c unmapped.bam) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")		

		#Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
		if library_type in ["PE", "SPE"]:
			outp.write("java -Xmx2g -jar /opt/software/picard/picard-tools-1.96/CollectInsertSizeMetrics.jar VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE="+curr_sample+"_InsertSizeHist.pdf INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_InsertSizeMetrics\n")	
		
		#Run cufflinks to count ERCC spike ins
		outp.write("mkdir "+out_dir+"/cufflinks_out_ERCC/\n")
		outp.write("cd "+out_dir+"/cufflinks_out_ERCC/\n")
		if library_type in ["DGE"]:
			outp.write("cufflinks --library-type fr-unstranded --no-length-correction -G "+ERCC_only+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
		elif library_type in ["SPE"]:
			outp.write("cufflinks --library-type fr-firststrand -G "+ERCC_only+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
		elif library_type in ["SSE"]:
			outp.write("cufflinks --library-type fr-secondstrand -G "+ERCC_only+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
		else:
			outp.write("cufflinks --library-type fr-unstranded -G "+ERCC_only+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")

		#Run cufflinks to count rRNA. Currently only works with hg38
		if ref_genome == "hg38":
			outp.write("mkdir "+out_dir+"/cufflinks_out_rRNA/\n")
			outp.write("cd "+out_dir+"/cufflinks_out_rRNA/\n")
			if library_type == "DGE":
				outp.write("cufflinks --library-type fr-unstranded --no-length-correction -G "+rRNA_gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
			elif library_type in ["SPE"]:
				outp.write("cufflinks --library-type fr-firststrand -G "+rRNA_gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
			elif library_type in ["SSE"]:
				outp.write("cufflinks --library-type fr-secondstrand -G "+rRNA_gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
			else:
				outp.write("cufflinks --library-type fr-unstranded -G "+rRNA_gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
		
		#Cufflinks to assemble and quantify transcripts
		outp.write("mkdir "+out_dir+"/cufflinks_out/\n")
		outp.write("cd "+out_dir+"/cufflinks_out/\n")
		if mask == "none":
			if library_type in ["DGE"]:
				outp.write("cufflinks --library-type fr-unstranded --no-length-correction -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
			elif library_type in ["SPE"]:
				outp.write("cufflinks --library-type fr-firststrand -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")			
			elif library_type in ["SSE"]:
				outp.write("cufflinks --library-type fr-secondstrand -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")			
			else:
				outp.write("cufflinks --library-type fr-unstranded -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
		else:
			if library_type in ["DGE"]:
				outp.write("cufflinks --library-type fr-unstranded --no-length-correction -M "+mask_gtf+" -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
			elif library_type in ["SPE"]:
				outp.write("cufflinks --library-type fr-firststrand -M "+mask_gtf+" -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")			
			elif library_type in ["SSE"]:
				outp.write("cufflinks --library-type fr-secondstrand -M "+mask_gtf+" -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")			
			else:
				outp.write("cufflinks --library-type fr-unstranded -M "+mask_gtf+" -G "+gtf+" -p 12 ../"+aligner+"_out/"+curr_sample+"_accepted_hits.sorted.bam \n")
		#outp.write("rm ../"+aligner+"_out/accepted_hits.bam \n")	#commented out for now
		outp.close()
	
		#subprocess.call("bsub < "+job_name+".lsf", shell=True)
		#subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project.")
	parser.add_argument("--standard_trim", default=0, type=int, help="Number of bases to be trimmed from leftmost end of all reads (default=0)")
	parser.add_argument("--discovery", default="yes", type=str, help="Should TopHat be run with options to discover novel transcripts (i.e. disable --no-novel-juncs --transcriptome-only)? "
		"(options: yes, no; default=yes) ")
	parser.add_argument("--mask", default="none", type=str, help="Should Cufflinks be run with options to mask transcripts (i.e. -M mask.gtf)? "
		"(options: rRNA, globin, rRNA_globin, none; default=none) ")
	parser.add_argument("--aligner", default="tophat", type=str, help="Should TopHat or STAR be used as aligner?"
		"(options: tophat, star)")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
	parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.samples_in, args.discovery, args.standard_trim, args.mask, args.aligner, args.path_start)

