#!/usr/bin/python
import sys
import subprocess
import os
from rnaseq_align_and_qc import *


def read_RnaSeqMetrics(fin):
	"""
	Read in output file created by Picardtools RnaSeqMetrics function 
	Reformat metrics summary and histogram data into to be put into table/plot in Rnw report
	"""
	f = open(fin,'r')
	c = f.read().split('\n\n')[1:]
	metrics = c[0].split('\n')[1:]
	hist = c[1].split('\n')[1:]
	metrics_out, hist_out = [], []
	for x in metrics:
		metrics_out.append( x.split('\t') )
	for x in hist:
		hist_out.append( x.split('\t') )
	return metrics_out, hist_out


def read_InsertSizeMetrics(fin):
	"""
	Read in output file created by Picardtools CollectInsertSizeMetrics function 
	Reformat metrics summary and histogram data to be put into table/plot in Rnw report
	"""
	f = open(fin,'r')
	c = f.read().split('\n\n')[1:]
	metrics = c[0].split('\n')[1:]
	hist = c[1].split('\n')[1:]
	metrics_out, hist_out = [], []
	for x in metrics:
		metrics_out.append( x.split('\t') )
	for x in hist:
		hist_out.append( x.split('\t') )
	return metrics_out, hist_out
	

def read_samtools_stats(fin, ref_genome):
	"""
	Read in output file created by Samtools stats function (of type accepted_hits.sorted.stats)
	Reformat and output:
	1) ERCC raw counts to be put into table in Rnw report
	2) ref_genome summary data to be put into table and plot in Rnw report
	Note:	hg19 rRNA summary based on sum of chrUn_gl000220 and chrM
			mm10 rRNA not included. Only MT appears to have very high counts that could be from rRNA
			Zv9 rRNA not included.	
	"""
	f = open(fin,'r')
	c = f.read().split('\n')
	if '' in c:
		c.remove('')
	ercc_out, rna_out = {}, []
	rrna, other = 0, 0
	for x in c:
		if "ERCC" in x:
			curr_line = x.split('\t')[:-1]
			ercc_out[curr_line[0]] = curr_line[2]
		else:
			if ref_genome == "hg19":
				if ("chrM" in x) or ("chrUn_gl000220" in x):
					rrna += int(x.split('\t')[2])
				elif len(x.split('\t')[0].split('_')) > 1:
					other += int(x.split('\t')[2])
				elif "*" not in x:
					curr_line = x.split('\t')[:-1]
					curr_line[0] = curr_line[0].strip("chr")
					rna_out.append( curr_line )
			if ref_genome == "mm10":
				if len(x.split('\t')[0].split('.')) > 1:
					other += int(x.split('\t')[2])
				elif len(x.split('\t')[0].split('_')) > 1:
					other += int(x.split('\t')[2])
				elif ".1" and "_" and "*" not in x:
					curr_line = x.split('\t')[:-1]
					curr_line[0] = curr_line[0]
					rna_out.append( curr_line )
				rrna = "NA"
	rna_out.append(['Other','',str(other)])
	rna_out.append(['rRNA','', str(rrna)])
	return ercc_out, rna_out


def read_ercc_cufflinks(fin):
	"""
	Read in Cufflinks output file for ERCC spike-in FPKM values and output
	FPKM values to be put into table and plot in Rnw report	
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	ercc_out = {}
	for x in c:
		curr_line = x.split('\t')
		ercc_out[curr_line[0]] = curr_line[9]
	return ercc_out


def read_tophat_logs(fin, library_type):
	"""
	Read in TopHat output files of type (prep_reads.info) to get number of raw reads passed to program and number that were aligned
	Returns output for Rnw report
	"""
	f = open(fin,'r')
	c = f.read().split('\n')
	if '' in c:
		c.remove('')
	#First entry is 'left_reads_in'; Second entry is 'left_reads_out'
	#For paired-end data, third entry is 'right_reads_in'; fourth entry is 'right_reads_out'
	if library_type in ["PE", "SPE"]:
		read_numbers = map(lambda x: x.split('=')[1], c[2:4]+c[6:])
	else:
		read_numbers = map(lambda x: x.split('=')[1], c[2:4])
	return read_numbers


def read_bamtools_stats(fin):
	"""
	Read in bamstats output file accepted_hits.sorted.bamstats
	Reformats and outputs portions of interest for Rnw report
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[5:]
	for x in range(c.count('')):
		c.remove('')
	bs = map(lambda x: x.split(':'), c)
	names = map(lambda x: x[0], bs)
	counts = map(lambda x: x[1].split('\t')[0].strip(' '), bs)
	bamstats = zip(names, counts)
	return bamstats


def get_unique_reads(fin, library_type):
	"""
	Read in file that has output from awk on number of reads, unique reads, and %unique reads from fastq files
	Output R1 number of reads, R1 unique reads, R1 %unique reads, R2 number of reads, R2 unique reads, and R2 %unique reads for Rnw report
	"""
	f = open(fin,'r')
	c = f.read().split('\n')
	if '' in c:
		c.remove('')
	#First row is R1. Second row is R2. 
	read_numbers = map(lambda x: x.split(' '), c)
	if library_type == "PE":
		read_numbers = read_numbers[0]+read_numbers[1]
	else:
		read_numbers = read_numbers[0]
	return read_numbers



def read_fastq_data(fin):
	"""
	Read fastqc_data.txt file from a FastQC report and extract the percentage of duplicates and histogram information
	"""
	f = open(fin,'r')
	c = f.read().split('#Total Duplicate Percentage')[1]
	c2 = c.split('>>END_MODULE')[0]
	c3 = c2.split('\n')
	if '' in c3:
		c3.remove('')
	duplicate_data =  [['Total Duplicate Percentage', c3[0].strip('\t')]]+map(lambda x: x.split('\t'), c3[2:])
	return duplicate_data	


def make_project_data_files(project_name, sample_names, sample_paths, path_out, ref_genome, library_type):
	"""
	Creates several text files to be loaded by R for Rnw report based on modified program outputs read in with preceding scripts
	These text files are matrices containing information for all samples in a project/batch 
	Currently, cycles through all samples multiple times, once to create each individual file type 
	Currently, there is no way to handle missing files. If an error is encountered the process will stop at that point and not complete
	Currently, assumes default naming convention of all programs used in rnaseq_align_and_qc.py
	"""
	
	#Read counts obtained from TopHat log files
	outp5 = open(path_out+project_name+"_read_counts.txt", "w")
	if library_type in ["PE", "SPE"]:
		name5 = ["Sample", "R1_Raw_Reads", "R1_Aligned_Reads", "R2_Raw_Reads", "R2_Aligned_Reads"]
	else:
		name5 = ["Sample", "Raw_Reads", "Aligned_Reads"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		total_reads = read_tophat_logs(curr_path+"tophat_out/prep_reads.info", library_type)
		if i == 0:
			e = [[curr_name] + total_reads]
		else:
			e.append( [curr_name] + total_reads )
	outp5.write("\t".join(name5)+"\n")
	outp5.write("\n".join(map("\t".join, e)))
	#To avoid warning in R about end of line not being present:
	outp5.write("\n")
	outp5.close()
	print "Created file "+path_out+project_name+"_read_counts.txt"

	#Unique read counts obtained by comprehensive count of fastq files
	outp9 = open(path_out+project_name+"_unique_counts.txt", "w")
	if library_type in ["PE", "SPE"]:
		name9 = ["Sample", "R1_Raw_Reads", "R1_Unique_Reads", "R1_Percent_Unique", "R2_Raw_Reads", "R2_Unique_Reads", "R2_Percent_Unique"]
	else:
		name9 = ["Sample", "Raw_Reads", "Unique_Reads", "Percent_Unique"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		unique_reads = get_unique_reads(curr_path+curr_name+"_ReadCount", library_type)
		if i == 0:
			unique_counts = [[curr_name] + unique_reads]
		else:
			unique_counts.append( [curr_name] + unique_reads)
	outp9.write("\t".join(name9)+"\n")
	outp9.write("\n".join(map("\t".join, unique_counts)))
	outp9.close()
	print "Created file "+path_out+project_name+"_unique_counts.txt"
	
	#Duplicate read info from fastq files
	outp10 = open(path_out+project_name+"_duplicates.txt", "w")
	name10 = ["Read_Number"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		#Move individual FastQC reports to report directory
		subprocess.call("cp "+curr_path+"*fastqc.zip "+path_out, shell=True)
		name10.append(curr_name+"_R1")		 
		if os.path.isfile(curr_path+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"):
			fastqc1 = curr_path+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"
		elif os.path.isfile(curr_path+"/"+curr_name+"_R1_fastqc/fastqc_data.txt"):
			fastqc1 = curr_path+"/"+curr_name+"_R1_fastqc/fastqc_data.txt"
		elif os.path.isfile(curr_path+"/"+curr_name+"_fastqc/fastqc_data.txt"):
			fastqc1 = curr_path+"/"+curr_name+"_fastqc/fastqc_data.txt"
		else:
			print "Missing FastQC report", curr_name
			break
		curr_dup_R1 = read_fastq_data(fastqc1)
		if library_type in ["PE", "SPE"]:
			name10.append(curr_name+"_R2")
			if os.path.isfile(curr_path+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"):
				fastqc2 = curr_path+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"
			elif os.path.isfile(curr_path+"/"+curr_name+"_R2_fastqc/fastqc_data.txt"):
				fastqc2 = curr_path+"/"+curr_name+"_R2_fastqc/fastqc_data.txt"
			else:
				print "Missing FastQC report for R2 ", curr_name
			curr_dup_R2 = read_fastq_data(fastqc2)
		if i == 0:
			duplicates = curr_dup_R1
			if library_type in ["PE", "SPE"]:
				for j in range(len(duplicates)):
					duplicates[j].append(curr_dup_R2[j][1])
		else:
			for j in range(len(duplicates)):
				duplicates[j].append(curr_dup_R1[j][1])
				if library_type in ["PE", "SPE"]:
					duplicates[j].append(curr_dup_R2[j][1])				
	outp10.write("\t".join(name10)+"\n")
	outp10.write("\n".join(map("\t".join, duplicates)))
	outp10.close()
	print "Created file "+path_out+project_name+"_duplicates.txt"
	
	#rnaseqmetrics output split into two files: 1) overall metrics (2) data to create normalized coverage histogram for all samples
	outp1 = open(path_out+project_name+"_rnaseqmetrics_summary.txt", "w")
	name1 = ["Type"]
	outp2 = open(path_out+project_name+"_rnaseqmetrics_hist.txt", "w")
	name2 = ["Normalized_Position"]	
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		rnaseqmetrics_summary, rnaseqmetrics_hist = read_RnaSeqMetrics(curr_path+"tophat_out/"+curr_name+"_RNASeqMetrics")
		name1.append(curr_name)
		name2.append(curr_name)
		if i == 0:
			a = zip(rnaseqmetrics_summary[0], rnaseqmetrics_summary[1])
			b = rnaseqmetrics_hist
		else:
			for j in range(len(a)):
				a[j] = list(a[j])+[rnaseqmetrics_summary[1][j]]
			for j in range(len(b)):
				b[j] = b[j]+[rnaseqmetrics_hist[j][1]]
	outp1.write("\t".join(name1)+"\n")
	outp1.write("\n".join(map("\t".join, a) ))
	outp1.close()
	outp2.write("\t".join(name2)+"\n")
	outp2.write("\n".join(map("\t".join, b[1:])))
	outp2.close()
	print "Created files "+path_out+project_name+"_rnaseqmetrics_summary.txt and "+path_out+project_name+"_rnaseqmetrics_hist.txt"
	
	
	#samtools stats counts of reads per chromosome
	#ercc transcript read counts - combo of samtools stats output and cufflinks run with ERCC gtf file
	outp3 = open(path_out+project_name+"_counts.txt", "w")
	name3 = ["Chromosome", "Length"]	
	outp4 = open(path_out+project_name+"_ERCC.txt", "w")
	name4 = ["SpikeIn"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]	
		name3.append(curr_name)
		name4.append(curr_name+"_raw")
		name4.append(curr_name+"_fpkm")
		ercc_raw, rna_out = read_samtools_stats(curr_path+"tophat_out/"+curr_name+"_accepted_hits.sorted.stats", ref_genome)
		ercc_fpkm = read_ercc_cufflinks(curr_path+"cufflinks_out_ERCC/genes.fpkm_tracking")
		ercc_full ={}
		for k,v in ercc_raw.iteritems():
			ercc_full[k] = [v, ercc_fpkm[k]]
		if i == 0:
			c = rna_out 
			d = ercc_full
		else:
			for j in range(len(c)):
				c[j] = c[j]+[rna_out[j][2]]
			for k,v in d.iteritems():
				d[k] = d[k] + ercc_full[k]
	outp3.write("\t".join(name3)+"\n")
	outp3.write("\n".join(map("\t".join, c)))
	outp3.close()	
	outp4.write("\t".join(name4)+"\n")
	outp4.write("\n".join( map("\t".join, zip(d.keys(), map("\t".join, d.values()))) ))
	outp4.close()
	print "Created files "+path_out+project_name+"_counts.txt and "+path_out+project_name+"_ERCC.txt"
	
	#bamstats output metrics on types of reads, including junction spanning reads
	outp6 = open(path_out+project_name+"_bamstats_counts.txt", "w")
	name6 = ["Type"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		name6.append(curr_name)
		bam_stats = read_bamtools_stats(curr_path+"tophat_out/"+curr_name+"_accepted_hits.sorted.bamstats")
		if i == 0:
			f = bam_stats
		else:
			for j in range(len(bam_stats)):
				f[j] = list(f[j])+[bam_stats[j][1]]
	outp6.write("\t".join(name6)+"\n")
	outp6.write("\n".join(map("\t".join, f)))
	outp6.close()
	print "Created file "+path_out+project_name+"_bamstats_counts.txt"

	#insertsizemetrics output on insert size statistics
	if library_type in ["PE", "SPE"]:
		outp7 = open(path_out+project_name+"_insertmetrics_summary.txt", "w")
		name7 = ["Type"]
		for i in range(len(sample_names)):
			curr_name = sample_names[i]
			curr_path = sample_paths[i]
			name7.append(curr_name)
			#Make individual insert size files for each sample
			outp8 = open(path_out+project_name+"_"+curr_name+"_insertmetrics_hist.txt", "w")
			name8 = ["Insert_Size", curr_name]
			outp8.write("\t".join(name8)+"\n")
			insert_summary, insert_hist = read_InsertSizeMetrics(curr_path+"tophat_out/"+curr_name+"_InsertSizeMetrics")
			if i == 0:
				g = zip(insert_summary[0], insert_summary[1])
				h = insert_hist
			else:
				for j in range(len(g)):
					g[j] = list(g[j])+[insert_summary[1][j]]
				h = insert_hist
			outp8.write("\n".join(map("\t".join, h[1:])))
			outp8.close()
		outp7.write("\t".join(name7)+"\n")
		outp7.write("\n".join(map("\t".join, g) ))
		outp7.close()
		print "Created file "+path_out+project_name+"_insertmetrics_summary.txt and sample-specific *_insertmetrics_hist.txt files"	

def make_rnw_html(rnw_template, project_name, path_start, sample_names, ercc_mixes, ref_genome, library_type):
	"""
	Creates Rnw report. The top of report is below and the rest concatenated from a separate text document (rnw_template).
	Runs Sweave to create html document with R output (R2HTML library necessary for this)
	"""
	outp = open(path_start+project_name+"_QC_RnaSeqReport.Rnw", "w")
	outp.write("<HTML>\n")
	outp.write("<Head>\n")
	outp.write("<Title>\n")
	outp.write("RNA-seq Summary QC for Project "+project_name+"\n")
	outp.write("</Title>\n")
	outp.write("</Head>\n")
	outp.write("<Body>\n\n")
	outp.write("<<input.data,echo=FALSE,results=hide>>=\n")
	outp.write("library(RColorBrewer)\n")
	outp.write("curr.batch=\""+project_name+"\"\n")
	if "./" in path_start:
		outp.write("path.start=\""+path_start.lstrip("./")+"\"\n")
	else:
		outp.write("path.start=\""+path_start+"\"\n")
	outp.write("ercc.mixes=c("+str(ercc_mixes)[1:-1]+")\n")
	outp.write("sample.names=c("+str(sample_names)[1:-1]+")\n")
	outp.write("genome=\""+ref_genome+"\"\n")
	outp.write("library.type=\""+library_type+"\"\n")
	outp.write(rnw_template)
	outp.close()
	subprocess.call("cd "+path_start+"; echo \"library(R2HTML); Sweave('"+project_name+"_QC_RnaSeqReport.Rnw', driver=RweaveHTML)\" | R --no-save --no-restore", shell=True)


def main(project_name, sample_info_file, path_start):
	"""
	Creates html report describing summary and QC statistics for a set of aligned RNA-Seq samples associated with a project
	Report is based on multiple output files created by rnaseq_align_and_qc.py
	Such files are first reformatted into matrices that are easily loaded into R
	Report is created using Sweave
	Input:
		project_name: name for report. Could be a batch number or another name if several batches were used
		sample_info_file: tab delimited txt file with sample information as described in rnaseq_align_and_qc.py
		rnw_template: txt file that contains most of the contents that will populate the Rnw file for the report
	Current ref_genome choices: hg19, mm10
	Current library_type choices: PE, SE, DGE, SPE
	"""
	new_dir = path_start+project_name+"/"+project_name+"_Alignment_QC_Report/"
	if not os.path.exists(new_dir):
		os.makedirs(new_dir)

	#Get list of sample info. Fields: [customer_id, gigpad_id, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type]
	runs = get_sample_info(sample_info_file)
	sample_names = map(lambda x: x[0], runs)
	sample_paths = map(lambda x: path_start+x[6]+"/"+x[0]+"/", runs)
	ercc_mixes = map(lambda x: x[4], runs)
	ref_genome_list = map(lambda x: x[8], runs)
	library_type_list = map(lambda x: x[8], runs)
	#Check whether all samples are of same reference genome
	if False in map(lambda y: y==ref_genome_list[0], ref_genome_list):
		print "Make sure all samples in project are of the same reference genome"
		sys.exit()
	else:
		ref_genome = ref_genome_list[0]
	#Check whether all samples are of same reference genome
	if False in map(lambda y: y==library_type_list[0], library_type_list):
		print "Make sure all samples in project are of the same library type"
		sys.exit()
	else:
		library_type = library_type_list[0]
	
	#Create the data files needed for report
	make_project_data_files(project_name, sample_names, sample_paths, new_dir, ref_genome, library_type)
	
	#Create the report
	if not os.path.exists("rnaseq_align_and_qc_report_Rnw_template.txt"):
		print "Cannot find rnaseq_align_and_qc_report_Rnw_template.txt"
		sys.exit()
	rnw_in = open("rnaseq_align_and_qc_report_Rnw_template.txt", "r")
	rnw_template = rnw_in.read()
	make_rnw_html(rnw_template, project_name, new_dir, sample_names, ercc_mixes, ref_genome, library_type)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Create HTML report of QC and alignment summary statistics for RNA-seq samples associated with a project.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where PCPGM batch-level directories are located and report directory will be written (default=./)")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to. Often a PCPGM batch, but it could correspond to a mixture of batches.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.path_start)
	
	
