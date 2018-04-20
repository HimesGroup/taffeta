#!/usr/bin/python
#
#Created by Blanca Himes to simultaneously download SRA files for a study from the command line

#04/18/18 - note that the instructions for downloading the "SRA Sample" list have changed from the previous iteration of this script - old instructions no longer work
#1) go to http://www.ncbi.nlm.nih.gov/Traces/sra/  
#2) click on "Search" tab and make sure "SRA Objects" tab is selected
#3) enter SRP (project) number (e.g. SRP081599) and click "Search"
#4) click on the number of SRA experiments that come up 
#5) on the results page, click "Send to" --> "File" --> Format: "Accession List" --> "Create File." This will give you a text file containing all the SRR (run) numbers for that SRP (project)
#6) use this to make sample info file

#The full study can be downloaded by selecting "Run Browser" tab, then "Download" tab for the object that is the full study.
#


#11/3/14 Script works but no good way to check that final files were downloaded properly. There should be a checksum somewhere on the SRA site.

#Structure of ftp address from NIH SRA:
#ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP005/SRP005411/
#
#Structure of ftp address for a study:
#ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP033/SRP033351/
#ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP033/SRP033351/SRR1039511/SRR1039511.sra

import subprocess
import os
import argparse
from rnaseq_align_and_qc import *


ftp_root = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/"
ascp_root = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"

def get_sra(accession_num, fin, transfer_type):
	f = open(fin,'r')
	c = f.read().split('\n')
	#Create file of sra sizes in bytes
	f_out_name = accession_num+"_sizes.txt"
	f_out = open(f_out_name, 'w')
	f_out.close()
	if '' in c:
		c.remove('')
	for i in c:
		out_f = "/Users/bhimes/Desktop/SRP043162/"+i+".sra"
		if transfer_type=="ascp":
			#use ascp for larger files but cannot get install to work properly. Holding this in case future install works. (Cannot locate openssh file and getting license error).
			ascp_a = ascp_root+i[0:3]+"/"+i[0:6]+"/"+i+"/"+i+".sra"
		else:
			#curl works very slowly, especially for larger dataset groups
			ftp_a = ftp_root+i[0:3]+"/"+i[0:6]+"/"+i+"/"+i+".sra"
			#To get size of files to be transferred:
			#subprocess.call("curl -sI "+ftp_a+" | awk \'/[Cc]ontent-[Ll]ength/ { print \""+i+"\",$2 }\' >> "+f_out_name, shell=True)
			#To get files transferred:
			subprocess.call("curl -s -o \""+out_f+"\" \""+ftp_a+"\" &\n", shell=True)


#get_sra("/Users/bhimes/Dropbox/R01_Feb_2015/EpitheliumExpression/SRP005411/SRR_Acc_List.txt", "ftp")
#get_sra("SRP033351", "/Users/bhimes/Dropbox/SRP033351/SRP033351_Acc_List.txt", "ftp")
#get_sra("SRP043162", "/Users/bhimes/Desktop/SRP043162/SRP043162_Acc_List3.txt", "ftp")


#On consign cluster, there is no need to get sra files first. These are downloaded to directory where fastq-dump is called and then converted to fastq files.
def main(sample_info_file, path_start, project):
	#Set up project and sample output directories
	if path_start == "./":
		path_start = os.getcwd()
	if path_start[-1] != "/":
		path_start = path_start+"/"		
	runs = get_sample_info(sample_info_file)
	for k in runs:
		sample_id = k[0]
		print sample_id
		project_dir = path_start+project+"/"
		if not os.path.exists(project_dir):
			os.makedirs(project_dir)
		print "fastq-dump "+sample_id
		subprocess.call("/project/bhimeslab/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump.2.9.0 --split-3 --gzip --outdir "+project_dir+" "+sample_id+" &\n", shell=True)
		#MS note: I changed this to use the most recent version of the sratoolkit
		#the one currently available through HPC as a module is sratoolkit-2.4.2
		#however, using it produces the following error: "fastq-dump.2.4.2 err: name incorrect while evaluating path within network system module - Scheme is 'https'," which is fixed by using the newer version


		
#MS new 04/19/2018

#def fastq_download(fin):
#	f = open(fin,'r')
#	c = f.read().split('\n')
#	for x in c: 
#		print x
#		print "Now starting prefetch"
#		#subprocess.call("/project/bhimeslab/sratoolkit.2.9.0-centos_linux64/bin/prefetch.2.9.0 -v "+x+"", shell=True)
#	for x in c: 
#		print x
#		print "Now starting fastq-dump"
#		subprocess.call("/project/bhimeslab/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump.2.9.0 '/project/bhimeslab/ncbi/sra/' "+x+" --gzip --outdir 'GSE85568' &\n", shell=True)

#fastq_download("SraAccList.txt")


#get_fastq("/Users/bhimes/Dropbox/R01_Feb_2015/EpitheliumExpression/SRP005411/SRR_Acc_List.txt")
#get_fastq("/Users/bhimes/Dropbox/SRP033351/SRP033351_Acc_List.txt")
#get_fastq("/project/bhimeslab/SRP005411/SRP005411_Acc_List.txt")
#get_fastq("/project/bhimeslab/SRP043162/SRP043162_Acc_List.txt")
#get_fastq("/project/bhimeslab/SRP033351/SRP033351_Acc_List.txt")
#get_fastq("SRP033335_Info_Sheet.txt")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Download individual sample fastq files from the SRA.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
	parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information.")
	parser.add_argument("project", help="Directory name for this project. Corresponds to GEO Series (Subseries) name - e.g. GSE85567")
	args = parser.parse_args()
	main(args.samples_in, args.path_start, args.project)
	
	
	
