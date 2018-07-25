#!/usr/bin/python

# Libraries
import argparse
import os
import sys
import subprocess

def make_rmd_html(rmd_template, geo_id, project_name, path_start, out_dir, pheno_info, rmd_file):
    """
    Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
    Two files generated: 1) geo_id+"_withoutQC.txt" raw phenotype data from GEO (if phenotype_fn is not specified) 2) project_name+"_sraFile.info" list ftp address of sra files
    """
    outp = open(rmd_file, "w")
    outp.write("---\n")
    outp.write("title: " + project_name + "SRA download\n")
    outp.write("author: 'Mengyuan Kan (mengykan@upenn.edu)'\n")
    outp.write("date: \"`r format(Sys.time(), '%d %B, %Y')`\"\n")
    outp.write("output:\n")
    outp.write("  html_document:\n")
    outp.write("    toc: TRUE\n")
    outp.write("    depth: 3\n")
    outp.write("editor_options:\n")
    outp.write("chunk_output_type: console\n")
    outp.write("---\n\n")

    outp.write("Assign the variables for GEO ID (geo_id), data directory (out_dir), phenotype file if user defined (pheno_fn).\n\n")

    outp.write("```{r var, echo=T}\n")
    outp.write("out_dir <- '" + out_dir + "'\n")
    outp.write("project_name <- '" + project_name + "'\n")
    if pheno_info is not None:
        outp.write("pheno_info <- '" + pheno_info + "'\n")
    outp.write("geo_id <- '" + geo_id + "'\n")
    outp.write("```\n\n")

    outp.write("\n")

    # insert contents in rmd_template
    outp.writelines(rmd_template)
    outp.write("\n")

    outp.close()

def lsf_file(job_name, cmd, memory=36000):
    """
    Creates .lsf files
    """

    outp = open(job_name+".lsf",'w')
    outp.write("#!/bin/bash\n")
    outp.write("#BSUB -L /bin/bash\n")
    outp.write("#BSUB -J "+job_name+"\n")
    outp.write("#BSUB -q normal\n")
    outp.write("#BSUB -o "+job_name+"_%J.out\n")
    outp.write("#BSUB -e "+job_name+"_%J.screen\n")
    outp.write("#BSUB -M "+str(memory)+"\n")
    outp.write("#BSUB -n 1\n")
    outp.write(cmd)
    outp.write("\n")
    outp.close()


def download_fastq(SRA_info, out_dir, path_start, fastqc):
    """
    Obtain ftp download link for samples in project_name + "_sraFile.info" from SRA
    Creates .lsf files for .fastq file download
    """

    sra=open(SRA_info, 'r')
    sra=sra.readlines()
    header=sra[0].rstrip()
    if header != "run\tsubmission\tstudy\tsample\texperiment\tftp":
        print "The column names in SRA_info file "+SRA_info+" should be run, submission, study, sample, experiment, ftp\n"
        sys.exit()

    SRA={} # key: Sample ID, value: ftp
    for line in sra[1:]:
        line = line.rstrip().split('\t')
        run = line[0]
        ftp = line[5]
        if run not in SRA:
            SRA[run]=[ftp]
        else:
            SRA[run]=SRA[run]+[ftp]

    # create .lsf file with download command
    for sample in SRA:
        fastq_names=map(lambda x: x.split('/')[-1] ,SRA[sample]) # exclude ftp path and obtain fastq file name

        if fastqc: # create new directory under current sample name for raw fastqc results
            fastqc_out=path_start+sample
            if not os.path.exists(fastqc_out):
                os.makedirs(fastqc_out)

        for i in range(len(fastq_names)):
            fastq_name=fastq_names[i]
            ftp=SRA[sample][i]
            if os.path.exists(out_dir+fastq_name):
                print(out_dir+fastq_name+" already exists. Skip download.") # if the fastq file already exists, skip download
                if fastqc:
                    download_cmd="fastqc "+out_dir+fastq_name+" -o "+fastqc_out+"\n"
                    lsf_file(sample+"_"+str(i+1), download_cmd) # create .lsf files for fastqc
            else:
                download_cmd="cd "+out_dir+"\n"
                download_cmd=download_cmd+"wget "+ ftp+"\n"
                if fastqc:
                    download_cmd=download_cmd+"fastqc "+out_dir+fastq_name+" -o "+fastqc_out+"\n"

                lsf_file(sample+"_"+str(i+1), download_cmd) # create .lsf files for download and/or fastqc
            

def main(geo_id, path_start, project_name, pheno_info, template_dir, fastqc):

    # Set up project and sample output directories
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
	path_start = path_start+"/"		

    out_dir=path_start+project_name+"_SRAdownload/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print("Create directory " + out_dir)

    if template_dir == "./":
        template_dir = os.getcwd()
    if template_dir[-1] != "/":
        template_dir = template_dir+"/"

    # check if QC template txt file exists
    rmd_template_fn=template_dir+"rnaseq_sra_download_Rmd_template.txt"
    if not os.path.exists(rmd_template_fn):
        print "Cannot find "+ rmd_template_fn
	sys.exit()

    # read in template file
    rmd_in = open(rmd_template_fn, "r")
    rmd_template = rmd_in.read()

    # create and run rmd report file to obtain sample infomation file from GEO and SRA
    rmd_file = out_dir + project_name + "_SRAdownload_RnaSeqReport.Rmd"
    make_rmd_html(rmd_template, geo_id, project_name, path_start, out_dir, pheno_info, rmd_file)
    ftpinfo_cmd="echo \"library(rmarkdown); rmarkdown::render('"+rmd_file+"')\" | R --no-save --no-restore"
    subprocess.call(ftpinfo_cmd, shell=True)

    # check if RSA info file is generated
    SRA_info = out_dir + project_name + "_sraFile.info"
    if not os.path.exists(SRA_info):
        print "The SRA file " + SRA_info +"does not exist. Please check!"
        sys.exit()

    # download .fastq file based on ftp address in generated project_name+sraFile.info file
    # assign TURE to fastqc variable to do fastqc for raw .fastq files. Fastqc results are saved in a new directory under the current sample's name
    download_fastq(SRA_info, out_dir, path_start, fastqc)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download RNA-Seq rawr reads .fastq files from SRA.")
    parser.add_argument("--geo_id", type=str, help="GEO accession id")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
    parser.add_argument("--project_name", type=str, help="Name of project that all samples correspond to.")
    parser.add_argument("--pheno_info", help="Use user defined phenotype file to download the .fastq files of corresponding samples from SRA. The name in 'SRA_ID' column should match the sample id in SRA database. If not defined, download samples based on GEO phenotype field 'relation.1'.")
    parser.add_argument("--template_dir", default="./", type=str, help="The template RMD file rnaseq_sra_download_Rmd_template.txt is located.")
    parser.add_argument("--fastqc", action='store_true', help="If specified, run fastqc for downloaded raw .fastq files.")
    args = parser.parse_args()

    if args.geo_id is None or args.project_name is None:
        parser.print_help()
        sys.exit()

    main(args.geo_id, args.path_start, args.project_name, args.pheno_info, args.template_dir, args.fastqc)
