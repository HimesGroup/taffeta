#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import re
import rnaseq_userdefine_variables as userdef 

def get_sample_info(fin):
    """
    Read in information from phenotype file
    Create a directory to store each column
    """
    f = open(fin,'r')
    f = f.read().split('\n')
    f = map(lambda x: x.rstrip(), f)

    if '' in f:
        f.remove('')
    header = f[0].split('\t') # list
    c = f[1:]
    
    # Obtain column index from file header
    if "Sample" not in header:
        print "Sample column is missing. Please check!"
        sys.exit()

    d = {} # key: column name, value: column

    for i in range(len(header)):
        colname=header[i]
        d[colname]=map(lambda x: x.split('\t')[i],c)

    return d


###
# Read in and process RnaSeqMetrics files
###

def read_RnaSeqMetrics(fin):
    """
    Read in output file created by Picardtools RnaSeqMetrics function 
    Reformat metrics summary and histogram data into to be put into table/plot in Rmd report
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


def make_RnaSeqMetrics_matrix(path_out, project_name, sample_names, sample_paths, aligner):

	outp1 = open(path_out+project_name+"_rnaseqmetrics_summary.txt", "w")
	name1 = ["Type"]
	outp2 = open(path_out+project_name+"_rnaseqmetrics_hist.txt", "w")
	name2 = ["Normalized_Position"]	
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		rnaseqmetrics_summary, rnaseqmetrics_hist = read_RnaSeqMetrics(curr_path+aligner+"_out/"+curr_name+"_RNASeqMetrics")
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


###
# Read in and process bamstat files
###

def read_InsertSizeMetrics(fin):
	"""
	Read in output file created by Picardtools CollectInsertSizeMetrics function 
	Reformat metrics summary and histogram data to be put into table/plot in Rmd report
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
	

def make_insertsize_matrix(path_out, project_name, sample_names, sample_paths, library_type, aligner):

	if library_type in ["PE"]:
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
			insert_summary, insert_hist = read_InsertSizeMetrics(curr_path+"/"+aligner+"_out/"+curr_name+"_InsertSizeMetrics")
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

###
# Read in and process samtools idxstat files
###


def read_samtools_stats(fin, ref_genome):
    """
    Read in output file created by Samtools stats function (of type accepted_hits.sorted.stats)
    Reformat and output:
    1) ERCC raw counts to be put into table in Rmd report
    2) ref_genome summary data to be put into table and plot in Rmd report
    Note: hg19 rRNA summary based on sum of chrUn_gl000220 and chrM
	  mm38 rRNA not included. Only MT appears to have very high counts that could be from rRNA
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
	    if ref_genome == "hg38":
		if ("chrM" in x) or ("chrUn_GL000220v1" in x):
		    rrna += int(x.split('\t')[2])
		elif len(x.split('\t')[0].split('_')) > 1:
		    other += int(x.split('\t')[2])
		elif ("*" not in x) and ("EBV" not in x):
		    curr_line = x.split('\t')[:-1]
		    curr_line[0] = curr_line[0].strip("chr")
		    rna_out.append( curr_line )
	    if ref_genome == "mm10" or ref_genome == "mm38":
		if len(x.split('\t')[0].split('.')) > 1:
		    other += int(x.split('\t')[2])
		elif len(x.split('\t')[0].split('_')) > 1:
		    other += int(x.split('\t')[2])
		elif ".1" and "_" and "*" not in x:
		    curr_line = x.split('\t')[:-1]
		    curr_line[0] = curr_line[0]
		    rna_out.append( curr_line )
		rrna = "NA"
	    if ref_genome == "rn6":
		if len(x.split('\t')[0].split('.')) > 1:
		    other += int(x.split('\t')[2])
		elif len(x.split('\t')[0].split('_')) > 1:
		    other += int(x.split('\t')[2])
		elif ".1" and "_" and "*" not in x:
		    curr_line = x.split('\t')[:-1]
		    curr_line[0] = curr_line[0]
		    rna_out.append( curr_line )
                rrna=NA
	    if ref_genome == "susScr3":
		if ("chrM" in x) or ("chrUn_gl000220" in x):
		    rrna += int(x.split('\t')[2])
		elif len(x.split('\t')[0].split('_')) > 1:
		    other += int(x.split('\t')[2])
		elif "*" not in x:
		    curr_line = x.split('\t')[:-1]
		    curr_line[0] = curr_line[0].strip("chr")
		    rna_out.append( curr_line )
	    if ref_genome == "Zv9":
		if len(x.split('\t')[0].split('_')) > 1:
		    other += int(x.split('\t')[2])
		elif "*" not in x:
		    curr_line = x.split('\t')[:-1]
		    curr_line[0] = curr_line[0]
		    rna_out.append( curr_line )
		rrna = "NA"
    rna_out.append(['Other','',str(other)])
    rna_out.append(['rRNA','', str(rrna)])
    return ercc_out, rna_out



def make_samidxstat_matrix(path_out, project_name, sample_names, sample_paths, ref_genome, aligner):
    """
    Read in individual output files created by curr_path+"/"+aligner+"_out/"+curr_name+"_accepted_hits.sorted.stats" according to sample_info_file
    Return as a single matrix for a whole batch
    """

    # path_out: directory for output files
    # project_name: file prefix
    # sample_names: a list stores all sample names
    # sample_paths: a list stores all directory paths for each sample
    # ref_genome: reference genome

    name1 = ["Chromosome", "Length"] # create file header
    c=[] # a list to store samtools idx stats for each sample by chromosome

    for i in range(len(sample_names)):
        curr_name = sample_names[i]
	curr_path = sample_paths[i]+aligner+"_out/"+curr_name+"_accepted_hits.sorted.stats"
	if not os.path.exists(curr_path):
            print "Missing samtools idxstat output file ", curr_path
	    break
 
        ercc_raw, rna_out = read_samtools_stats(curr_path, ref_genome)
        #print "Read in samtools idxstat output for sample "+curr_name+": Done"
        name1.append(curr_name)

	if i == 0:
            c = rna_out # for the first sample
	else:
	    for j in range(len(c)):
                c[j] = c[j]+[rna_out[j][2]]

        #print "Append sample "+curr_name+" to matrix: Done"

    # output samtools idx stat file
    outp1 = open(path_out+project_name+"_counts.txt", "w")
    outp1.write("\t".join(name1)+"\n")
    outp1.write("\n".join(map("\t".join, c)))
    outp1.write("\n")
    outp1.close()
    print "Created samtools idx stat matrix file "+path_out+project_name+"_counts.txt"


###
# Read in and process htseq files
###

def read_htseq_file(fin):
    """
    Read in output file created by htseq-count
    Return as a matrix
    """
    f = open(fin,'r')
    c = f.read().split('\n')
    if '' in c:
        c.remove('')

    # Obtain nofeature counts
    nofeature_ct = (filter(lambda x:'__no_feature' in x, c)[0]).split('\t')[1]
    # Sum up all counts
    total_ct=sum(map(lambda x: int(x.split('\t')[1]),c))
    # Create list for total count and no feature count
    list_ct=[total_ct, nofeature_ct]

    # Removing final 5 rows of htseq output that contains summary, e.g., __no_feature __ambiguous etc
    c = [x for x in c if not re.match('__*', x)]

    # Create list for htseq counts mapped to gtf annotations including genes, rRNAs, and ERCCs
    htseq_out=map(lambda x: x.split('\t'),c) # extract gene_id and count from htseq results

    return htseq_out, list_ct

def make_htseq_count_matrix(path_out, project_name, sample_names, sample_paths, ref_genome):
    """
    Read in individual output files created by htseq sample_path+"/htseq_out/"+curr_sample+"_counts.txt according to sample_info_file
    Return as a single matrix for a whole batch
    """
    # path_out: directory for output files
    # project_name: file prefix
    # sample_names: a list stores all sample names
    # sample_paths: a list stores all directory paths for each sample
    # ref_genome: reference genome

    # start header names
    ct_init = "\t".join(["Sample", "Total_count", "Nofeature_count", "Nofeature_count_Percentage"])
    name1 = ["Gene"]
    name2 = ["count"]

    a=[] # a list to store count for each sample by each gene

    for i in range(len(sample_names)):
        curr_name = sample_names[i]
	curr_path = sample_paths[i]+"htseq_out/"+curr_name+"_counts.txt"
	if not os.path.exists(curr_path):
            print "Missing htseq output file ", curr_path
	    break

	curr_htseq, ct_list = read_htseq_file(curr_path) # read in htseq outputs
        #print "Read in htseq output for sample "+sample_names[i]+": Done"

        # obtain count statistics
        ct_list=map(float, ct_list)
        total_ct=round(ct_list[0]/1000000,3)
        nofeature_ct=round(ct_list[1]/1000000,3)
        # compute nofeature count proportion
        nofeature_prop=round(ct_list[1]/ct_list[0]*100,2)
        # list of items for each line
        ct_line=[curr_name]+map(str,[total_ct, nofeature_ct, nofeature_prop])
        # add to ct_init
        ct_init = ct_init+"\n"+"\t".join(ct_line)

        if len(curr_htseq) == 0:
            print "Current sample will not be included ", curr_name
	    print "Empty htseq file ", curr_path
	    print "Rerun with complete sample set"
	    break
    
        name2_curr=map(lambda x:x+"_"+curr_name, name2)
        name1=name1+name2_curr # create header in a list: [Gene, samp1_count, samp1_tpm, sampe1_fpkm]
        
	if i == 0:
            a = curr_htseq # for the first sample
	else:
	    for j in range(len(a)):
                if a[j][0] == curr_htseq[j][0]: # check if gene name matches
                    a[j] = a[j]+curr_htseq[j][1:]
		else:
		    print "Error, gene names don't match:", curr_name, a[j], curr_htseq[j]
        #print "Append sample "+curr_name+" to matrix: Done"

    # join count, tpm, fpkm in each line by tab
    a_str =map(lambda x:'\t'.join(x), a)
    # extract ercc
    ercc_out=filter(lambda x:'ERCC' in x, a_str)

    # extract rRNA (hg38 only)
    if ref_genome=="hg38":
        # obtain hg38 rRNA reference genome
        # read in user-defined variable python script: import hg38_rRNA_gtf
        rrna_gtf=userdef.hg38_rRNA_gtf

        rrna_list=rrna_extract(rrna_gtf)
        rrna_out=map(lambda x:'\t'.join(x),(filter(lambda x:x[0] in rrna_list, a)))

        # filter out ERCC and rRNA from gene list
        gene_out=filter(lambda x:x not in ercc_out+rrna_out, a_str)

    else:
        gene_out=filter(lambda x:x not in ercc_out, a_str)

    # gene output file
    outp1 = open(path_out+project_name+"_htseq_gene.txt", "w")
    outp1.write("\t".join(name1)+"\n")
    outp1.write("\n".join(gene_out))
    outp1.close()
    print "Created gene htseq quantification matrix file "+path_out+project_name+"_htseq_gene.txt"

    # ERCC output file
    outp2 = open(path_out+project_name+"_htseq_ercc.txt", "w")
    outp2.write("\t".join(name1)+"\n")
    outp2.write("\n".join(ercc_out))
    outp2.close()
    print "Created ERCC htseq quantification matrix file "+path_out+project_name+"_htseq_ercc.txt"

    # rRNA output file (hg38 only)
    if ref_genome=="hg38":
        outp3 = open(path_out+project_name+"_htseq_rrna.txt", "w")
        outp3.write("\t".join(name1)+"\n")
        outp3.write("\n".join(rrna_out))
        outp3.close()
        print "Created rRNA htseq quantification matrix file "+path_out+project_name+"_htseq_rrna.txt"

    # output nofeature count statistics
    outp4 = open(path_out+project_name+"_htseq_nofeature.txt", "w")
    outp4.write(ct_init+"\n")
    outp4.close()
    print "Created htseq no_feature count statistics file "+path_out+project_name+"_htseq_nofeature.txt"

def rrna_extract(rrna_gtf):
    """
    Extract rRNA gene_id from gtf
    Return as a list
    """

    rrna=[]
    with open(rrna_gtf,'r') as f:
        for line in f:
            c=line.split(';')
            c=filter(lambda x:'gene_id' in x, c)[0]
            rrna=rrna+re.findall('"([^"]*)"',c)

    return rrna

###
# Read in and process cufflinks files
###

def read_cufflinks_file(fin):
    """
    Read in Cufflinks output file for ERCC spike-in FPKM values and output
    FPKM values to be put into table and plot in Rmd report
    """
    f = open(fin,'r')
    ercc_out = {}

    c = f.read().split('\n')[1:] # exclude header
    if '' in c:
        c.remove('')

    # Create directory for cufflinks counts mapped to gtf annotations
    cufflinks_out={} # key: gene_id; value: fpkm
    for x in c:
        cufflinks_out[x.split('\t')[0]]=[x.split('\t')[9]]

    return cufflinks_out

def make_cufflinks_matrix(path_out, project_name, sample_names, sample_paths, suffix):
    """
    Read in individual output files created by cufflinks sample_path+"/cufflinks_out"+suffix+"/genes.fpkm_tracking" according to sample_info_file
    Return as a single matrix for a whole batch
    """
    # path_out: directory for output files
    # project_name: file prefix
    # sample_names: a list stores all sample names
    # sample_paths: a list stores all directory paths for each sample
    # suffix: directory "cufflinks_out" suffix ("_ERCC" for ercc, or "" for gene)

    name1 = ["Gene"]
    name2 = ["fpkm"]

    # Create directory for cufflinks counts for all samples
    a = {} # key: gene_id, value: fpkm for all samples

    for i in range(len(sample_names)):
        curr_name = sample_names[i]
	curr_path = sample_paths[i]+"cufflinks_out"+suffix+"/genes.fpkm_tracking"
	if not os.path.exists(curr_path):
	    print "Missing cufflinks ercc output file ", curr_path
	    break

        curr_cufflinks = read_cufflinks_file(curr_path) # read in cufflinks outputs
        print "Read in cufflinks output for sample "+sample_names[i]+": Done"

        if len(curr_cufflinks) == 0:
            print "Current sample will not be included ", curr_name
	    print "Empty cufflinks file ", curr_path
	    print "Rerun with complete sample set"
	    break
    
        name2_curr=map(lambda x:x+"_"+curr_name, name2)
        name1=name1+name2_curr # create header in a list: [Gene, samp1_count, samp1_tpm, sampe1_fpkm]

        if i == 0:
            a = curr_cufflinks # for the first sample
        else:
            for j in a:
                a[j]=a[j]+curr_cufflinks[j]

        #print "Append sample "+curr_name+" to matrix: Done"

    # output cufflinks fpkm for all samples
    gene_out=[]
    for key in a:
        gene_out=gene_out+[key+"\t"+"\t".join(a[key])]

    # output cufflinks fpkm file
    if "_ERCC" in suffix: # ERCC
        outp1 = open(path_out+project_name+"_cufflinks_ercc.txt", "w")
	outp1.write("\t".join(name1)+"\n")
	outp1.write("\n".join(gene_out))
	outp1.close()
	print "Created ercc cufflinks quantification matrix file "+path_out+project_name+"_cufflinks_ercc.txt"
    else: # gene
        outp1 = open(path_out+project_name+"_cufflinks_gene.txt", "w")
	outp1.write("\t".join(name1)+"\n")
	outp1.write("\n".join(gene_out))
	outp1.close()
	print "Created gene cufflinks quantification matrix file "+path_out+project_name+"_cufflinks_gene.txt"

###
#  Read in and process tophat log files
###

def read_tophat_logs(fin, library_type):
    """
    Read in TopHat output files of type (prep_reads.info) to get number of raw reads passed to program and number that were aligned
    Returns output for Rmd report
    """
    f = open(fin,'r')
    c = f.read().split('\n')
    #if '' in c:
	#c.remove('')
    #First entry is 'left_reads_in'; Second entry is 'left_reads_out'
    #For paired-end data, third entry is 'right_reads_in'; fourth entry is 'right_reads_out'
    if library_type in ["PE"]:
	read_numbers = map(lambda x: x.split('=')[1], c[2:4]+c[6:])
    else:
	read_numbers = map(lambda x: x.split('=')[1], c[2:4])
    return read_numbers

def make_tophat_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type):
    outp5 = open(path_out+project_name+"_read_counts.txt", "w")
    if library_type in ["PE"]:
	name5 = ["Sample", "R1_Raw_Reads", "R1_Aligned_Reads", "R2_Raw_Reads", "R2_Aligned_Reads"]
    else:
	name5 = ["Sample", "Raw_Reads", "Aligned_Reads"]
    for i in range(len(sample_names)):
	curr_name = sample_names[i]
	curr_path = sample_paths[i]
	total_reads = read_tophat_logs(curr_path+aligner+"_out/prep_reads.info", library_type)
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

###
# Read in and process bamstats files
###

def read_bamtools_stats(fin):
	"""
	Read in bamstats output file accepted_hits.sorted.bamstats
	Reformats and outputs portions of interest for Rmd report
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

def make_bamstats_matrix(path_out, project_name, sample_names, sample_paths, aligner):
	outp6 = open(path_out+project_name+"_bamstats_counts.txt", "w")
	name6 = ["Type"]
	for i in range(len(sample_names)):
		curr_name = sample_names[i]
		curr_path = sample_paths[i]
		name6.append(curr_name)
		bam_stats = read_bamtools_stats(curr_path+aligner+"_out/"+curr_name+"_accepted_hits.sorted.bamstats")
		if i == 0:
			f = bam_stats
		else:
			for j in range(len(bam_stats)):
				f[j] = list(f[j])+[bam_stats[j][1]]
	outp6.write("\t".join(name6)+"\n")
	outp6.write("\n".join(map("\t".join, f)))
        outp6.write("\n")
	outp6.close()
	print "Created file "+path_out+project_name+"_bamstats_counts.txt"


###
# Read in and process _ReadCount from raw .fastq files
###

def get_unique_reads(fin, library_type):
	"""
	Read in file that has output from awk on number of reads, unique reads, and %unique reads from fastq files
	Output R1 number of reads, R1 unique reads, R1 %unique reads, R2 number of reads, R2 unique reads, and R2 %unique reads for Rmd report
	"""
	f = open(fin,'r')
	c = f.read().split('\n')
	if '' in c:
		c.remove('')
	#First row is R1. Second row is R2. 
	read_numbers = map(lambda x: x.split(' '), c)
	if library_type in ["PE"]:
		read_numbers = read_numbers[0]+read_numbers[1]
	else:
		read_numbers = read_numbers[0]
	return read_numbers

def make_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type):
    outp9 = open(path_out+project_name+"_unique_counts.txt", "w")
    if library_type in ["PE"]:
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
    outp9.write("\n")
    outp9.close()
    print "Created file "+path_out+project_name+"_unique_counts.txt"


###
# Read in and process fastqc files
###

def read_fastq_data(fin):
	"""
	Read fastqc_data.txt file from a FastQC report and extract the percentage of duplicates and histogram information
	"""
	f = open(fin,'r')
	c = f.read().split('#Total Deduplicated Percentage')[1]
	c2 = c.split('>>END_MODULE')[0]
	c3 = c2.split('\n')
	if '' in c3:
		c3.remove('')
	#Each entry contains Duplication Level, Percentage of deduplicated, Percentage of total
        duplicate_data =  [['Total Deduplicated Percentage', c3[0].strip('\t')]]+map(lambda x: x.split('\t'), c3[2:]) # c3[0]: Total Deduplicated Percentage; c3[1]: header
	return duplicate_data

def make_duplicate_matrix(path_out, project_name, sample_names, sample_paths, library_type):

    #Duplicate read info from fastq files
    outp10 = open(path_out+project_name+"_duplicates.txt", "w")
    name10 = ["Read_Number"]
    for i in range(len(sample_names)):
        curr_name = sample_names[i]
        curr_path = sample_paths[i]
	#Move individual FastQC reports to report directory
        cmd_cp="cp "+curr_path+curr_name+"_R1_Trimmed_fastqc.zip "+path_out
	subprocess.Popen(cmd_cp, shell=True).wait()
        cmd_unzip="unzip -o -q -d "+path_out+" "+path_out+curr_name+"_R1_Trimmed_fastqc.zip"
	subprocess.Popen(cmd_unzip, shell=True).wait()
        cmd_rm="rm  "+path_out+curr_name+"_R1_Trimmed_fastqc.zip"
	subprocess.Popen(cmd_rm, shell=True).wait()

	name10.append(curr_name+"_R1")
	if os.path.isfile(path_out+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"
	else:
	    print "Missing FastQC report", curr_name
	    sys.exit()

	curr_dup_R1 = read_fastq_data(fastqc1)

	if library_type in ["PE"]:
            cmd_cp="cp "+curr_path+curr_name+"_R2_Trimmed_fastqc.zip "+path_out
    	    subprocess.Popen(cmd_cp, shell=True).wait()
            cmd_unzip="unzip -o -q -d "+path_out+" "+path_out+curr_name+"_R2_Trimmed_fastqc.zip"
	    subprocess.Popen(cmd_unzip, shell=True).wait()
            cmd_rm="rm  "+path_out+curr_name+"_R2_Trimmed_fastqc.zip"
	    subprocess.Popen(cmd_rm, shell=True).wait()

            name10.append(curr_name+"_R2")
	    if os.path.isfile(path_out+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"):
                fastqc2 = path_out+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"
	    else:
		print "Missing FastQC report for R2 ", curr_name
                sys.exit()
	    curr_dup_R2 = read_fastq_data(fastqc2)

	if i == 0:
	    duplicates = curr_dup_R1
	    for j in range(1, len(duplicates)):
                duplicates[j] = [duplicates[j][0], duplicates[j][2]]
	    if library_type in ["PE"]:
                duplicates[0].append(curr_dup_R2[0][1])
		for j in range(1, len(duplicates)):
                    duplicates[j].append(curr_dup_R2[j][2])
	else:
	    duplicates[0].append(curr_dup_R1[0][1])
	    for j in range(1, len(duplicates)):
	        duplicates[j].append(curr_dup_R1[j][2])
	    if library_type in ["PE"]:
	        duplicates[0].append(curr_dup_R2[0][1])
		for j in range(1, len(duplicates)):
                    duplicates[j].append(curr_dup_R2[j][2])

    outp10.write("\t".join(name10)+"\n")
    outp10.write("\n".join(map("\t".join, duplicates)))
    outp10.write("\n")
    outp10.close()
    print "Created file "+path_out+project_name+"_duplicates.txt"


def make_project_data_files(project_name, sample_names, sample_paths, path_out, ref_genome, library_type, aligner):
	"""
	Creates several text files to be loaded by R for Rmd report based on modified program outputs read in with preceding scripts
	These text files are matrices containing information for all samples in a project/batch
	Currently, cycles through all samples multiple times, once to create each individual file type
	Currently, there is no way to handle missing files. If an error is encountered the process will stop at that point and not complete
	Currently, assumes default naming convention of all programs used in rnaseq_align_and_qc.py
	"""
	#Read counts obtained from TopHat log files -- outp5
	if aligner == "tophat":
                make_tophat_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type)

	#Unique read counts obtained by comprehensive count of fastq files -- outp9
        make_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type)
	
	#Duplicate read info from fastq files -- outp10
        make_duplicate_matrix(path_out, project_name, sample_names, sample_paths, library_type)
	
	#rnaseqmetrics output split into two files: 1) overall metrics (2) data to create normalized coverage histogram for all samples
        # --outp1 and outp2
        make_RnaSeqMetrics_matrix(path_out, project_name, sample_names, sample_paths, aligner)

	#samtools stats counts of reads per chromosome
	#ercc transcript read counts - combo of samtools stats output and cufflinks run with ERCC gtf file
        make_samidxstat_matrix(path_out, project_name, sample_names, sample_paths, ref_genome, aligner)

        if aligner=="tophat":
            make_cufflinks_matrix(path_out, project_name, sample_names, sample_paths, "_ERCC")
        if aligner=="star":
            make_htseq_count_matrix(path_out, project_name, sample_names, sample_paths, ref_genome)

	#bamstats output metrics on types of reads, including junction spanning reads
        make_bamstats_matrix(path_out, project_name, sample_names, sample_paths, aligner)

	#insertsizemetrics output on insert size statistics -- outp6
        make_insertsize_matrix(path_out, project_name, sample_names, sample_paths, library_type, aligner)

def make_rmd_css(path_start):
    """
    create custom.css for rmarkdown
    """

    css_outp = open(path_start+"custom.css", "w")
    css_outp.write("blockquote {\n")
    css_outp.write("    padding: 10px 20px;\n")
    css_outp.write("    margin: 0 0 20px;\n")
    css_outp.write("    font-size: 14px;\n")
    css_outp.write("    background-color: #eee;\n")
    css_outp.write("    border-left: 5px solid #eee;\n")
    css_outp.write("}\n\n")
    css_outp.write(".main-container {\n")
    css_outp.write("    max-width: 2000px !important;\n") # "!important" overrides other rules
    css_outp.write("}\n")
    css_outp.close()


def make_rmd_title(path_start, project_name, aligner, library_type, strand, ref_genome):
    """
    Rmarkdown creation: Create title and major description
    """

    # read in user-defined variable python script
    # import software version
    trimmomatic_version=userdef.trimmomatic_version
    fastqc_version=userdef.fastqc_version
    star_version=userdef.star_version
    samtools_version=userdef.samtools_version
    bamtools_version=userdef.bamtools_version
    picard_version=userdef.picard_version
    # import author information
    author=userdef.author

    rmd="---\ntitle: 'RNA-Seq Report of Sample QC and Alignment Summary Statistics for "+project_name+"'\n"
    rmd=rmd+"author: "+author+"\n"
    rmd=rmd+"date: \"`r format(Sys.time(), '%d %B, %Y')`\"\n"
    rmd=rmd+"output: \n"
    rmd=rmd+"  html_document:\n"
    rmd=rmd+"    css: custom.css\n"
    rmd=rmd+"    toc: true\n"
    rmd=rmd+"    toc_float: true\n---\n\n"

    # Variable used
    rmd=rmd+"**Project:** "+project_name+"\n\n"
    if aligner == "star":
        align_abbrev = "STAR ("+star_version+")"
    elif aligner=="tophat":
	align_abbrev = "Tophat"
    rmd=rmd+"**Aligner:** "+align_abbrev+"\n\n"
    if ref_genome == "hg19":
        rmd=rmd+"**Genome:** For human, the hg19 assembly was used. We estimate the number of rRNA reads as those mapped to chrM plus chrUn_gl000220, corresponding to 12S, 16S and 5.8S rRNA. The 'Other' category contains all other chr*_random and chrUn_* available. If using the 2014 updated version of the hg19 files, these categories are no longer present.\n"
    elif ref_genome == "hg38":
        rmd=rmd+"**Genome:** For human, the hg38 assembly was used. We estimate the number of rRNA reads as those mapped to chrM plus chrUn_GL000220v1, corresponding to 12S, 16S and 5.8S rRNA. The 'Other' category contains all other chr*_random and chrUn_* available.\n"
    elif ref_genome == "mm38":
	rmd=rmd+"**Genome:** For mouse, the ENSEMBL GRCm38 assembly available in iGenomes was used.\n"
    elif ref_genome == "mm10":
	rmd=rmd+"**Genome:** For mouse, the UCSC mm10 assembly available in iGenomes was used.\n"
    elif ref_genome == "rn6":
	rmd=rmd+"**Genome:** For rat, the rn6 assembly was used.\n"
    elif ref_genome == "susScr3":
	rmd=rmd+"**Genome:** For pig, the susScr3 assembly was used.\n"
    elif ref_genome == "Zv9":
	rmd=rmd+"**Genome:** For zebrafish, the Zv9 assembly comprises a sequence length of 1.4 Gb in 26 chromosomes (labels 1-25 and MT) and 1,107 scaffolds (merged into label 'Other').\n"
    rmd=rmd+"\n\n"

    # Bioinformatics tools
    rmd=rmd+"**Informatics tools used:**\n\n"
    rmd=rmd+"* Trimmomatic ("+trimmomatic_version+")\n"
    rmd=rmd+"* FastQC ("+fastqc_version+")\n"
    if aligner=="star":
        rmd=rmd+"* STAR ("+star_version+")\n"
    rmd=rmd+"* samtools ("+samtools_version+")\n"
    rmd=rmd+"* bamtools ("+bamtools_version+")\n"
    rmd=rmd+"* Picard Tools ("+picard_version+")\n"
    rmd=rmd+"\n\n"

    # Sequencing parameters
    rmd=rmd+"**Sequencing parameters:**\n\n"
    rmd=rmd+"* library_type = "+library_type+"\n"
    rmd=rmd+"* strand = "+strand+"\n"
    rmd=rmd+"* ref_genome = "+ref_genome+"\n\n"

    return rmd

def make_rmd_var(path_start, project_name, ref_genome, library_type, aligner, sample_names, ercc_mixes, sample_info_file, template_dir):
    """
    Rmarkdown creation: Define variables and files
    """

    rmd="```{r vars, echo=F}\n"
    rmd=rmd+"project_name=\""+project_name+"\"\n"
    if "./" in path_start:
        rmd=rmd+"path.start=\""+path_start.lstrip("./")+"\"\n"
    else:
	rmd=rmd+"path.start=\""+path_start+"\"\n"

    rmd=rmd+"sample.names.orig <- c("+str(sample_names)[1:-1]+")\n" # create original sample names
    rmd=rmd+"sample.names <- sample.names.orig\n"
    if ercc_mixes is not None:
        rmd=rmd+"ercc.mixes=c("+str(ercc_mixes)[1:-1]+")\n" # convert list to string
    rmd=rmd+"genome=\""+ref_genome+"\"\n"
    rmd=rmd+"library_type=\""+library_type+"\"\n"
    rmd=rmd+"aligner=\""+aligner+"\"\n"
    rmd=rmd+"sample_info_file='"+sample_info_file+"'\n"
    if aligner=="star":
        rmd=rmd+"count_data_file='"+path_start+project_name+"_htseq_gene.txt'\n"

    # define files
    if ercc_mixes is not None:
        rmd=rmd+"ambion_file <- '"+template_dir+"ERCC_SpikeIn_Controls_Analysis.txt'\n"
        if aligner=="tophat":
            rmd=rmd+"ercc.data <- read.table('"+project_name+"_cufflinks_ercc.txt', sep='\\t', header=T, as.is=T)\n"
        elif aligner=="star" and ercc_mixes is not None:
            rmd=rmd+"ercc.data <- read.table('"+project_name+"_htseq_ercc.txt', sep='\\t', header=T, as.is=T)\n"    
    rmd=rmd+"```\n\n"

    return rmd



def make_rmd_featurestat(project_name):
    """
    Rmarkdown creation: Obtain proportion of no feature counts from htseq-count results
    Separate here because of the htseq-count-specific feature
    """
    rmd="\n"
    rmd=rmd+"## HTSeq-count: No feature counts statistics\n\n"
    rmd=rmd+"Numbers of reads that can not mapped to any feature (Nofeature count) are shown by per million reads from htseq-count quantification results\n\n"
    rmd=rmd+"```{r nofeature, eval=T, echo=F, message=F, warning=F, results='asis'}\n"
    rmd=rmd+"nofeature.data <- read.table('"+project_name+"_htseq_nofeature.txt', sep='\\t', header=T, as.is=T)\n"
    rmd=rmd+"DT::datatable(nofeature.data, rownames=FALSE, options = list(pageLength = 25))\n"
    rmd=rmd+"```\n\n"
    return rmd

def lsf_file(job_name, cmd, memory=36000, thread=1, queue=userdef.queue):
    """
    Creates .lsf files
    """

    outp = open(job_name+".lsf",'w')
    outp.write("#!/bin/bash\n")
    outp.write("#BSUB -L /bin/bash\n")
    outp.write("#BSUB -J "+job_name+"\n")
    outp.write("#BSUB -q "+queue+"\n")
    outp.write("#BSUB -o "+job_name+"_%J.out\n")
    outp.write("#BSUB -e "+job_name+"_%J.screen\n")
    outp.write("#BSUB -M "+str(memory)+"\n")
    outp.write("#BSUB -n "+str(thread)+"\n")
    outp.write(cmd)
    outp.write("\n")
    outp.close()

def make_rmd_html(sample_info_file, project_name, path_start, sample_names, ercc_mixes, ref_genome, library_type, strand, aligner, template_dir):
    """
    Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
    Runs rmarkdown to create html document
    Also makes a custom css format file - improvement on rmarkdown defaults
    """
    # create custom.css for rmarkdown
    make_rmd_css(path_start)

    # create contents for rmarkdown RMD file
    rmd = ""

    # create title and description
    rmd_title =make_rmd_title(path_start, project_name, aligner, library_type, strand, ref_genome)
    rmd = rmd + rmd_title

    # create varialbes
    rmd_var = make_rmd_var(path_start, project_name, ref_genome, library_type, aligner, sample_names, ercc_mixes, sample_info_file, template_dir)
    rmd = rmd + rmd_var

    # read in contents in template Rmdfile
    rmd_in = open(template_dir+"rnaseq_align_and_qc_report_Rmd_template.txt", "r")
    rmd_template = rmd_in.read()
    rmd = rmd + rmd_template
    rmd = rmd + "\n\n"

    # output no feature count stat table if htseq-count is used
    if aligner=="star":
       rmd_featurestat = make_rmd_featurestat(project_name)
       rmd = rmd + rmd_featurestat

    # write in project_name+"_QC_RnaSeqReport.Rmd" file
    outp = open(path_start+project_name+"_QC_RnaSeqReport.Rmd", "w")
    outp.write(rmd)
    # output session info
    outp.write("```{r sessioninfo, echo=F}\n")
    outp.write("pander(sessionInfo())\n")
    outp.write("```\n\n")
    outp.close()
    print "Created file "+path_start+project_name+"_QC_RnaSeqReport.Rmd"

    # create .lsf file. Run it on HPC becuase rlog for all sample counts in pca step is computationally demanding.
    cmd="cd "+path_start+"; echo \"library(rmarkdown); rmarkdown::render('"+project_name+"_QC_RnaSeqReport.Rmd')\" | R --no-save --no-restore\n"
    lsf_file(project_name+"_qc", cmd)
    print "Created LSF script "+project_name+"_qc.lsf in current directory"

def main(project_name, sample_info_file, path_start, aligner, ref_genome, library_type, strand, template_dir):
    """
    Creates html report describing summary and QC statistics for a set of aligned RNA-Seq samples associated with a project
    Report is based on multiple output files created by rnaseq_align_and_qc.py
    Such files are first reformatted into matrices that are easily loaded into R
    Input:
        project_name: name for report.
	sample_info_file: tab delimited txt file with sample information as described in rnaseq_align_and_qc.py
	rmd_template: txt file that contains most of the contents that will populate the Rmd file for the report
    Current ref_genome choices: hg19, mm38
    Current library_type choices: PE, SE, DGE, SPE
    """
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
	path_start = path_start+"/"
    new_dir = path_start+project_name+"_Alignment_QC_Report_"+aligner+"/"
    if not os.path.exists(new_dir):
	os.makedirs(new_dir)
    if template_dir == "./":
        template_dir = os.getcwd()
    if template_dir[-1] != "/":
        template_dir = template_dir+"/"

    # check if QC template txt file exists
    if not os.path.exists(template_dir+"rnaseq_align_and_qc_report_Rmd_template.txt"):
        print "Cannot find "+template_dir+"rnaseq_align_and_qc_report_Rmd_template.txt"
	sys.exit()

    #Get list of dictionary of sample information. Keys: [sample_id, ercc_mix, index]
    info_dict = get_sample_info(sample_info_file)
    sample_names = info_dict["Sample"]
    sample_paths = map(lambda x: path_start+x+"/", sample_names)

    if "ERCC_Mix" in info_dict:
        ercc_mixes = info_dict["ERCC_Mix"]
        # ERCC mix exists. Use ERCC concentration file
        if not os.path.exists(template_dir+"ERCC_SpikeIn_Controls_Analysis.txt"):
            print "ERCC spike-in used. Cannot find ERCC concentration file "+template_dir+"ERCC_SpikeIn_Controls_Analysis.txt"
            sys.exit()

    else:
        ercc_mixes=None

    # create stats files used for QC report
    make_project_data_files(project_name, sample_names, sample_paths, new_dir, ref_genome, library_type, aligner)

    #Create the report
    make_rmd_html(sample_info_file, project_name, new_dir, sample_names, ercc_mixes, ref_genome, library_type, strand, aligner, template_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create HTML report of QC and alignment summary statistics for RNA-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of for output directory and files")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
    parser.add_argument("--aligner", default="star", type=str, help="Should TopHat or STAR be used as aligner (options: star)")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome (options: hg38, hg19, mm38, mm10, rn6, susScr3)")
    parser.add_argument("--library_type", default="PE", type=str, help="Specify library type (options: PE (paired-end), SE (single-end))")
    parser.add_argument("--strand", type=str, default="nonstrand", help="Whether data is from a strand-specific assay."
        "(options: nonstrand, reverse, forward)")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put template RMD file rnaseq_align_and_qc_report_Rmd_template.txt for QC report")
    args = parser.parse_args()

    if args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.project_name, args.samples_in, args.path_start, args.aligner, args.ref_genome, args.library_type, args.strand, args.template_dir)
