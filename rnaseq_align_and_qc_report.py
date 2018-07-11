#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import re

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
	v7: library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE", "SSE", "SPE-PrepX",
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end",
							"stranded single-end", "stranded paired-end from PrepX")
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


#ERCC gtf file
ERCC_only = "/project/bhimeslab/Reference/ERCC92.gtf"

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

	if library_type in ["PE", "SPE", "SPE-PrepX"]:
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
        print "Read in samtools idxstat output for sample "+curr_name+": Done"
        name1.append(curr_name)

	if i == 0:
            c = rna_out # for the first sample
	else:
	    for j in range(len(c)):
                c[j] = c[j]+[rna_out[j][2]]

        print "Append sample "+curr_name+" to matrix: Done"

    # output samtools idx stat file
    outp1 = open(path_out+project_name+"_counts.txt", "w")
    outp1.write("\t".join(name1)+"\n")
    outp1.write("\n".join(map("\t".join, c)))
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
    ct_init = "\t".join(["Sample", "Total_count", "Nofeature_count", "Proportion"])
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
        print "Read in htseq output for sample "+sample_names[i]+": Done"

        # obtain count statistics
        ct_list=map(float, ct_list)
        total_ct=round(ct_list[0]/1000000,2)
        nofeature_ct=round(ct_list[1]/1000000,2)
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
        print "Append sample "+curr_name+" to matrix: Done"

    # join count, tpm, fpkm in each line by tab
    a_str =map(lambda x:'\t'.join(x), a)
    # extract ercc
    ercc_out=filter(lambda x:'ERCC' in x, a_str)

    # extract rRNA
    rrna_list=rrna_extract(ref_genome)
    rrna_out=map(lambda x:'\t'.join(x),(filter(lambda x:x[0] in rrna_list, a)))

    # filter out ERCC and rRNA from gene list
    gene_out=filter(lambda x:x not in ercc_out+rrna_out, a_str)

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

    # rRNA output file
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

def rrna_extract(ref_genome):
    """
    Extract rRNA gene_id from gtf
    Return as a list
    """

    # obtain rRNA reference genome
    ref_index, fa, gtf, ref, ERCC_gtf, mask_gtf, genome_dir = get_genome_ref_files(ref_genome, "rRNA")

    rrna=[]

    with open(mask_gtf,'r') as f:
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
	curr_path = sample_paths[i]+"/cufflinks_out"+suffix+"/genes.fpkm_tracking"
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

        print "Append sample "+curr_name+" to matrix: Done"

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
    if library_type in ["PE", "SPE", "SPE-PrepX"]:
	read_numbers = map(lambda x: x.split('=')[1], c[2:4]+c[6:])
    else:
	read_numbers = map(lambda x: x.split('=')[1], c[2:4])
    return read_numbers

def make_tophat_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type):
    outp5 = open(path_out+project_name+"_read_counts.txt", "w")
    if library_type in ["PE", "SPE", "SPE-PrepX"]:
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
		bam_stats = read_bamtools_stats(curr_path+"/"+aligner+"_out/"+curr_name+"_accepted_hits.sorted.bamstats")
		if i == 0:
			f = bam_stats
		else:
			for j in range(len(bam_stats)):
				f[j] = list(f[j])+[bam_stats[j][1]]
	outp6.write("\t".join(name6)+"\n")
	outp6.write("\n".join(map("\t".join, f)))
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
	if library_type in ["PE", "SPE", "SPE-PrepX"]:
		read_numbers = read_numbers[0]+read_numbers[1]
	else:
		read_numbers = read_numbers[0]
	return read_numbers

def make_readcount_matrix(path_out, project_name, sample_names, sample_paths, library_type):
    outp9 = open(path_out+project_name+"_unique_counts.txt", "w")
    if library_type in ["PE", "SPE", "SPE-PrepX"]:
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
	duplicate_data =  [['Total Deduplicated Percentage', c3[0].strip('\t')]]+map(lambda x: x.split('\t'), c3[2:])
	return duplicate_data

def make_duplicate_matrix(path_out, project_name, sample_names, sample_paths, library_type):

    #Duplicate read info from fastq files
    outp10 = open(path_out+project_name+"_duplicates.txt", "w")
    name10 = ["Read_Number"]
    for i in range(len(sample_names)):
        curr_name = sample_names[i]
        curr_path = sample_paths[i]
	#Move individual FastQC reports to report directory
	subprocess.call("cp "+curr_path+curr_name+"_R1_Trimmed_fastqc.zip "+path_out, shell=True)
	subprocess.call("unzip -o -q -d "+path_out+" "+path_out+curr_name+"_R1_Trimmed_fastqc.zip", shell=True)
	subprocess.call("rm  "+path_out+curr_name+"_R1_Trimmed_fastqc.zip", shell=True)
	name10.append(curr_name+"_R1")	 
	if os.path.isfile(path_out+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_R1_Trimmed_fastqc/fastqc_data.txt"
	elif os.path.isfile(path_out+"/"+curr_name+"_R1_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_R1_fastqc/fastqc_data.txt"
	elif os.path.isfile(path_out+"/"+curr_name+"_fastqc/fastqc_data.txt"):
	    fastqc1 = path_out+"/"+curr_name+"_fastqc/fastqc_data.txt"
	else:
	    print "Missing FastQC report", curr_name
	    break

	curr_dup_R1 = read_fastq_data(fastqc1)

	if library_type in ["PE", "SPE", "SPE-PrepX"]:
	    subprocess.call("cp "+curr_path+curr_name+"_R2_Trimmed_fastqc.zip "+path_out, shell=True)
	    subprocess.call("unzip -o -q -d "+path_out+" "+path_out+curr_name+"_R2_Trimmed_fastqc.zip", shell=True)
	    subprocess.call("rm  "+path_out+curr_name+"_R2_Trimmed_fastqc.zip", shell=True)
	    name10.append(curr_name+"_R2")
	    if os.path.isfile(path_out+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"):
                fastqc2 = path_out+"/"+curr_name+"_R2_Trimmed_fastqc/fastqc_data.txt"
	    elif os.path.isfile(path_out+"/"+curr_name+"_R2_fastqc/fastqc_data.txt"):
		fastqc2 = path_out+"/"+curr_name+"_R2_fastqc/fastqc_data.txt"
	    else:
		print "Missing FastQC report for R2 ", curr_name
	    curr_dup_R2 = read_fastq_data(fastqc2)

	if i == 0:
	    duplicates = curr_dup_R1
	    for j in range(1, len(duplicates)):
                duplicates[j] = [duplicates[j][0], duplicates[j][2]]
	    if library_type in ["PE", "SPE", "SPE-PrepX"]:
                duplicates[0].append(curr_dup_R2[0][1])
		for j in range(1, len(duplicates)):
                    duplicates[j].append(curr_dup_R2[j][2])
	else:
	    duplicates[0].append(curr_dup_R1[0][1])
	    for j in range(1, len(duplicates)):
	        duplicates[j].append(curr_dup_R1[j][2])
	    if library_type in ["PE", "SPE", "SPE-PrepX"]:
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

def make_rmd_html(rmd_template, project_name, path_start, sample_names, ercc_mixes, ref_genome, library_type, aligner, template_dir):
    """
    Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
    Runs rmarkdown to create html document
    Also makes a custom css format file - improvement on rmarkdown defaults
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
    outp = open(path_start+project_name+"_QC_RnaSeqReport.Rmd", "w")
    outp.write("---\ntitle: \'RNA-Seq Report of Sample QC and Alignment Summary Statistics for Project "+project_name+"'\n")
    outp.write("author: 'Blanca Himes (bhimes@upenn.edu)'\n")
    outp.write("output: \n")
    outp.write("  html_document:\n")
    outp.write("    css: custom.css\n")
    outp.write("    toc: true\n")
    outp.write("    toc_float: true\n---\n\n")
    outp.write("**Project:** "+project_name+"\n\n")
    if aligner == "star":
        align_abbrev = "STAR (v. 2.5.2b)"
    else:
	align_abbrev = "Tophat"
    outp.write("**Aligner:** "+align_abbrev+"\n\n")
    if ref_genome == "hg19":
        outp.write("**Genome:** For human, the hg19 assembly was used. We estimate the number of rRNA reads as those mapped to chrM plus chrUn_gl000220, corresponding to 12S, 16S and 5.8S rRNA. The 'Other' category contains all other chr*_random and chrUn_* available. If using the 2014 updated version of the hg19 files, these categories are no longer present.\n")
    elif ref_genome == "hg38":
        outp.write("**Genome:** For human, the hg38 assembly was used. We estimate the number of rRNA reads as those mapped to chrM plus chrUn_GL000220v1, corresponding to 12S, 16S and 5.8S rRNA. The 'Other' category contains all other chr*_random and chrUn_* available.\n")
    elif ref_genome == "mm38":
	outp.write("**Genome:** For mouse, the ENSEMBL GRCm38 assembly available in iGenomes was used.\n")
    elif ref_genome == "mm10":
	outp.write("**Genome:** For mouse, the UCSC mm10 assembly available in iGenomes was used.\n")
    elif ref_genome == "rn6":
	outp.write("**Genome:** For rat, the rn6 assembly was used.\n")
    elif ref_genome == "susScr3":
	outp.write("**Genome:** For pig, the susScr3 assembly was used.\n")
    elif ref_genome == "Zv9":
	outp.write("**Genome:** For zebrafish, the Zv9 assembly comprises a sequence length of 1.4 Gb in 26 chromosomes (labels 1-25 and MT) and 1,107 scaffolds (merged into label 'Other').\n")
    outp.write("```{r vars, eval=T, echo=F, message=F, warning=F}\n")
    outp.write("project_name=\""+project_name+"\"\n")
    if "./" in path_start:
        outp.write("path.start=\""+path_start.lstrip("./")+"\"\n")
    else:
	outp.write("path.start=\""+path_start+"\"\n")
    outp.write("ercc.mixes=c("+str(ercc_mixes)[1:-1]+")\n") # convert list to string
    outp.write("sample.names.orig=c("+str(sample_names)[1:-1]+")\n") # create original sample names
    outp.write("sample.names <- sample.names.orig\n")
    outp.write("genome=\""+ref_genome+"\"\n")
    outp.write("library.type=\""+library_type+"\"\n")
    outp.write("aligner=\""+aligner+"\"\n")
    outp.write("\n```\n\n")

    # define files
    outp.write("```{r files, eval=T, echo=F}\n")
    if aligner=="tophat":
        outp.write("ercc.data <- read.table('"+project_name+"_cufflinks_ercc.txt', sep='\\t', header=T, as.is=T)\n")
    elif aligner=="star":
        outp.write("ercc.data <- read.table('"+project_name+"_htseq_ercc.txt', sep='\\t', header=T, as.is=T)\n")
    outp.write("ambion_file <- '"+template_dir+"ERCC_SpikeIn_Controls_Analysis.txt'\n")
    outp.write("```\n\n")

    outp.write(rmd_template)

    # output no feature count stat table if htseq-count is used
    if aligner=="star":
        outp.write("\n")
        outp.write("## HTSeq-count: No feature counts statistics\n\n")
        outp.write("No feature count (per million reads) statistics from htseq-count quantification results\n\n")
        outp.write("```{r nofeature, eval=T, echo=F, message=F, warning=F, results='asis'}\n")
        outp.write("nofeature.data <- read.table('"+project_name+"_htseq_nofeature.txt', sep='\\t', header=T, as.is=T)\n")
        outp.write("DT::datatable(nofeature.data)\n")
        outp.write("```\n\n")

    # output session info
    outp.write("```{r sessioninfo, eval=T, echo=F}\n")
    outp.write("pander(sessionInfo())\n")
    outp.write("```\n\n")

    outp.close()

    #subprocess.call("cd "+path_start+"; echo \"library(R2HTML); Sweave('"+project_name+"_QC_RnaSeqReport.Rnw', driver=RweaveHTML)\" | R --no-save --no-restore", shell=True)
    #subprocess.call("cd "+path_start+"; echo \"library(knitr); library(markdown); knit2html('"+project_name+"_QC_RnaSeqReport.Rmd', force_v1 = TRUE, options = c('toc', markdown::markdownHTMLOptions(TRUE)))\" | R --no-save --no-restore", shell=True)
    subprocess.call("cd "+path_start+"; echo \"library(rmarkdown); rmarkdown::render('"+project_name+"_QC_RnaSeqReport.Rmd')\" | R --no-save --no-restore", shell=True)

def main(project_name, sample_info_file, path_start, aligner, template_dir):
    """
    Creates html report describing summary and QC statistics for a set of aligned RNA-Seq samples associated with a project
    Report is based on multiple output files created by rnaseq_align_and_qc.py
    Such files are first reformatted into matrices that are easily loaded into R
    Report is created using Sweave
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
    new_dir = path_start+project_name+"/"+project_name+"_Alignment_QC_Report_"+aligner+"/"
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

    #Get list of sample info. Fields: [curr_sample, index, ercc_mix, top_dir, project, label, ref_genome, library_type]
    runs = get_sample_info(sample_info_file)
    sample_names = map(lambda x: x[0], runs)
    sample_paths = map(lambda x: path_start+x[4]+"/"+x[0]+"/", runs)
    ercc_mixes = map(lambda x: x[2], runs)
    ref_genome_list = map(lambda x: x[6], runs)
    library_type_list = map(lambda x: x[7], runs)
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

    # create stats files used for QC report
#    make_project_data_files(project_name, sample_names, sample_paths, new_dir, ref_genome, library_type, aligner)

    #Create the report
    if not os.path.exists(template_dir+"rnaseq_align_and_qc_report_Rmd_template.txt"):
        print "Cannot find rnaseq_align_and_qc_report_Rmd_template.txt"
	sys.exit()

    rmd_in = open(template_dir+"rnaseq_align_and_qc_report_Rmd_template.txt", "r")
    rmd_template = rmd_in.read()
    make_rmd_html(rmd_template, project_name, new_dir, sample_names, ercc_mixes, ref_genome, library_type, aligner, template_dir)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Create HTML report of QC and alignment summary statistics for RNA-seq samples associated with a project.")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
	parser.add_argument("--aligner", default="tophat", type=str, help="Should TopHat or STAR be used as aligner?"
		"(options: tophat, star)")
        parser.add_argument("--template_dir", default="./", type=str, help="directory to put template RMD file rnaseq_align_and_qc_report_Rmd_template.txt for QC report")
	parser.add_argument("project_name", type=str, help="Name of project that all samples correspond to.")
	parser.add_argument("samples_in", help="A tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.project_name, args.samples_in, args.path_start, args.aligner, args.template_dir)
	
	
