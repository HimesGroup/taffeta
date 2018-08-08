#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import re
import fnmatch

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


def lsf_file(job_name, cmd, memory=36000, thread=1):
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
    outp.write("#BSUB -n "+str(thread)+"\n")
    outp.write(cmd)
    outp.write("\n")
    outp.close()


def make_deseq2_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file):
    """
    Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
    """

    ###
    # Create global design variable based on comparison file
    ###

    #load text file containing all comparisons of interest
    comps_file = open(comp_file)
    comps = comps_file.readlines()[1:] # exclude header line

    # obtain paired variable -- control_vars
    if "paired:" in ''.join(comps): # if paired variable exists in comparison file
        control_vars=map(lambda x:x.rstrip().split(':')[1],filter(lambda x:"paired:" in x, comps)) # obtain all paired variables
        control_vars=list(set(control_vars)) # obtain unique variables

        # check if paired variables are within sample_info_file
        samps_file=open(sample_info_file)
        samps = (samps_file.readlines()[0].rstrip()).split('\t') # obtain column names

        if len(list(set(control_vars)-set(samps)))>0:
            print "The paired variable(s) must be columns in the sample_info_file "+sample_info_file+". The following specified variables are not in sample_info_file: "+" ,".join(list(set(control_vars)-set(samps)))
            sys.exit()

    # assign global variable

    if "unpaired" in ''.join(comps) and "paired:" not in ''.join(comps):
        design="unpaired"
    elif "unpaired" not in ''.join(comps) and "paired:" in ''.join(comps):
        design="paired"
    else:
        design="mixed"

    # Create out directory
    out_dir = path_start+project_name+"_deseq2_out/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ###
    # Create RMD file for DESeq2
    ###

    # create custom.css for rmarkdown
    make_rmd_css(out_dir)

    outp = open(out_dir+'/'+project_name+"_DESeq2_Report.Rmd", "w")

    import rnaseq_userdefine_variables as userdef # read in user-defined variable python script
    # import software version
    star_version=userdef.star_version
    htseq_version=userdef.htseq_version
    deseq2_version=userdef.deseq2_version
    # import author information
    author=userdef.author
    # import favorite genes
    fav_gene=userdef.fav_gene

    # title
    outp.write("---\ntitle: 'Differential Expression Results for "+project_name+"'\n")
    outp.write("author: "+author+"\n")
    outp.write("date: \"`r format(Sys.time(), '%d %B, %Y')`\"\n")
    outp.write("output: \n")
    outp.write("  html_document:\n")
    outp.write("    css: custom.css\n")
    outp.write("    toc: true\n")
    outp.write("    toc_float: true\n---\n\n")

    # description
    outp.write("Reads were aligned to the "+ref_genome+" assembly using STAR ("+star_version+").  The following alignment QC report was produced:<br>\n\n")
    outp.write("> "+project_name+"_QC_RnaSeqReport.html<br>\n\n")
    outp.write("HTSeq ("+htseq_version+") function htseq-count was used to count reads. Counts for all samples were concatenated into the following text file:<br>\n\n")
    outp.write("> "+project_name+"_htseq_gene.txt<br>\n\n")
    outp.write("DESeq2 ("+deseq2_version+") was used for differential gene expression analaysis, based on the HTSeq counts matrix and the phenotype file provided.  Normalized counts from DESeq2 are saved in the following text file:<br>\n\n")
    outp.write("> "+project_name+"_counts_normalized_by_DESeq2.txt<br>\n\n")
    outp.write("Normalized counts are obtained from DESeq2 function estimateSizeFactors(), which divides counts by the geometric mean across samples; this function does not correct for read length. The normalization method is described in detail here: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106<br>\n\n")
    outp.write("Differential gene expression analysis was done for all comparisons provided in the comparisons file.  The following design was used:<br>\n\n")

    if design=="unpaired":
        outp.write("> design = ~ Status<br>\n\n")
    elif design=="paired":
        for control_var in control_vars:
            outp.write("> design = ~ "+" + "+control_var+" + Status<br>\n\n")
    elif design=="mixed":
        outp.write("For unpaired comparisons:\n")
        outp.write("> design = ~ Status<br>\n\n")
        outp.write("For paired comparisons:\n")
        for control_var in control_vars:
            outp.write("> design = ~ "+" + "+control_var+" + Status<br>\n\n")

    outp.write("If desired, the design can be modified to include more independent variables. In addition to the partial results displayed in this report, the full set of DESeq2 results for each comparison was saved down in separate text files, with names of the form:<br>\n\n")
    outp.write("> "+project_name+"_CASE_vs_CONTROL_DESeq2_results.txt<br>\n\n")
    outp.write("where CASE and CONTROL are pairs of conditions specified in the comparisons file.<br>\n\n")
    outp.write("\n\n```{r lib, echo=F, message=F, warning=F}\n")
    outp.write("library(gplots)\nlibrary(reshape2)\nlibrary(RColorBrewer)\nlibrary(plyr)\nlibrary(lattice)\nlibrary(genefilter)\nlibrary(ggplot2)\nlibrary(viridis)\nlibrary(DESeq2)\nlibrary(DT)\nlibrary(tidyr)\nlibrary(biomaRt)\nlibrary(pander)\noptions(width = 1000)\n```\n")
    outp.write("\n")
    outp.write("\n\n```{r vars, eval=T, echo=F}\n")
    outp.write("project_name=\""+project_name+"\"\n")
    outp.write("path.start='"+path_start+"'\n")
    if ref_genome=="hg38":
        outp.write("housekeeping_genes <- c('ACTB','GAPDH','B2M','RPL19','GABARAP')\n")
        outp.write("house_list=list() # create an empty list to save results of house-keeping genes\n")
    elif ref_genome=="mm38" or ref_genome=="mm10" or ref_genome=="rn6":
	outp.write("housekeeping_genes <- c('Actb','Gapdh','B2m','Rpl19','Gabarap')\n") #also had added 'Hprt','Gnb2l1','Rpl32' for mouse_macrophages
        outp.write("house_list=list() # create an empty list to save results of house-keeping genes\n")
    else:
	print("Housekeeping genes only available for human, mouse and rat.")

    # user-defined favorite genes
    outp.write("fav_genes <- c('"+"', '".join(fav_gene)+"')\n")

    #Read in phenofile; phenofile should contain columns"Sample" and "Status"; if desired, it can have other columns (e.g. "Individual" or "Cell_Line", etc.) & design can be modified accordingly
    outp.write("coldata <- read.table('"+sample_info_file+"', sep='\\t', header=TRUE)\n")
    outp.write("coldata <- subset(coldata, QC_Pass==1)\n")

    #load counts from HTSeq output
    outp.write("countdata <- read.table(paste0(path.start, project_name,'_Alignment_QC_Report_star/', project_name, '_htseq_gene.txt'), sep='\\t', header=TRUE, check.names=FALSE)\n") # need check.names=FALSE -- else dashes in colnames converted to periods
    outp.write("names(countdata) <- gsub('count_','',names(countdata))\n") # remove "count_" from column names
    outp.write("countdata <- countdata[,c('Gene',as.character(coldata$Sample))]\n\n") # subset counts with samples who passed QC
    outp.write("row.names(countdata) <- countdata$Gene\n")
    outp.write("is_ensg <- if (substr(rownames(countdata)[1], 1, 4) == 'ENSG') {TRUE} else {FALSE}\n")
    outp.write("countdata$Gene <- sapply(strsplit(as.character(countdata$Gene), '\\\.'), '[[', 1) # remove .5 in transcript ENSG00000000005.5\n")

    if ref_genome=="hg38":
        outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='hsapiens_gene_ensembl')\n")
	outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_gene_id', 'hgnc_symbol'), values=countdata$Gene, mart=mart) # 30 duplicates transcripts from ensembles. They are microRNA and SRP RNAs. Just use the unique transcripts.\n")
        outp.write("genes <- genes[!duplicated(genes$ensembl_gene_id),]\n")
	outp.write("if (is_ensg) {countdata <- merge(countdata, genes, by.x='Gene', by.y='ensembl_gene_id')} else {countdata <- countdata}\n")
	outp.write("if (is_ensg) {countdata <- rename(countdata, c('hgnc_symbol'='gene_symbol'))} else {countdata$gene_symbol <- countdata$Gene}\n")

    elif ref_genome=="mm38" or ref_genome=="mm10":
	outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='mmusculus_gene_ensembl')\n")
	outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_gene_id', 'mgi_symbol'), values=countdata$Gene, mart=mart)\n")
        outp.write("genes <- genes[!duplicated(genes$ensembl_gene_id),]\n")
	outp.write("if (is_ensg) {countdata <- merge(countdata, genes, by.x='Gene', by.y='ensembl_gene_id')} else {countdata <- countdata}\n")
	outp.write("if (is_ensg) {countdata <- rename(countdata, c('mgi_symbol'='gene_symbol'))} else {countdata$gene_symbol <- countdata$Gene}\n")

    elif ref_genome=="rn6":
	outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='rnorvegicus_gene_ensembl')\n")
	outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_gene_id', 'rgd_symbol'), values=countdata$Gene, mart=mart)\n")
        outp.write("genes <- genes[!duplicated(genes$ensembl_gene_id),]\n")
	outp.write("if (is_ensg) {countdata <- merge(countdata, genes, by.x='Gene', by.y='ensembl_gene_id')} else {countdata <- countdata}\n")
	outp.write("if (is_ensg) {countdata <- rename(countdata, c('rgd_symbol'='gene_symbol'))} else {countdata$gene_symbol <- countdata$Gene}\n")

    else:
	print("This code can only append official gene symbols for hg38, mm38, rn6.")

    outp.write("if (is_ensg) {row.names(countdata) <- countdata$Gene} else {row.names(countdata) <- row.names(countdata)}\n\n")
    outp.write("```\n")

    ###
    # DESeq2 analysis
    ###

    outp.write("\n\n```{r normcount, eval=T, echo=F, message=F, warnings=F, results=\"hide\"}\n")
    #this part is only for getting normalized counts (this way all samples normalized together)
    outp.write("\n# this part is only for getting normalized counts (this way all samples normalized together)\n")


    outp.write("# unpaired design, testing for effect of Status\n")
    outp.write("ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata[,2:(ncol(countdata)-1)], colData = coldata, design = ~ Status)\n")

    outp.write("dds <- DESeq(ddsFullCountTable)\n") 
    # not actually doing analysis here, so reference level does not matter.
    outp.write("norm.counts <- counts(dds, normalized=TRUE)\n")
    outp.write("norm.counts <- merge(norm.counts, countdata[,c('Gene','gene_symbol')], by='row.names')\n")
    outp.write("norm.counts <- norm.counts[,2:ncol(norm.counts)] #else end up with a column called 'Row.names'\n")
    outp.write("```\n\n")
    outp.write("```{r normcount_save, eval=T, echo=F}\n")
    outp.write("write.table(norm.counts, paste0(project_name,'_counts_normalized_by_DESeq2.txt'), sep='\\t', quote=F, row.names=F, col.names=T)\n")
    outp.write("```\n\n")

    #create and paste the portion of the report that is unique to each comparison
    for line in comps:
        line=line.rstrip()
        case=line.split('\t')[0]
        ctrl=line.split('\t')[1]

        # define local design variable for each comparison
        if "unpaired" in line.split('\t')[2]:
            design="unpaired"
        elif "paired:" in line.split('\t')[2]:
            design="paired"
            control_var=(line.split('\t')[2]).split(':')[1]

        outp.write("```{r, eval=T, echo=F}\n")
        outp.write("case <- '"+case+"'\n")
        outp.write("ctrl <- '"+ctrl+"'\n")
        #outp.write("res <- results(dds, contrast=c('Status','"+case+"','"+ctrl+"'))\n")
        outp.write("```\n\n")
        outp.write("## "+case+" vs. "+ctrl+"\n")
        outp.write("\n")
        outp.write("### Samples in this comparison\n")
        outp.write("```{r, eval=T, echo=F, message=F}\n")
        outp.write("#conditions from file - select the portion of the info sheet relevant to the two conditions being tested\n")
        outp.write("coldata_curr <- coldata[which(coldata$Status==case | coldata$Status==ctrl),]\n")


        if design=="paired": # if there is a sample without a pair in this comparison, remove it from the comparison
            outp.write("#for paired design, want each unique value of control variable to correspond to an equal number of case and control samples -- else toss all samples with that value of control variable.\n")
	    outp.write("for (i in unique(coldata_curr$"+control_var+")) {\n")
	    outp.write("	if (nrow(coldata_curr[which(coldata_curr$Status=='"+case+"' & coldata_curr$"+control_var+"==i),]) != nrow(coldata_curr[which(coldata_curr$Status=='"+ctrl+"' & coldata_curr$"+control_var+"==i),])) {\n")
	    outp.write("		coldata_curr <- coldata_curr[-which(coldata_curr$"+control_var+"==i),]\n")
	    outp.write("		cat('Samples with "+control_var+"==', i, 'were removed from the "+case+" vs. "+ctrl+" comparison because there was an unequal number of case and control samples')\n}\n}\n")


        outp.write("coldata_curr <- coldata_curr[order(as.character(coldata_curr$Sample)), ]\n")
        outp.write("coldata_curr$Status <- factor(coldata_curr$Status, levels = c(ctrl,case)) # make sure that control is being used as reference level (else DESeq2 does it alphabetically)\n")
        outp.write("DT::datatable(coldata_curr, rownames=FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all')))) # dom = 't' removes search box\n")
        outp.write("rownames(coldata_curr) <- coldata_curr$Sample\n")
        outp.write("coldata_curr$Sample <- NULL\n")
        outp.write("```\n\n")

        outp.write("### DE analysis\n")
        outp.write("```{r, eval=T, echo=F, message=F}\n")
        outp.write("#select the portion of the HTSeq output matrix relevant to the two conditions being tested\n")
        outp.write("a <- colnames(countdata)\n")
        outp.write("b <- rownames(coldata_curr)\n")
        outp.write("countdata_curr <- countdata[,which(a %in% b)] # subset countdata and coldata\n")
        outp.write("countdata_curr <- countdata_curr[ ,order(as.character(names(countdata_curr)))] #columns of countdata & rows of coldata must be ordered in the same way\n")
        outp.write("#pre-filter low count genes before running the DESeq2 functions. Keep only genes (rows) that have at least 10 reads total.\n")
        outp.write("keep <- rowSums(countdata_curr)>=10\n")
        outp.write("countdata_curr <- countdata_curr[keep,]\n")
        outp.write("#combine HTSeq counts and info from info sheet\n")
        outp.write("#in spedifying design, order matters: test for the effect of condition (the last factor), controlling for the effect of individual (first factor)\n")

        if design=="unpaired":
            outp.write("# unpaired design, testing for effect of Status\n")
	    outp.write("ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata_curr, colData = coldata_curr, design = ~ Status)\n")
        elif design=="paired":
	    outp.write("# paired design, testing for effect of Status while controlling for "+control_var+"\n")
	    outp.write("ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata_curr, colData = coldata_curr, design = ~ "+control_var+" + Status)\n")

	outp.writelines(rmd_template)
	outp.write("\n")

    #Housekeeping gene expression barplots - once per report
    outp.write("## Housekeeping genes\n")
    outp.write("Counts have been normalized by estimated size factors using DESeq2. Obtain the count matrix using function DESeq2::counts.\n\n")
    outp.write("The table shows p-values of house-keeping genes for each comparison. Generally, house-keeping gene expressions do not change significantly in different conditions.\n")
    outp.write("\n")
    outp.write("```{r, eval=T, echo=F, cache=F, warning=F, message=F}\n")
    outp.write("if (exists(\"housekeeping_genes\")) {\n")
    outp.write("  sel_ids=as.character() # save selected Ensembl ID for house-keeping genes\n")
    outp.write("  for (i in housekeeping_genes) {\n")
    outp.write("    gene_symbol <- i\n")
    outp.write("    gene_ids <- norm.counts[which(norm.counts$gene_symbol==i),'Gene']\n")
    outp.write("    curr_data <- norm.counts[which(norm.counts$gene_symbol==i),names(norm.counts)[names(norm.counts)%in%coldata$Sample]]\n")
    outp.write("    curr_data <- curr_data[which.max(rowSums(curr_data)),]\n")
    outp.write("    gene_id <- gene_ids[which.max(rowSums(curr_data))] # choose Ensembl gene with the most counts\n")
    outp.write("    sel_ids <- c(sel_ids, gene_id)\n")
    outp.write("    curr_data <- data.frame(Sample=colnames(curr_data),Gene=rep(gene_id,ncol(curr_data)),gene_symbol=rep(gene_symbol,ncol(curr_data)),count=as.numeric(curr_data))\n")
    outp.write("    curr_data <- merge(curr_data,coldata[which(coldata$Sample%in%curr_data$Sample),c('Sample','Status')],by='Sample')\n")
    outp.write("    print(boxplot_func(df=curr_data))")
    outp.write("  }\n")
    outp.write("  dat <- do.call(rbind,lapply(1:length(house_list),function(x){dat=house_list[[x]];dat$Comparison=names(house_list)[x];dat[which(dat$Gene%in%sel_ids),c('gene_symbol','pvalue','Comparison')]})) # obtain all p-value for house-keeping genes\n")
    outp.write("  dat[,c('pvalue')] <- formatC(dat[,c('pvalue')], format = \"e\", digits = 2)\n")
    outp.write("  dat <- dat %>% spread(gene_symbol, pvalue)\n")
    outp.write("  DT::datatable(dat, rownames=FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = \"_all\"))))\n")    
    outp.write("}\n```\n\n")

    outp.write("```{r session_info, eval=T, echo=F}\n")
    outp.write("pander(sessionInfo())\n")
    outp.write("```\n\n")
    outp.close()

    # create .lsf file for HPC use
    lsf_cmd="cd "+out_dir+"; echo \"library(rmarkdown); rmarkdown::render('"+project_name+"_DESeq2_Report.Rmd')\" | R --no-save --no-restore\n"
    lsf_file(project_name+"_deseq2", lsf_cmd)

	
def make_sleuth_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file):
	"""
	Creates Rmd report. The top of report is below and the rest concatenated from a separate text document (rmd_template).
	"""
	out_dir = path_start+"/"+project_name+"/"+project_name+"_DE_Report/"
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	sleuth_dir = path_start+"/"+project_name+"/sleuth_out/"

        # create custom.css for rmarkdown
        make_rmd_css(out_dir)

	outp = open(out_dir+project_name+"_Sleuth_Report.Rmd", "w")
	outp.write("---\ntitle: \"Sleuth results - based on Kallisto TPM\"\n")
	outp.write("output: html_document \ntoc: true \ntoc_depth: 2 \n---\n")
	outp.write("\n")
	outp.write("## Sleuth results - based on Kallisto TPM for "+project_name+"\n\n")
	outp.write("Kallisto was used to quantify transcritpt abundances.  Kallisto TPMs for each comparison provided in the comparisons file were saved down in separate text files, with names of the form: <blockquote>"+project_name+"_CASE_vs_CTRL_kallisto_output.txt</blockquote>  Note that each comparison's kallisto output file only contains kallisto results for samples/conditions relevant to that comparison.\n\n")
	outp.write("Sleuth was used to compute transcript-level differential expression results based on kallisto TPMs.  Sleuth results for each comparison provided in the comparisons file were saved down in separate text files, with names of the form: <blockquote>"+project_name+"_CASE_vs_CTRL_sleuth_output.txt</blockquote>\n")
	outp.write("\nProject "+project_name+" consists of the following samples:\n")
	outp.write("```{r set-options, echo=F, cache=F, warning=F, message=F}\n")
	outp.write("library(reshape2)\nlibrary(sleuth)\nlibrary(biomaRt)\nlibrary(dplyr)\noptions(width = 2000)\n")
	outp.write("\n")
	outp.write("curr.batch=\""+project_name+"\"\n")
	outp.write("path.start='"+path_start+"'\n")	
	outp.write("kallisto_output <- data.frame()\n")
	if ref_genome=="hg38":
		outp.write("housekeeping_genes <- c('ACTB','GAPDH','B2M','RPL19','GABARAP')\n")
                outp.write("house_list=list() # create an empty list to save results of house-keeping genes\n")
	elif ref_genome=="mm38" or ref_genome=="mm10" or ref_genome=="rn6":
		outp.write("housekeeping_genes <- c('Actb','Gapdh','B2m','Rpl19','Gabarap')\n")
                outp.write("house_list=list() # create an empty list to save results of house-keeping genes\n")
	else:
		print("Housekeeping genes only available for hg38, mm38, rn6.")

	#load info sheet 
	outp.write("info_sheet <- read.table('"+sample_info_file+"', header = TRUE, stringsAsFactors=FALSE)\n") #outp.write("info_sheet <- read.table(paste0(path.start,'"+sample_info_file+"'), header = TRUE, stringsAsFactors=FALSE)\n")
	#outp.write("info_sheet <- subset(info_sheet, select=c('Sample','Status'))\n") #regular
	outp.write("info_sheet <- subset(info_sheet, select=c('Sample','Status','Day'))\n") #paired
	outp.write("colnames(info_sheet) <- c('run_accession', 'condition')\n")
	outp.write("info_sheet <- info_sheet[order(info_sheet$condition),]\n")
        outp.write("info_sheet <- subset(info_sheet, QC_Pass==1)\n")
	outp.write("print(info_sheet, row.names=FALSE)\n")
	outp.write("\n```\n")

	#load text file containing all comparisons of interest
	comps_file = open(comp_file)
	comps = comps_file.readlines()

	#create and paste the portion of the report that is unique to each comparison
	for str in comps:
		#note sleuth needs case and control to be in alphabetical order
		#i.e. the thing that comes first alphabetically will become the control
		#the first string in comps_file is the one you want to be 'case'; hence, add 'zz_' to it

		str1 = 'zz_' + re.split('_vs_', str)[0]
		str2 = re.split('_vs_', str)[1]
		mylist = [str1, str2]
		mylist.sort()
		ctrl = mylist[0].rstrip() #control will be the one that comes first alphabetically
		case = mylist[1].rstrip() #case wil be the one that starts with 'zz_' and hence is second in the sorted list
		outp.write("```{r, echo=F}\n")
		outp.write("case <- paste0('condition','"+case+"')\n")
		outp.write("ctrl <- paste0('condition','"+ctrl+"')\n")
		outp.write("\n```\n")
		outp.write("### "+case+" vs. "+ctrl+" comparison\n")
		outp.write("```{r, echo=F}\n")
		outp.write("so <- readRDS(paste0('"+sleuth_dir+"', 'so_', '"+case+"', '_vs_', '"+ctrl+"', '.rds'))\n")
		outp.write("k_curr <- kallisto_table(so)\n")
		outp.write("k_transc <- k_curr$target_id #save down before string split b/c sleuth object has transcripts in original form\n") 
		outp.write("k_curr$target_id <- sapply(strsplit(as.character(k_curr$target_id), '\\\.'), '[[', 1)\n")

		if ref_genome=="hg38":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='hsapiens_gene_ensembl')\n")
			outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), values=k_curr$target_id, mart=mart)\n")
			outp.write("genes <- dplyr::rename(genes, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)\n")
		elif ref_genome=="mm38" or ref_genome=="mm10":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='mmusculus_gene_ensembl')\n")
			outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_transcript_id', 'ensembl_gene_id', 'mgi_symbol'), values=k_curr$target_id, mart=mart)\n")
			outp.write("genes <- dplyr::rename(genes, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = mgi_symbol)\n")
		elif ref_genome=="rn6":
			outp.write("mart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', host='mar2016.archive.ensembl.org', path='/biomart/martservice' ,dataset='rnorvegicus_gene_ensembl')\n")
			outp.write("genes <- biomaRt::getBM(attribute=c('ensembl_transcript_id', 'ensembl_gene_id', 'rgd_symbol'), values=k_curr$target_id, mart=mart)\n")
			outp.write("genes <- dplyr::rename(genes, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = rgd_symbol)\n")
		else:
			print("This code can only append official gene symbols for hg38, mm38, mm10 or rn6.")
		outp.write("k_curr <- merge(k_curr, genes, by.x='target_id', by.y='target_id')\n")		
		outp.write("kallisto_output <- rbind(kallisto_output, k_curr)\n")
		outp.write("kallisto_output$sample <- gsub('zz_','',as.character(kallisto_output$sample))\n") 
		outp.write("kallisto_output$condition <- gsub('zz_','',as.character(kallisto_output$condition))\n") 
		outp.write("kallisto_output <- unique(kallisto_output)\n")
		outp.write("```\n")
		outp.write("\n")
		outp.writelines(rmd_template)
		outp.write("\n")

	#Housekeeping gene expression barplots - once per report
	outp.write("### Housekeeping genes\n")
	outp.write("```{r, eval=T, echo=F}\n")
	outp.write("house_ensg <- unique(kallisto_output$ens_gene[which(kallisto_output$ext_gene %in% housekeeping_genes)])\n")
	outp.write("for (i in 1:length(house_ensg)) {\n")
  	outp.write("  curr_gene <- house_ensg[i]\n")
  	outp.write("  curr_data <- kallisto_output[which(kallisto_output$ens_gene==curr_gene),]\n")
 	outp.write("  curr_data <- subset(curr_data, select=c(target_id, tpm, condition, ext_gene))\n")
	outp.write("  curr_data$condition <- factor(curr_data$condition)\n") #, levels=c('Control_Baseline', 'Asthma_Baseline'))\n")
  	outp.write("  gene_symbol <- unique(curr_data$ext_gene)\n")
	outp.write("  print({\n")
    	outp.write("	ggplot(curr_data, aes(x = condition, y = tpm, fill=condition)) + \n")
        outp.write("	geom_boxplot(outlier.colour=NA, lwd=0.2, color='grey18') + \n")
       	outp.write("	stat_boxplot(geom ='errorbar', color='grey18') + \n")
       	outp.write("	geom_jitter(size=0.8, width=0.2) + \n")
      	outp.write("	facet_wrap(~target_id) + \n")
       	outp.write("	guides(fill=FALSE) + \n")
       	outp.write("	theme_bw() +  \n")
       	outp.write("	labs(title=gene_symbol) + \n")
       	outp.write("	labs(x='condition') + labs(y='TPM') + \n")
       	outp.write("	theme(text = element_text(size=9), \n")
        outp.write("	  strip.text.x = element_text(size = 10), \n")
        outp.write("	  axis.text.x = element_text(angle = 90, hjust = 1, size=12),\n")
        outp.write("	  axis.text.y = element_text(size=9),\n")
        outp.write("	  title = element_text(size=12),\n")
        outp.write("	  axis.title.x = element_text(size=12),\n")
        outp.write("	  axis.title.y = element_text(size=12))})\n")
	outp.write("}\n")
	outp.write("```\n")
	outp.close()
	#subprocess.call("cd "+out_dir+"; echo \"library(knitr); library(markdown); knit2html('"+project_name+"_Sleuth_Report.Rmd', force_v1 = TRUE, options = c('toc', markdown::markdownHTMLOptions(TRUE)))\" | R --no-save --no-restore", shell=True)

        # create .lsf file for HPC use
        lsf_cmd="cd "+out_dir+"; echo \"library(rmarkdown); rmarkdown::render('"+project_name+"_Sleuth_Report.Rmd')\" | R --no-save --no-restore\n"
        lsf_file(project_name+"_sleuth", lsf_cmd)

def main(project_name, sample_info_file, de_package, path_start, comp_file, template_dir, ref_genome):
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
        path_start = path_start+"/"

    if template_dir == "./":
	template_dir = os.getcwd()
    if template_dir[-1] != "/":
        template_dir = template_dir+"/"

    # check if sample info exists
    if not os.path.exists(sample_info_file):
        print "Cannot find sample_info_file: "+sample_info_file
        sys.exit()

    # check if compare file exists
    if not os.path.exists(comp_file):
        print "Cannot find the comparison file: "+comp_file
        sys.exit()

    if de_package == "cummerbund":
        # check if rnw template file exists
        if not os.path.exists(template_dir+"rnaseq_de_report_Rnw_template.txt"):
            print "Cannot find rnaseq_de_report_Rnw_template.txt"
	    sys.exit()
	
        rnw_in = open(template_dir+"rnaseq_de_report_Rnw_template.txt", "r")
	rnw_template = rnw_in.read()
	make_de_rnw_html(rnw_template, project_name, path_start, gtf, ref_genome)

    elif de_package == "deseq2":
        # check if count matrix file exists
        if not os.path.exists(path_start+project_name+"_Alignment_QC_Report_star/"+project_name+"_htseq_gene.txt"):
            print "Cannot find gene count matrix file: "+path_start+project_name+"_Alignment_QC_Report_star/"+project_name+"_htseq_gene.txt"
            sys.exit()

	# check if deseq2 template file exists
	if not os.path.exists(template_dir+"rnaseq_deseq2_Rmd_template.txt"):
	    print "Cannot find rnaseq_deseq2_Rmd_template.txt"
	    sys.exit()

        rmd_in = open(template_dir+"rnaseq_deseq2_Rmd_template.txt", "r")
        rmd_template = rmd_in.readlines()
	make_deseq2_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file)

    if de_package == "sleuth":
	# check if sleuth template file exists
        if not os.path.exists(template_dir+"rnaseq_sleuth_Rmd_template.txt"):
            print "Cannot find rnaseq_sleuth_Rmd_template.txt"
	    sys.exit()

        rmd_in = open(template_dir+"rnaseq_sleuth_Rmd_template.txt", "r")
	rmd_template = rmd_in.readlines()
	make_sleuth_html(rmd_template, project_name, path_start, sample_info_file, ref_genome, comp_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create HTML report of differential expression results for RNA-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of project for all output files.")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt, but add an additional QC_Pass column")
    parser.add_argument("--comp", help="A tab-delimited txt file containing sample comparisons to be made. One comparison per line, columns are Condition1, Condition0, Design. "
            "Design: specify paired or unpaired. For paired design, specify condition to correct for, matching the column name in the 'coldata' file - e.g. paired:Donor.")
    parser.add_argument("--de_package", default="deseq2", type=str, help="Should be DESeq2 or sleuth be used for diferential expression (DE) analysis? If sleuth, a larger memory ~36 Mb is required."
            "(options: deseq2, sleuth)")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome (options: hg38, mm38, mm10, rn6)")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path to where project directory created by rnaseq_de.py is located (default=./)")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put template RMD file rnaseq_sleuth_Rmd_template.txt for QC report")

    args = parser.parse_args()

    if args.comp is None or args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.project_name, args.samples_in, args.de_package, args.path_start, args.comp, args.template_dir, args.ref_genome)


