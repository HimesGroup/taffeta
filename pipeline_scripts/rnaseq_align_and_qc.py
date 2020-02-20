#!/usr/bin/python
import argparse
import sys
import subprocess
import os
import glob
import re
import rnaseq_userdefine_variables as userdef

def check_exist(file):
    if not os.path.exists(file):
        print "The file "+file+" does not exist"
        sys.exit()

def get_sample_info(fin):
    """
    Read in information from phenotype file
    Create a dictionary to store each column
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



def get_genome_ref_files(genome):
    """
    Location of all reference files needed for a given genome.
    The ERCC gtf files were appended separately to each species own gtf file
    Current choice: "hg38"
    """

    # read in user-defined variable python script: improt reference files

    if genome == "hg38":
	fa = userdef.hg38_fa
	gtf = userdef.hg38_gtf
	ref = userdef.hg38_ref
	ERCC_gtf = userdef.hg38_ERCC_gtf
	star_index_dir = userdef.hg38_star_index_dir

    elif genome == "hg19":
	fa = userdef.hg19_fa
	gtf = userdef.hg19_gtf
	ref = userdef.hg19_ref
	ERCC_gtf = userdef.hg19_ERCC_gtf
	star_index_dir = userdef.hg19_star_index_dir

    elif genome == "mm38":
	fa = userdef.mm38_fa
	gtf = userdef.mm38_gtf
	ref = userdef.mm38_ref
	ERCC_gtf = userdef.mm38_ERCC_gtf
	star_index_dir = userdef.mm38_star_index_dir

    elif genome == "mm10":
	fa = userdef.mm10_fa
	gtf = userdef.mm10_gtf
	ref = userdef.mm10_ref
	ERCC_gtf = userdef.mm10_ERCC_gtf
	star_index_dir = userdef.mm10_star_index_dir

    elif genome == "rn6":
	fa = userdef.rn6_fa
	gtf = userdef.rn6_gtf
	ref = userdef.rn6_ref
	ERCC_gtf = userdef.rn6_ERCC_gtf
	star_index_dir = userdef.rn6_star_index_dir

    elif genome == "susScr3":
	fa = userdef.susScr3_fa
	gtf = userdef.susScr3_gtf
	ref = userdef.susScr3_ref
	ERCC_gtf = userdef.susScr3_ERCC_gtf
	star_index_dir = userdef.susScr3_star_index_dir

    else:
	print 'Unknown genome selected: ', genome
	sys.exit()

    # check if reference files exist
    map(lambda x:check_exist, [fa, gtf, ref, ERCC_gtf, star_index_dir])

    return(fa, gtf, ref, ERCC_gtf, star_index_dir)

def get_adapter_info(fin, index_type):
    """
    Read in information from adapter sequence file
    Obtain adapter sequences of corresponding index type
    """
    f = open(fin,'r')
    f = f.read().split('\n')
    f = map(lambda x: x.rstrip(), f)

    if '' in f:
        f.remove('')

    header = f[0].split('\t') # list
    c = f[1:]

    # check if Type, Index, Description, and Sequence are in the file
    if set(["Type", "Index", "Description", "Sequence"]) != set(header):
        print "Not all the following colunms Type, Index, Description, and Sequence are in the file. Please check!"
        sys.exit()

    d = {} # key: column name, value: column

    for i in range(len(header)):
        colname=header[i]
        d[colname]=map(lambda x: x.split('\t')[i],c)

    # check if index_type within provided types
    if index_type not in d["Type"]:
        print index_type+" is not in the provided index types: "+', '.join(map(str,list(set(d["Type"]))))
        print "Please check!"
        sys.exit()

    # obtain idx of corresponding index type
    idx=filter(lambda x:index_type in d["Type"][x], range(len(c)))

    # Create dictionary to store index sequences of corresponding index type
    d_seq = {} # key: Index; value: Sequence

    for i in idx:
        index_type=d["Type"][i]
        index_seq=d["Index"][i]
        index_des=">"+d["Description"][i] # add '>' before description
        index_allseq=d["Sequence"][i]

        if index_seq not in d_seq:
            d_seq[index_seq]=index_des+"\n"+index_allseq
        else:
            d_seq[index_seq]=d_seq[index_seq]+"\n"+index_des+"\n"+index_allseq

    return d_seq

def index_check(indexes, index_type, template_dir):
    """
    For user specified index, check whether they are in the provided adapter sequence files.
    """	

    indexes=[x for x in indexes if x !='NA'] # remove NA values

    if index_type is None:
        print "Index type (--index_type) is not specified. Please check."
        sys.exit()

    # if unique dual (UD) index adapters are used, the two indexes i7 and i5 combined by "+" are needed
    if index_type in ["illumina_ud_sys1", "illumina_ud_sys1"]:
        if not all(map(lambda x:'+' in x, indexes)):
            print "Unique dual (UD) index type is specified. Provide i7 and i5 sequences using i7+i5 in Index column. Please check!"
            sys.exit()

    index_fn=template_dir+"rnaseq_adapter_primer_sequences.txt"
    print "index_type = "+index_type
    print "Use provided adapter and primer sequence file: "+index_fn
    check_exist(index_fn)

    # obtain adapter and primer sequences based on index type and index sequence
    index_dict=get_adapter_info(index_fn, index_type)
    index_seqs=index_dict.keys()

    # check if index sequences within specified index type
    # check if all indexes sequences from phenotype file are in provided reference file
    index_diff=list(set(indexes)-set(index_seqs))
    if len(index_diff) >0:
        print "The index(es) from phenotype file are not in the provided adapter sequence file: "+', '.join(map(str,index_diff))
        print "Please check!"
        sys.exit()

    return index_dict

def make_adapter_fa(curr_sample, out_dir, curr_index, index_dict):
    """
    Make a .fa file with adapter and primer sequences of corresponding sample
    """
    fa_out=index_dict[curr_index] # output index description and sequences from index dictionary
    primer_keys=filter(lambda x:"Primer" in x, index_dict.keys()) # obtain Read1 and/or Read2 primer names from dictionary
    for key in primer_keys:
        fa_out=fa_out+"\n"+index_dict[key]
    outp=open(out_dir+curr_sample+"_adapter.fa",'w') # output in curr_sample+"_adapter.fa"
    outp.write(fa_out)
    outp.write("\n")
    outp.close()

def name_fastqc(filename):
    """
    Create filename prefix of fastqc output
    """

    # extract basename
    filename=os.path.basename(filename)
    # remove extension
    filename=re.split("\.fastq\.|\.fastq$|\.txt$",filename)[0]
    return(filename)


def trim_and_fastqc(curr_sample, curr_index, out_dir, R1, R2, R1_trim, R2_trim):
    """
    Use Trimmomatic to trim adaptor and perform fastqc
    """
    
    # Trim adapter
    # read in user-defined variable python script: import trimmomatic java path
    trimmomatic=userdef.trimmomatic

    cmd="" # create command variable
    if curr_index == 'NA':
        print curr_sample+": no index defined. Skip trimming."
    else:
        if R2=="":
            if os.path.isfile(R1_trim):
                print curr_sample+" trimmed file already exists. Skip trimming."
            else:
                cmd = "java -Xmx1024m  -classpath "+trimmomatic+" org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n"

        else:
            R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"
            if os.path.isfile(R1_trim) and os.path.isfile(R1_trim):
                print curr_sample+" R1 and R2 trimmed files already exist. Skip trimming."
            else:
                cmd = "java -Xmx1024m  -classpath " +trimmomatic+" org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+R1+" "+R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:40\n"

    # Fastqc after trimming
    # create standard fastqc .zip file name
    R1_fastqc_fn = out_dir+curr_sample+"_R1_Trimmed_fastqc.zip"
    R2_fastqc_fn = out_dir+curr_sample+"_R2_Trimmed_fastqc.zip"

    if (R2=="" and os.path.isfile(R1_fastqc_fn)) or (os.path.isfile(R1_fastqc_fn) and os.path.isfile(R2_fastqc_fn)):
        print curr_sample+" fastqc results already exist. Skip fastqc."
    else:
        # Run fastqc
        cmd=cmd+"fastqc -o "+out_dir+" --extract "+R1_trim+" "+R2_trim+"\n" # since R2 could be "", no harm to add an empty string directly here
        if curr_index=='NA': # use no-trimmed files
            # retrieve original sample name without path and fastq extension
            R1_org_name=out_dir+name_fastqc(R1_trim)+"_fastqc" # change the path for original fastq file
            # rename fastqc result folder to current name
            cmd=cmd+"cd "+out_dir+"; mv "+R1_org_name+" "+curr_sample+"_R1_Trimmed_fastqc\n"
            # create zip file
            cmd=cmd+"zip -rm "+R1_fastqc_fn+" "+curr_sample+"_R1_Trimmed_fastqc\n"
            if R2_trim!="": # if R2 exists
                # retrieve original sample name without path and fastq extension
                R2_org_name=out_dir+name_fastqc(R2_trim)+"_fastqc" # change the path for original fastq file
                # rename fastqc result folder to current name
                cmd=cmd+"cd "+out_dir+"; mv "+R2_org_name+" "+curr_sample+"_R2_Trimmed_fastqc\n"
                # create zip file
                cmd=cmd+"zip -rm "+R2_fastqc_fn+" "+curr_sample+"_R2_Trimmed_fastqc\n"

    #Get total number of reads, unique reads, % unique reads from trimmed file(s).
    if ".gz" in R1_trim: # for gzip .fastq file
        cmd=cmd+"zcat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+out_dir+curr_sample+"_ReadCount\n"
    else:
        cmd=cmd+"cat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+out_dir+curr_sample+"_ReadCount\n"
    if R2_trim!="":
        if ".gz" in R2_trim:
            cmd=cmd+"zcat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+out_dir+curr_sample+"_ReadCount\n"
        else:
            cmd=cmd+"cat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+out_dir+curr_sample+"_ReadCount\n"

    return cmd

def star_align(curr_sample, out_dir, R1_trim, R2_trim, star_index_dir, strand):
    """
    Use STAR for alignment
    """

    cmd="mkdir "+out_dir+"star_out\n"
    cmd=cmd+"cd "+out_dir+"star_out\n"
    if R2_trim!="": # paired-end
        if strand!="nonstrand": # strand-specific assay
            cmd=cmd+"STAR --genomeDir "+star_index_dir+" --runThreadN 12 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim+" "+R2_trim
        else: # non-strand-specific assay
	    cmd=cmd+"STAR --genomeDir "+star_index_dir+" --runThreadN 12 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim+" "+R2_trim

    else: # single-end
        if strand!="nonstrand": # strand-specific assay
            cmd=cmd+"STAR --genomeDir "+star_index_dir+" --runThreadN 12 --sjdbOverhang 100 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim
        else: # non-strand-specific assay
            cmd=cmd+"STAR --genomeDir "+star_index_dir+" --runThreadN 12 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --readFilesIn "+R1_trim

    if ".gz" in R1_trim: # if .gz file used
        cmd=cmd+" --readFilesCommand zcat"
    cmd=cmd+"\n"

    return cmd


def get_bamstat_metrics(curr_sample, out_dir, ref, strand, library_type):
    """
    Obtain QC metrics from bam file
    """
    # bam stat
    # read in user-defined variable python script: import picard directory
    picard_dir=userdef.picard_dir

    #Create sorted bam file:
    cmd="samtools sort Aligned.sortedByCoord.out.bam -@12 -T "+curr_sample+".tmp -o "+curr_sample+"_accepted_hits.sorted.bam\n"
    #Create indexed bam file:
    cmd=cmd+"samtools index -@12 "+curr_sample+"_accepted_hits.sorted.bam\n"
    #Write out index stats of where reads align to by chr:
    cmd=cmd+"samtools idxstats "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.stats\n"
    #Write out bamtools summary stats:
    cmd=cmd+"bamtools stats -in "+curr_sample+"_accepted_hits.sorted.bam > "+curr_sample+"_accepted_hits.sorted.bamstats\n"

    #Run CollectRnaSeqMetrics
    if strand=="forward":
        cmd=cmd+"java -Xmx2g -jar "+picard_dir+"CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND VALIDATION_STRINGENCY=LENIENT INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n"
    elif strand=="reverse":
        cmd=cmd+"java -Xmx2g -jar "+picard_dir+"CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND VALIDATION_STRINGENCY=LENIENT INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n"
    elif strand=="nonstrand":
        cmd=cmd+"java -Xmx2g -jar "+picard_dir+"CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=NONE VALIDATION_STRINGENCY=LENIENT INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n"

    #Get number of reads spanning junctions by getting "N"s in CIGAR field of bam file. Create cigarN.script
    if not os.path.isfile(out_dir+"cigarN.script"):
        cigar = open(out_dir+"cigarN.script",'w')
        cigar.write("{\n")
        cigar.write("\"cigar\" : \"*N*\"\n")
        cigar.write("}\n")
        cigar.close()

    #Be sure that Junction Spanning Reads are Added First then Unmapped Reads for proper ordering of fields in report
    cmd=cmd+"echo \"Junction Spanning Reads: \" $(bamtools filter -in "+curr_sample+"_accepted_hits.sorted.bam -script "+out_dir+"cigarN.script | bamtools count ) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n"

    #Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
    if library_type in ["PE"]:
	cmd=cmd+"java -Xmx2g -jar "+picard_dir+"CollectInsertSizeMetrics.jar VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE="+curr_sample+"_InsertSizeHist.pdf INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_InsertSizeMetrics\n"

    return cmd

def htseq_quant(curr_sample, out_dir, ERCC_gtf, strand):
    """
    Use htseq-count for quantification
    """

    cmd = ""
    cmd=cmd+"mkdir "+out_dir+"htseq_out/\n"

    if strand=="nonstrand": # non-strand-specific assay: capture coding transcriptome without strand information
        cmd=cmd+"samtools view "+curr_sample+"_accepted_hits.sorted.bam | htseq-count -r pos --stranded=no - "+ERCC_gtf+" > "+out_dir+"htseq_out/"+curr_sample+"_counts.txt\n"
    elif strand=="reverse": # strand-specific assay: sequence first strand (reverse)
        cmd=cmd+"samtools view "+curr_sample+"_accepted_hits.sorted.bam | htseq-count -r pos --stranded=reverse - "+ERCC_gtf+" > "+out_dir+"htseq_out/"+curr_sample+"_counts.txt\n"
    elif strand=="forward": # strand-specific assay: sequence second strand (forward)
        cmd=cmd+"samtools view "+curr_sample+"_accepted_hits.sorted.bam | htseq-count -r pos --stranded=yes - "+ERCC_gtf+" > "+out_dir+"htseq_out/"+curr_sample+"_counts.txt\n"

    return cmd

def bam2bw_and_track(curr_sample, curr_color, out_dir, template_dir, track_fn, bigdata_path, len_fn):
    """
    Convert bam to bw file
    """

    cmd=""
    cmd=cmd+"genomeCoverageBed -split -bg -ibam "+curr_sample+"_accepted_hits.sorted.bam"+" -g "+len_fn+" > "+curr_sample+".bdg\n"
    cmd=cmd+"LC_COLLATE=C sort -k1,1 -k2,2n "+curr_sample+".bdg > "+curr_sample+".sorted.bdg\n"
    cmd=cmd+"bedGraphToBigWig "+curr_sample+".sorted.bdg "+len_fn+" "+curr_sample+".bw\n"

    """
    Create UCSC track file
    """

    outp=open(track_fn, 'a')
    outp.write("track type=bigWig name="+"\""+curr_sample+"\" color="+curr_color+" gridDefault=on maxHeightPixels=50 visibility=full autoScale=off viewLimits=5:100 description=\""+curr_sample+"\" bigDataUrl="+bigdata_path+curr_sample+".bw\n")
    outp.close()

    return cmd


def lsf_file(job_name, cmd, memory=36000, thread=12, queue=userdef.queue):
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


def main(sample_info_file, project_name, aligner, ref_genome, library_type, index_type, strand, path_start, template_dir, bam2bw):
    """
    Read in phenotype sample info provided by the users and perform the following steps:
    1) Perform adapter trimming - if Illumina indexes are not available in phenotype files
    2) run fastqc if they are not available
    3) Get unique reads 
    4) Run star to align reads to reference genome
    5) Run htseq to quantify mRNA, ERCC spike-in and rRNA (if applicable)
    6) Obtain various QC metrics on aligned files
    """

    ####
    # Set up and check
    ####

    # Set up project directory
    if path_start == "./":
        path_start = os.getcwd()
    if path_start[-1] != "/":
	path_start = path_start+"/"

    # Set up template directory
    if template_dir == "./":
        template_dir = os.getcwd()
    if template_dir[-1] != "/":
	template_dir = template_dir+"/"


    # Set up reference genome files
    fa, gtf, ref, ERCC_gtf, star_index_dir = get_genome_ref_files(ref_genome)
    print "ref_genome = "+ref_genome
    print "fa = "+fa
    print "gtf = "+gtf
    print "ref = "+ref
    print "ERCC_gtf = "+ERCC_gtf
    print "star_index_dir = "+ star_index_dir

    # Obtain sample information from phenotype file
    print "sample_info_file = " + sample_info_file
    check_exist(sample_info_file)
    info_dict = get_sample_info(sample_info_file)
    sample_names = info_dict["Sample"]

    # Check if R1 and/or R2 column for .fastq file path exists
    print "library_type = "+library_type

    if "R1" not in info_dict:
        print "R1 column for Read 1 .fastq paths does not exist. Please check!"
        sys.exit()
    if library_type in ["PE"]:
        if "R2" not in info_dict:
            print "Paired-end library is specified. R2 column for Read 2 .fastq paths does not exist. Please check!"
            sys.exit()

    # Print alignment and quantification parameters
    print "aligner = "+aligner
    print "strand = "+strand

    # Check if Index column exists and create index list
    # If index does not exist, assign NA
    if "Index" not in info_dict:
        print "No Index column in phenotype file. All samples skip adapter trimming."
        indexes = ['NA']*len(sample_names)
    else:
        indexes=info_dict['Index']

        if all(map(lambda x: 'NA' in x, indexes)): # if all indexes are NA
            print "All specified indexes are NA. Skip adapter trimming."

        # Obtain dictionary with adapter sequences of corresponding index type
        index_dict=index_check(indexes, index_type, template_dir)


    # If perform bam to bigwig file conversion, create url and colors
    if bam2bw:
        # overwrite ucsc track file if it already exists
        track_fn=path_start+project_name+"_ucsc_track.txt"
        if os.path.exists(track_fn):
            print "Warning: ucsc track file already exists: "+track_fn+". Overwrite it."
            outp=open(track_fn,'w')
            outp.write("")
            outp.close()

        # read in url
        bigdata_url=userdef.bigdata_url
        bigdata_path=bigdata_url+"/"+project_name+"/"
        print "Convert .bam files to .bw files. Use user-provided URL: "+bigdata_url
        print "Need to upload all the generated .bw files to this path: "+bigdata_path

        # obtain genome length file
        if ref_genome == "hg38":
    	    len_fn = userdef.hg38_len

        elif ref_genome == "hg19":
    	    len_fn = userdef.hg19_len

        check_exist(len_fn)

        # create color list for ucsc track display based on treatment condition
        rgb_colors=["27,158,119", "217,95,2", '117,112,179', '231,41,138', '102,166,30', '230,171,2', '166,118,29', '102,102,102'] # Dark2 color set "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666" convert to RGB (https://www.rapidtables.com/web/color/RGB_Color.html)
        if "Status" not in info_dict:
            colors=track_colors[0]*len(sample_names) # use one color for all samples
        else:
            conds=info_dict["Status"]
            cond_uniq=list(set(conds))
            # Assign color to each Status condition
            RGBS={}
            # repeat colors if the status is greater than the rgb_colors
            rgb_colors=(len(cond_uniq)/len(rgb_colors)+1)*rgb_colors
            for i in range(len(cond_uniq)):
                RGBS[cond_uniq[i]]=rgb_colors[i]
            colors=map(lambda x: RGBS[x], conds)

    ####
    # Run by each sample
    ####

    for i in range(len(sample_names)):
        curr_sample=sample_names[i]
        curr_index=indexes[i]

        # Obtain read1 and read2 .fastq file paths
        R1=info_dict["R1"][i]
        check_exist(R1)
        if library_type in ["PE"]:
            R2=info_dict["R2"][i]
            check_exist(R2)
        elif library_type in ["SE"]:
            R2=""

        # Create output directory
        out_dir = path_start+curr_sample+"/"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Create cmd varialble to save linux commands output in .lsf files
        cmd = "cd "+out_dir+"\n"

        ###
        # Trim Adaptor and fastqc
        ###

        if curr_index=="NA": # no index specified. skip trim
            R1_trim = R1
            R2_trim = R2
        else:
            # create adapter .fa file
            make_adapter_fa(curr_sample, out_dir, curr_index, index_dict)
            # create trimmed filename
            R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
            if R2=="":
                R2_trim = ""
            else:
                R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"
            
        # trim_and_fastqcL 1. trim, 2. fastqc, 3 total unique counts from .fastq file
        trim_and_fastqc_cmd=trim_and_fastqc(curr_sample, curr_index, out_dir, R1, R2, R1_trim, R2_trim)
        cmd = cmd + trim_and_fastqc_cmd

        ###
        # Alignment
        ###

        # STAR alignment
        if aligner == "star": 
            star_align_cmd = star_align(curr_sample, out_dir, R1_trim, R2_trim, star_index_dir, strand)
            cmd = cmd + star_align_cmd

        ###
        # Obtain QC metrics from .bam file
        ###
        bamstat_cmd = get_bamstat_metrics(curr_sample, out_dir, ref, strand, library_type)
        cmd = cmd + bamstat_cmd

        ###
        # Quantification
        ###
        
        htseq_quant_cmd = htseq_quant(curr_sample, out_dir, ERCC_gtf, strand)
        cmd = cmd + htseq_quant_cmd

        ###
        # Convert bam to bw
        ###

        if bam2bw:
            curr_color=colors[i]
            bam2bw_cmd=bam2bw_and_track(curr_sample, curr_color, out_dir, template_dir, track_fn, bigdata_path, len_fn)
            cmd=cmd+bam2bw_cmd


        ###
        # Create .lsf files
        ###
        lsf_file(curr_sample+"_align", cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a project.")
    parser.add_argument("--project_name", type=str, help="Prefix name of for output directory and files")
    parser.add_argument("--samples_in", help="A tab-delimited txt file containing sample information with full path. See example file: sample_info_file.txt")
    parser.add_argument("--aligner", default="star", type=str, help="Aligner name. We used to have tophat which is retired now."
	"(options: star)")
    parser.add_argument("--ref_genome", default="hg38", type=str, help="Specify reference genome (options: hg38, hg19, mm38, mm10, rn6, susScr3)")
    parser.add_argument("--library_type", default="PE", type=str, help="Specify library type (options: PE (paired-end), SE (single-end))")
    parser.add_argument("--index_type", type=str, help="If Index column is in phenotype file, specify index type for adapter trim."
        "(options: truseq_single_index, illumin_ud_sys1, illumin_ud_sys2 or user specified in the user-defined adapter reference file.)")
    parser.add_argument("--strand", type=str, default="nonstrand", help="Whether data is from a strand-specific assay."
        "(options: nonstrand, reverse, forward)")
    parser.add_argument("--path_start", default="./", type=str, help="Directory path where project-level directories are located and report directory will be written (default=./)")
    parser.add_argument("--template_dir", default="./", type=str, help="directory to put provided or user defined reference index files")
    parser.add_argument("--bam2bw", action='store_true', help="If specified, generate bigwig files (.bw) and create ucsc track file.")
    args = parser.parse_args()

    if args.project_name is None or args.samples_in is None:
        parser.print_help()
        sys.exit()

    main(args.samples_in, args.project_name, args.aligner, args.ref_genome, args.library_type, args.index_type, args.strand, args.path_start, args.template_dir, args.bam2bw)

