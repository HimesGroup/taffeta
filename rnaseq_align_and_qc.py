#!/usr/bin/python
import subprocess
import os
import argparse


def get_sample_info(fin):
	"""
	Open tab-delimited txt file containing the following columns:
	v0: customer_ID		| ID given to sample by customer
	v1: gigpad_ID		| A sample ID that may differ from "Customer_ID" and that is associated with library prior to pooling
	v2: lane			| Lane of sequencer where this sample's pool was run, needed to locate file output by HiSeq
	v3: index			| Six digit sequence of the index for this library
	v4: ercc_mix		| Mix of ERCC spike used for library construction (options: "1", "2", "-")
	v5: file_directory	| Directory where sample's files were written to by Casava
	v6: batch			| GIGPAD batch number associated with sample
	v7: label			| Biological condition associated with the sample, provided by customer
	v8: ref_genome		| Rerence genome associated with sample. (options: "hg19", "Zv9", "mm10")
	v9: library_type	| Type of library for sample (options: "PE", "SE", "DGE", "SPE",
							corresponding to: "paired-end", "single-end", "digital gene expression", "stranded paired-end")
	"""
	f = open(fin,'r')
	c = f.read().split('\n')[1:]
	if '' in c:
		c.remove('')
	d = []
	for x in c:
		customer_id = x.split('\t')[0]
		gigpad_id = x.split('\t')[1]
		lane = x.split('\t')[2]
		index = x.split('\t')[3]
		ercc_mix = x.split('\t')[4]
		top_dir = x.split('\t')[5]
		if top_dir[-1] != "/":
			top_dir = top_dir+"/"
		batch = x.split('\t')[6]
		label = x.split('\t')[7]
		ref_genome = x.split('\t')[8]
		library_type = x.split('\t')[9]
		d.append([customer_id, gigpad_id, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type])
	return d


def get_genome_ref_files(genome):
	"""
	Location of all reference files needed for a given genome.
	For human, hg19 is preferred because GENCODE does not have gtf formatting needed for full use of Tuxedo Tools
	The ERCC gtf files were appended separately to each species own gtf file
	Current choices are: "hg19", "Zv9", "mm10"
	"""
	if genome == "hg19":
		ref_index = "/data/pcpgm/rnaseq/Indexes/hg19/hg19_ERCC"
		fa = "/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
		gtf = "/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
		ref = "/data/pcpgm/rnaseq/Indexes/hg19/refFlat.txt"
		ERCC_gtf = "/data/pcpgm/rnaseq/Indexes/hg19/hg19_ERCC_tuxedo.gtf"
	elif genome == "Zv9":
		ref_index = "/data/pcpgm/rnaseq/Indexes/Zv9/Zv9_ERCC"
		fa = "/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa"
		gtf = "/pub/genome_references/Zv9/Danio_rerio.Zv9.69.gtf"
		ref = "/data/pcpgm/rnaseq/Indexes/Zv9/refFlat_nochr.txt"
		ERCC_gtf = "/data/pcpgm/rnaseq/Indexes/Zv9/Zv9_ERCC.gtf"
	elif genome == "mm10":
		ref_index = "/data/pcpgm/rnaseq/Indexes/mm10/mm10_ERCC"
		fa = "/pub/genome_references/ensemble/Mmusculus/Mus_musculus.GRCm38.70.dna.toplevel.fa"
		gtf = "/pub/genome_references/ensemble/Mmusculus/Mus_musculus.GRCm38.70.gtf"
		ref = "/data/pcpgm/rnaseq/Indexes/mm10/refFlat_nochr.txt"
		ERCC_gtf = "/data/pcpgm/rnaseq/Indexes/mm10/mm10_ERCC.gtf"
	else:
		print 'Unknown genome selected: ', genome
	return(ref_index, fa, gtf, ref, ERCC_gtf)

	
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
	outp.write(">"+index_dictionary[index][0]+"\n"+index_dictionary[index][1]+"\n")
	outp.close()


def main(sample_info_file, discovery, standard_trim, path_start):
	"""
	Dispatches an lsf job to locate fastq files that were output by Casava (GIGPAD A1 routine) and then:
	1) Perform adapter trimming
	2) Run fastqc
	3) Get unique reads 
	4) Run tophat to align reads to reference genome
	5) Obtain various QC metrics on aligned files
	6) Run cufflinks to quantify ERCC spike ins (if applicable)
	7) Run cufflinks to quantify all other transcripts in sample
	
	The directory structure below is specific to the Partners HPC cluster
	"""
	runs = get_sample_info(sample_info_file)
	#for curr_sample, k in runs.iteritems():
	for k in runs:
		#Get sample information
		curr_sample, gigpad, lane, index, ercc_mix, top_dir, batch, label, ref_genome, library_type = k
		
		#Get genome reference files
		ref_index, fa, gtf, ref, ERCC_gtf = get_genome_ref_files(ref_genome)
		
		batch_dir = path_start+batch+"/"
		if not os.path.exists(batch_dir):
			os.makedirs(batch_dir)

		out_dir = batch_dir+curr_sample+"/"
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		
		job_name = curr_sample
				
		#This directory structure and naming convention is (obviously!) unique to PCPGM
		R1 = top_dir+"Unaligned-"+lane+"/Project_pcpgm/Sample_"+gigpad+"/filtered/"+gigpad+"_"+index+"_L00"+lane+"_R1.fastq"
		R2 = top_dir+"Unaligned-"+lane+"/Project_pcpgm/Sample_"+gigpad+"/filtered/"+gigpad+"_"+index+"_L00"+lane+"_R2.fastq"
		local_R1 = out_dir+curr_sample+"_R1.fastq"
		local_R2 = out_dir+curr_sample+"_R2.fastq"
		
		#Make lsf file		
		outp = open(job_name+".lsf", "w")
		outp.write("#!/bin/bash \n")
		outp.write("#BSUB -L /bin/bash\n")
		outp.write("#BSUB -J "+job_name+"\n")
		outp.write("#BSUB -q big-multi \n")
		outp.write("#BSUB -o "+job_name+"_%J.out\n")
		outp.write("#BSUB -e "+job_name+"_%J.screen\n")
		outp.write("#BSUB -R 'rusage[mem=24000]'\n")
		
		#Check whether unaligned fastq files that were processed by Casava are zipped and make local unzipped copies
		if not os.path.isfile(local_R1):
			if os.path.isfile(R1):
				outp.write("cp "+R1+" "+local_R1+"\n")
			elif os.path.isfile(R1+".gz"):
				outp.write("zcat "+R1+".gz > "+local_R1+"\n")
			else:
				print "R1 file not found ", R1
		if not os.path.isfile(local_R2):		
			if os.path.isfile(R2):
				outp.write("cp "+R2+" "+local_R2+"\n")
			elif os.path.isfile(R2+".gz"):
				outp.write("zcat "+R2+".gz > "+local_R2+"\n")
			else:
				print "R2 file not found", R2
		
		outp.write("cd "+out_dir+"\n")
		
		#Perform adapter trimming with trimmomatic
		#May perform a standard trimming of bases from reads by amount given above if standard_trim variable is greater than 0. Most will use standard_trim=0 
		#Create fa file of adapters specific to file
		if library_type in ["PE", "SE", "SPE"]:
			make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
		elif library_type == "DGE":
			make_adapter_fa(index, nextflex_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
		if standard_trim == 0:
			R1_trim = out_dir+curr_sample+"_R1_Trimmed.fastq"
			R2_trim = out_dir+curr_sample+"_R2_Trimmed.fastq"		
			if library_type in ["PE", "SPE"]:
				make_adapter_fa(index, illumina_indexes, out_dir+curr_sample+"_adapter.fa", library_type)
				outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+local_R1+" "+local_R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")
			elif library_type in ["SE", "DGE"]:				
				outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+local_R1+" "+R1_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")
		else:
			R1_trim = out_dir+curr_sample+"_R1_Trim"+str(standard_trim)+".fastq"
			R2_trim = out_dir+curr_sample+"_R2_Trim"+str(standard_trim)+".fastq"
			if library_type in ["PE", "SPE"]:
				outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticPE -phred33 "+local_R1+" "+local_R2+" "+R1_trim+" R1_Trimmed_Unpaired.fastq "+R2_trim+" R2_Trimmed_Unpaired.fastq HEADCROP:"+standard_trim+" ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")
			elif library_type in ["SE", "DGE"]:
				outp.write("java -Xmx1024m  org.usadellab.trimmomatic.TrimmomaticSE -phred33 "+local_R1+" "+R1_trim+" HEADCROP:"+standard_trim+"ILLUMINACLIP:"+out_dir+curr_sample+"_adapter.fa:2:30:10 MINLEN:50\n")					
								
		#Run fastqc on trimmed files.  
		#In some cases fastqc should be run on original files, but we have dropped this as a routine practice because the reports haven't changed much after trimming - adapter contamination has been minimal.
		if library_type in ["PE", "SPE"]:
			outp.write("fastqc -o "+out_dir+" "+R1_trim+" "+R2_trim+"\n")
		elif library_type in ("SE", "DGE"):
			outp.write("fastqc -o "+out_dir+" "+R1_trim+"\n")
		
		#Get total number of reads, unique reads, % unique reads from trimmed file(s). 
		outp.write("cat "+R1_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' > "+curr_sample+"_ReadCount\n")
		if library_type in ["PE", "SPE"]:
			outp.write("cat "+R2_trim+" | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}};print total,unique,unique*100/total}' >> "+curr_sample+"_ReadCount\n")
		
		#Run TopHat with options specific to library type
		if discovery == "no":
			if library_type == "PE":
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type == "SPE":
				outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type in ["DGE", "SE"]:
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" --no-novel-juncs --transcriptome-only -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
		
		elif discovery == "yes":
			if library_type == "PE":
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type == "SPE":
				outp.write("tophat --library-type fr-firststrand -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+" "+R2_trim+"\n")
			elif library_type in ["DGE", "SE"]:
				outp.write("tophat --library-type fr-unstranded -G "+ERCC_gtf+" -r 50 -p 12 "+ref_index+" "+R1_trim+"\n")
				
		#Get samtools mapping stats
		outp.write("cd "+out_dir+"/tophat_out/\n")
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
			outp.write("java -Xmx2g -jar /source/picardtools/picard-tools-1.58/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		else:
			outp.write("java -Xmx2g -jar /source/picardtools/picard-tools-1.58/CollectRnaSeqMetrics.jar REF_FLAT="+ref+" STRAND_SPECIFICITY=NONE INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_RNASeqMetrics\n")	
		
		#Get number of reads spanning junctions by getting "N"s in CIGAR field of bam file
		#Be sure that Junction Spanning Reads are Added First then Unmapped Reads for proper ordering of fields in report
		outp.write("echo \"Junction Spanning Reads: \" $(bamtools filter -in "+curr_sample+"_accepted_hits.sorted.bam -script /data/pcpgm/rnaseq/cigarN.script | bamtools count ) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")
		#Get number of unmapped reads
		outp.write("echo Unmapped Reads: $(samtools view -c unmapped.bam) >> "+curr_sample+"_accepted_hits.sorted.bamstats \n")		

		#Gather metrics unique to paired-end samples using CollectInsertSizeMetrics
		if library_type in ["PE", "SPE"]:
			outp.write("java -Xmx2g -jar /source/picardtools/picard-tools-1.58/CollectInsertSizeMetrics.jar HISTOGRAM_FILE="+curr_sample+"_InsertSizeHist.pdf INPUT="+curr_sample+"_accepted_hits.sorted.bam OUTPUT="+curr_sample+"_InsertSizeMetrics\n")	
		
		#Run cufflinks to count ERCC spike ins
		outp.write("mkdir "+out_dir+"/cufflinks_out_ERCC/\n")
		outp.write("cd "+out_dir+"/cufflinks_out_ERCC/\n")
		if library_type == "DGE":
			outp.write("cufflinks --library-type fr-unstranded --no-length-correction -G /data/pcpgm/rnaseq/ERCC/Ambion_Documents/ERCC92.gtf -p 12 ../tophat_out/accepted_hits.bam \n")
		elif library_type == "SPE":
			outp.write("cufflinks --library-type fr-firststrand -G /data/pcpgm/rnaseq/ERCC/Ambion_Documents/ERCC92.gtf -p 12 ../tophat_out/accepted_hits.bam \n")
		else:
			outp.write("cufflinks --library-type fr-unstranded -G /data/pcpgm/rnaseq/ERCC/Ambion_Documents/ERCC92.gtf -p 12 ../tophat_out/accepted_hits.bam \n")
		
		#Cufflinks to assemble and quantify transcripts
		outp.write("mkdir "+out_dir+"/cufflinks_out/\n")
		outp.write("cd "+out_dir+"/cufflinks_out/\n")
		if library_type == "DGE":
			outp.write("cufflinks --library-type fr-unstranded --no-length-correction -M /data/pcpgm/rnaseq/ERCC/Ambion_Documents/ERCC92.gtf -p 12 ../tophat_out/accepted_hits.bam \n")
		elif library_type == "SPE":
			outp.write("cufflinks --library-type fr-firststrand -M /data/pcpgm/rnaseq/ERCC/Ambion_Documents/ERCC92.gtf -p 12 ../tophat_out/accepted_hits.bam \n")			
		else:
			outp.write("cufflinks --library-type fr-unstranded -M /data/pcpgm/rnaseq/ERCC/Ambion_Documents/ERCC92.gtf -p 12 ../tophat_out/accepted_hits.bam \n")
		outp.close()
	
		subprocess.call("bsub < "+job_name+".lsf", shell=True)
		subprocess.call("mv "+job_name+".lsf "+out_dir, shell=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Write and execute an lsf job to perform QC and read alignment for RNA-seq samples associated with a PCPGM project.")
	parser.add_argument("--standard_trim", default=0, type=int, help="Number of bases to be trimmed from leftmost end of all reads (default=0)")
	parser.add_argument("--discovery", default="yes", type=str, help="Should TopHat be run with options to discover novel transcripts (i.e. disable --no-novel-juncs --transcriptome-only)? "
		"(options: yes, no; default=yes) "
		"Note: the 'no' option only works with hg19 at the moment")
	parser.add_argument("--path_start", default="./", type=str, help="Directory path where PCPGM batch-level directories are located and report directory will be written (default=./)")
	parser.add_argument("samples_in", help="Path to a tab-delimited txt file containing sample information. See example file: sample_info_file.txt")
	args = parser.parse_args()
	main(args.samples_in, args.discovery, args.standard_trim, args.path_start)

