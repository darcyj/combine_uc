# combine_uc.r

  	# Written by John L. Darcy, June 2017
	# See my website at http://amo.colorado.edu/jld/
	
	# Script to combine uc tables produced by dereplicating sequence reads and running 
	# unoise. This speeds up the process of creating a zOTU table using unoise. 

	# This program can be used to speed up making an OTU table using
	# zero-radius OTUs (zOTUs) picked by unoise. Making a zOTU table
	# can be done by mapping all raw sequence reads back onto zOTU seeds,
	# but this is very time consuming. Since dereplicated reads are 
	# clustered at 100% identity, it's fair to use centroid sequences
	# in lieu of all raw reads. This speeds up the process greatly. 
	
	# The R packages optparse and data.table are required to use this program.
	
	# Pipeline:
	
	# 0. Format raw read fasta headers to be in usearch format. 
		# this looks like >seqID;barcodelabel=sampleID
	# 1. Dereplicate raw reads
		# vsearch --derep_fulllength rawreads_renamed.fasta --output derep_seqs.fasta --uc derep_uctable.txt --sizeout
	# 2. Run unoise on dereplicated reads
		# usearch -unoise3 derep_seqs.fasta -minsize 4 -zotus zotu_seeds.fasta
	# 3. Map dereplicated reads onto zOTU seeds
		# vsearch --usearch_global derep_seqs.fasta --db zotu_seeds.fasta --id 0.97 --threads 24 --maxaccepts 0 --maxrejects 0 --uc unoise_uctable.txt
	# 4. Run this script
		# combine_uc.r -d derep_uctable.txt -u unoise_uctable.txt -t 20 -o zOTUtable.txt
