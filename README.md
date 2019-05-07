# combine_uc.r

  	# Written by John L. Darcy, May 2019
	# See my website at jldarcy.tk
	
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
		vsearch --derep_fulllength rawreads_renamed.fasta --output derep_seqs.fasta --uc seq2derep_uctable.txt --sizeout
	# 2. Run unoise on dereplicated reads to make zOTU seeds
		usearch -unoise3 derep_seqs.fasta -minsize 4 -zotus esv_seqs.fasta
	# 3. Map dereplicated reads onto zOTU seeds
		vsearch --usearch_global derep_seqs.fasta --db esv_seqs.fasta --id 0.97 --threads 12 --maxaccepts 0 --maxrejects 0 --uc derep2esv_uctable.txt
	# 4. Run this script
		./combine_uc.r --lo seq2derep_uctable.txt.gz --hi derep2esv_uctable.txt.gz -t 12 -o esv_table.txt
