#!/usr/bin/env Rscript

## This program combines two UC tables, where one is nested within the other.
	# output is an otu table.

## load in packages
	suppressPackageStartupMessages(require(optparse))
	suppressPackageStartupMessages(require(parallel))
	suppressPackageStartupMessages(require(data.table))

## make options for script
	option_list <- list(
		make_option(c("--lo"), action="store", default=NA, type='character',
			help="Lower-level uc table, produced first. Larger file."), 
		make_option(c("--hi"), action="store", default=NA, type='character',
			help="Higher-level uc table, produced second. Smaller file."), 
		make_option(c("-o", "--output"), action="store", default=NA, type='character',
			help="Output filename for summary table (OTU table)"),
		make_option(c("-t", "--threads"), action="store", default=4, type='integer',
			help="Number of threads to use for parallel processing.")
	)
	
	# parse option arguments
	opt = parse_args(OptionParser(option_list=option_list))

## just here for testing
	if(FALSE){
		opt=list(lo="seq2derep_uctable.txt.gz", hi="derep2esv_uctable.txt.gz", threads=12, output="esv_table.txt")
	}


## read in input files
	# in each uc table, col 9 is the sequence that is grouped into the ESV/centroid in col 10.
	
	# uc_lo is the lower-level uctable, i.e. uc_hi will be nested within it.
	# it will have been generated first, and be a much larger file than uc_hi.
	message("Reading in lower-level uc file...")
	uc_lo <- fread(opt$lo, sep='\t', header=F, stringsAsFactors=F)
	message("   ...done.")
	
	# uc_hi is the higher-level uctable
	# it will have been generated second, and be a much smaller file than uc_hi.
	message("Reading in higher-level uc file...")
	uc_hi <- fread(opt$hi, sep='\t', header=F, stringsAsFactors=F)
	message("   ...done.")

## prune useless data from uc tables
	message("Pruning inputs...")

	# remove S and N rows. These are "no-hit" and "summary" rows.
	# thus, KEEP "C" and "H" rows. These are "centroid" and "hit" rows. 
	uc_lo <- uc_lo[uc_lo[[1]] %in% c("C", "H"),]
	uc_hi <- uc_hi[uc_hi[[1]] %in% c("C", "H"),]
	# make centroids self-referential so they can be counted.
	# otherwise, the count of each lower-level cluster will be off-by-one.
	uc_lo[[10]][uc_lo[[10]]=="*"] <- uc_lo[[9]][uc_lo[[10]]=="*"]

	# remove useless columns:
	uc_lo <- uc_lo[, 9:10]
	uc_hi <- uc_hi[, 9:10]
	message("   ...done.")\

## check inputs before computationally intensive steps
	message("Quickly checking inputs...")
	if(grepl(pattern=";barcodelabel=", x=uc_lo[[1]][1]) && grepl(pattern=";barcodelabel=", x=uc_lo[[2]][1])){
		message("...lower-level uc file looks OK.")
	}else{
		stop("ERROR: no barcodelabel found in lower-level uc file.")
	}
	if(grepl(pattern=";barcodelabel=", x=uc_hi[[1]][1])){
		message("...higher-level uc file looks OK.")
	}else{
		stop("ERROR: no barcodelabel found in higher-level uc file.")
	}




## digest labels into component parts (seqid, barcodelabel)
	# to enable quick boolean operations
	get_bclabel <- function(h){
		h_sep <- unlist(strsplit(h, split=";"))
		return(sub("barcodelabel=", "", h_sep[grep("barcodelabel=", h_sep)]))
	}
	get_seqid <- function(h){
		return(unlist(strsplit(h, split=";"))[1])
	}
	# uc_lo columns become sample and seq, uc_hi cols become seq and esv
	message("Extracting sample info...")
	colnames(uc_hi) <- c("seq", "esv")
	colnames(uc_lo) <- c("sample", "seq")
	uc_lo$sample <- as.factor(simplify2array(mclapply(X=uc_lo$sample, FUN=get_bclabel, mc.cores=opt$threads)))
	message("   ...1/3 columns done.")
	uc_lo$seq <- simplify2array(mclapply(X=uc_lo$seq, FUN=get_seqid, mc.cores=opt$threads))
	message("   ...2/3 columns done.")
	uc_hi$seq <- simplify2array(mclapply(X=uc_hi$seq, FUN=get_seqid, mc.cores=opt$threads))
	message("   ...3/3 columns done.")

## for each zOTU, calculate abundances across samples.
	message("Making final table...")
	unique_samples <- unique(uc_lo$sample)
	unique_esvs <- unique(uc_hi$esv)
	n_esvs <- length(unique_esvs)

	# break up esvs into chunks of bmult*ncores and do 'em
	# i wrote this without putting enough comments in my code and now i don't know how it works
	# don't do this, kiddos
	bmult <- 3
	esv_bins <- rep(1:(ceiling(n_esvs / opt$threads /bmult) ), each=opt$threads * bmult)[1:n_esvs]
	process_bin <- function(bin_number){
		get_row <- function(esv){ return(table(uc_lo$sample[uc_lo$seq %in% uc_hi$seq[uc_hi$esv == esv]])) }

			get_ESV_row <- function(ESV){
				# find which derep centroids are within ESV_i, get their sampleids too
				ESV_rowsTF <- unoise_uc[[10]] == ESV
				centroids_i <- unoise_uc[[9]] [ESV_rowsTF]
				
				# get all sampleids of original sequences within those centroids (incl. centroids themselves)
				sampids_i <- c(derep_uc_sampleids[derep_uc[[10]] %in% centroids_i], centroid_sampleids [ESV_rowsTF])
				
				# make sampids_i a factor so that table returns empty fields
				# this makes it work when an OTU isn't observed in a sample.
				# note that same derep_uc_sampleids_unique is used later as the colnames for the ESV table
				sampids_i <- factor(sampids_i, levels=derep_uc_sampleids_unique)
				
				# make otu table row. table() will be sorted in the same order as derep_uc_sampleids_unique.
				return(as.numeric(table(sampids_i)))
			}


		esvs_i <- unique_esvs[esv_bins == bin_number]
		submat <- simplify2array(mclapply(X=esvs_i, FUN=get_row, mc.cores=opt$threads)) 
		colnames(submat) <- esvs_i
		return(submat)
	}


## process each bin in parallel, but do bins serially. this is a compromise that
	# saves memory use and allows for a progress bar with no overhead.
	esvlist <- vector("list", max(esv_bins))
	message("...processing bins of ESVs in parallel...")
	pb <- txtProgressBar(min=0, max=max(esv_bins), style=3)
	for(b in 1:max(esv_bins)){
		esvlist[[b]] <- process_bin(b)
		setTxtProgressBar(pb, b)

	}
	message("")



	# esvlist <- lapply(X=1:max(esv_bins), FUN=process_bin)
	message("...writing out table...")
	otutable <- t(do.call("cbind", esvlist))

	# turn otutable into list of cols so it can be more quickly written out with fwrite
	# cols are used because fwrite transposes output for whatever reason.
	cols2write <- lapply(X=1:ncol(otutable), FUN=function(i){as.vector(c(colnames(otutable)[i], otutable[,i]))})
	# add OTU_ID row
	cols2write <- append( list(c("OTU_ID", rownames(otutable))), cols2write)



	fwrite(file=opt$output, x=cols2write, row.names=FALSE, sep='\t', quote=FALSE)

	message("   ...all done.")
