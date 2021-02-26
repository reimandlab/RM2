#' Generates bins and calls prepare_sites such that site coordinates, number of trinucleotide contexts and indels, and all possible quadnucleotide contexts are returned for each mutation rate bin and site + flank
#'
#' @param sites Data frame of sites
#' @param window_size Integer indicating the half-width of sites and flanking regions
#' @param maf Data frame of mutations prepared by get_mut_trinuc_strand
#' @param n_bin Integer (default 10) indicating the number of bins to which to split sites
#' @param n_min_mut Integer (default 100) indicating minimum number of mutations to proceed with analysis
#' @param global_mut_rate_window Integer (default 1e6) indicating the window for mutation rate binning
#'
#' @return List of length n_bin each containing two GRanges of coordinates (site and flanks) and two dataframes of trinucleotide counts and mutation contexts
.prepare_bins_of_sites = function(sites, window_size, maf, n_bin=NA, n_min_mut = 100, global_mut_rate_window = 1e6) {

	# first count mutations in sites + flanks in total
	# if too low, then we will skip this analysis
	# do this early to avoid computing small problems
	window_size_flanked = 3 * window_size
	sites_mid = floor((sites$start + sites$end) / 2)
	gr_sites = GenomicRanges::GRanges(sites$chr, 
	                                  IRanges::IRanges(sites_mid - window_size_flanked, 
	                                  sites_mid + window_size_flanked - 1))
	
	# unique all sites to avoid counting muts multiple times
	gr_sites = GenomicRanges::reduce(gr_sites)
	gr_maf = GenomicRanges::GRanges(maf$chr, IRanges::IRanges(maf$start, maf$end))
	total_muts = sum(GenomicRanges::countOverlaps(gr_sites, gr_maf))
	
	if (total_muts < n_min_mut) {
	  stop("Too few mutations in sites+flanks.", call. = FALSE)
	}

	
	# prepare all sites together
	if (is.na(n_bin)) {
		sites_prepared = list("mutrate__1"= .prepare_sites2(sites, window_size))
	} else {

		# prepare sites for every bin based on average mutation rate
		dfr_1mb = data.frame(chr = sites$chr,
				start = sites_mid - global_mut_rate_window / 2,
				end = sites_mid + global_mut_rate_window / 2,
				stringsAsFactors = FALSE)

		# make sure windows don't exceed chr boundaries
		chr_ends = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[dfr_1mb$chr]
		chr_ends_index = which(dfr_1mb$end >= chr_ends)
		dfr_1mb$end[chr_ends_index] = chr_ends[chr_ends_index] - 1
		dfr_1mb$start[dfr_1mb$start < 1] = 1

		# count mutations per bin
		gr_1mb = GenomicRanges::GRanges(dfr_1mb$chr, IRanges::IRanges(dfr_1mb$start, dfr_1mb$end))
		site_mb_mut_counts = GenomicRanges::countOverlaps(gr_1mb, gr_maf)
		
		# some bins on chromosome boundaries and thus shorter; correct them
		bin_length_coefficient = 
		  ( dfr_1mb$end - dfr_1mb$start ) / global_mut_rate_window
		site_mb_mut_counts_corrected = 
		  site_mb_mut_counts / bin_length_coefficient

		# split sites into equal bins based on 1MB mut rates
		# return NULL; to be handled in higher-level functions
		bin_index = try(ggplot2::cut_number(site_mb_mut_counts, n_bin), silent = T)
		if (any(class(bin_index) == "try-error")) {
		  stop("Too few mutations to produce bins.", call. = FALSE)
		}

		site_bins = split(sites, bin_index)
		bin_mean_mut_rate = c(by(site_mb_mut_counts, bin_index, mean, na.rm=T))

		# prepare sites for every bin
		sites_prepared = lapply(site_bins, .prepare_sites2, window_size)
		names(sites_prepared) = paste("mutrate__", bin_mean_mut_rate, sep="")
	}
	return(sites_prepared)
}




#' Get merged and trimmed sites and flanking regions and trinucleotide counts
#'
#' @param this_sites Data frame of sites
#' @param window_size Integer indicating the half-width of sites and flanking regions
#'
#' @return List containing two GRanges of sites and flanks, and two dataframes for their trinucleotide counts and mutation contexts
.prepare_sites2 = function(this_sites, window_size) {
  
  window_size_flanked = 3 * window_size
  
  sites_mid = floor((this_sites$start + this_sites$end)/2)
  gr_sites = GenomicRanges::GRanges(this_sites$chr, 
                             IRanges::IRanges(sites_mid - window_size, 
                             sites_mid + window_size - 1))
  
  gr_sites_with_flank = GenomicRanges::GRanges(this_sites$chr, 
                                        IRanges::IRanges(sites_mid - window_size_flanked, 
                                        sites_mid + window_size_flanked - 1))
  
  # merge consecutive regions
  gr_sites = GenomicRanges::reduce(gr_sites)
  gr_sites_with_flank = GenomicRanges::reduce(gr_sites_with_flank)
  
  # remove sites from sites+flank to get only flank
  gr_flank = GenomicRanges::setdiff(gr_sites_with_flank, gr_sites)
  
  # trim sites and flanks that exceed chromosomal boundaries
  gr_sites = .trim_chrom_boundaries(gr_sites)
  gr_flank = .trim_chrom_boundaries(gr_flank)
  
  # count trinucleotides in both sites and flanks
  trinuc_sites = .get_trinucs2(gr_sites)
  trinuc_flank = .get_trinucs2(gr_flank)
  
  result = list(gr_sites, gr_flank, trinuc_sites, trinuc_flank)
  names(result) = c("gr_sites", "gr_flank", "trinuc_sites", "trinuc_flank")
  result
}



#' Clean sites near chromosome boundaries
#'
#' @param gr_sites_this GRanges object of sites
#'
#' @return GRanges object containing cleaned and trimmed sites
.trim_chrom_boundaries = function(gr_sites_this) {
  
  chr_ends = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[as.character(GenomeInfoDb::seqnames(gr_sites_this))]
  
  # these regions only both coordinates beyond boundaries, remove
  which_both_end_exceeds = 
    which(GenomicRanges::end(gr_sites_this) >= chr_ends & GenomicRanges::start(gr_sites_this) >= chr_ends)
  
  # these regions only both coordinates beyond boundaries, remove
  which_both_start_exceeds = 
    which(GenomicRanges::start(gr_sites_this) < 2 & GenomicRanges::end(gr_sites_this) < 2)
  sites_keep = setdiff(1:length(gr_sites_this), 
                       c(which_both_end_exceeds, which_both_start_exceeds))
  gr_sites_this = gr_sites_this[sites_keep]
  chr_ends = chr_ends[sites_keep]
  
  # these regions only have one coordinate beyond boundaries, trim
  which_end_exceeds = which(GenomicRanges::end(gr_sites_this) >= chr_ends)
  which_start_exceeds = which(GenomicRanges::start(gr_sites_this) < 2)
  
  # constrain by more than one because trinucleotide computing will grab one extra nt
  GenomicRanges::end(gr_sites_this)[which_end_exceeds] = chr_ends[which_end_exceeds] - 1
  GenomicRanges::start(gr_sites_this)[which_start_exceeds] = 2
  gr_sites_this
}
