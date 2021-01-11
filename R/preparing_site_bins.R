#' Generates bins and calls prepare_sites such that site coordinates, number of trinucleotide contexts and indels, and all possible quadnucleotide contexts are returned for each mutation rate bin and site + flank
#'
#' @param sites Data frame of sites
#' @param window_size Integer indicating the half-width of sites and flanking regions
#' @param maf Data frame of mutations prepared by get_mut_trinuc_strand
#' @param n_bin Integer (10) indicating the number of bins to which to split sites
#' @param n_min_mut Integer (100) indicating minimum number of mutations to proceed with analysis
#' @param global_mut_rate_window Integer (1e6) indicating the window for mutation rate binning
#'
#' @return List of length n_bin each containing 3 GRanges of coordinates (site and flanks) and 3 dataframes of trinucleotide counts and mutation contexts

.prepare_bins_of_sites = function(sites, window_size, maf, n_bin=NA, n_min_mut = 100, global_mut_rate_window = 1e6) {

	# first count mutations in sites + flanks in total. Skip analysis if too low
	gr_maf = GenomicRanges::GRanges(maf$chr, IRanges::IRanges(maf$start, maf$end))
	sites_mid = floor((sites$start+sites$end)/2)
	gr_sites = GenomicRanges::GRanges(sites$chr,
		IRanges::IRanges(sites_mid - 3 * window_size - 1, sites_mid + 3 * window_size))
	total_muts = sum(GenomicRanges::countOverlaps(gr_sites, gr_maf))

	if (total_muts < n_min_mut) {
		print(paste0("Too few mutations in sites+flanks (", total_muts, "), skipping."))
		return(NULL)
	}

	# remove sites that together with windows exceed chromosomal coordinates
	chr_ends = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[as.character(GenomeInfoDb::seqnames(gr_sites))]
	chr_ends_index = which(end(gr_sites) >= chr_ends)
	chr_starts_index = which(start(gr_sites) <= 1)
	keep_site_index = setdiff(1:nrow(sites), c(chr_ends_index, chr_starts_index))
	sites = sites[keep_site_index,, drop = FALSE ]
	sites_mid = sites_mid[keep_site_index]

	# prepare all sites together
	if (is.na(n_bin)) {
		sites_prepared = list("mutrate__1"=.prepare_sites(sites, window_size))
	} else {

		# prepare sites for every bin based on average mutation rate
		dfr_1mb = data.frame(chr = sites$chr,
				start = sites_mid - global_mut_rate_window / 2,
				end = sites_mid + global_mut_rate_window / 2 - 1,
				stringsAsFactors = FALSE)

		# make sure windows don't exceed chr boundaries
		chr_ends = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[dfr_1mb$chr]
		chr_ends_index = which(dfr_1mb$end >= chr_ends)
		dfr_1mb$end[chr_ends_index] = chr_ends[chr_ends_index] - 1
		dfr_1mb$start[dfr_1mb$start < 1] = 1

		# count mutations per bin
		gr_1mb = GenomicRanges::GRanges(dfr_1mb$chr, IRanges::IRanges(dfr_1mb$start, dfr_1mb$end))
		site_mb_mut_counts = GenomicRanges::countOverlaps(gr_1mb, gr_maf)

		# split sites into equal bins based on 1MB mut rates
		# return NULL; to be handled in higher-level functions
		bin_index = try(ggplot2::cut_number(site_mb_mut_counts, n_bin), silent = T)
		if (any(class(bin_index) == "try-error")) {
					stop(paste0("Too few mutations to produce ", n_bin, " bins."))
		}

		site_bins = split(sites, bin_index)
		bin_mean_mut_rate = c(by(site_mb_mut_counts, bin_index, mean, na.rm=T))

		# prepare sites for every bin
		sites_prepared = lapply(site_bins, .prepare_sites, window_size)
		names(sites_prepared) = paste("mutrate__", bin_mean_mut_rate, sep="")
	}
	return(sites_prepared)
}




#' For each mutation rate bin and site + flanks, calculates site coordinates, number of trinucleotide contexts and indels, and all possible mutation contexts
#'
#' @param sites Data frame of sites
#' @param window_size Integer indicating the half-width of sites and flanking regions
#'
#' @return List containing 3 GRanges of coordinates (site and flanks) and 3 dataframes of trinucleotide counts and mutation contexts
.prepare_sites = function(sites, window_size) {

	sites_mid = floor((sites$start+sites$end)/2)
	gr_sites = GenomicRanges::GRanges(sites$chr, IRanges::IRanges(sites_mid - window_size - 1, sites_mid + window_size))

	site_ids = paste(GenomicRanges::seqnames(gr_sites), GenomicRanges::start(gr_sites), GenomicRanges::end(gr_sites), sep="_")

	# create unique site IDs
	dupl_site_index = which(duplicated(site_ids))
	site_ids[dupl_site_index] = paste(site_ids[dupl_site_index], 1:length(dupl_site_index), sep="_")

	# compare each site with immediate background of same length
	gr_sites_upstream = .get_up_flank(gr_sites)
	gr_sites_downstream = .get_down_flank(gr_sites)

	trinuc_sites = .get_nucs(gr_sites, site_ids)
	trinuc_upstream = .get_nucs(gr_sites_upstream, site_ids)
	trinuc_downstream = .get_nucs(gr_sites_downstream, site_ids)

	list(gr_sites, gr_sites_upstream, gr_sites_downstream, trinuc_sites, trinuc_upstream, trinuc_downstream)
}
