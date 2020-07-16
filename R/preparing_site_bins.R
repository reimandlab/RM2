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

	# prepare all sites together
	if (is.na(n_bin)) {
		sites_prepared = list("mutrate__1"=.prepare_sites(sites, window_size))
	} else {

		# prepare sites for every bin based on average mutation rate and set fixed background if boundaries exceeded
		legal_chr = paste0("chr", c(1:22, "M", "X", "Y"))
		chr_max = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[legal_chr]

		sites_1mb = do.call(rbind, parallel::mclapply(1:length(sites_mid), .add_fixed_background, sites_mid, sites$chr, chr_max, global_mut_rate_window/2, mc.cores=4))
		gr_1mb = GenomicRanges::GRanges(sites_1mb$site_chr, IRanges::IRanges(sites_1mb$start, sites_1mb$end))

		site_mb_mut_counts = GenomicRanges::countOverlaps(gr_1mb, gr_maf)

		# split sites into equal bins based on 1MB mut rates
		# return NULL; to be handled in higher-level functions
		bin_index = try(ggplot2::cut_number(site_mb_mut_counts, n_bin), silent = T)
		if (any(class(bin_index)=="try-error")) {
			print(paste0("Insufficient data values to produce ", n_bin, " bins."))
			return(NULL)
		}

		site_bins = split(sites, bin_index)
		bin_mean_mut_rate = c(by(site_mb_mut_counts, bin_index, mean, na.rm=T))

		# prepare sites for every bin
		sites_prepared = lapply(site_bins, .prepare_sites, window_size)
		names(sites_prepared) = paste("mutrate__", bin_mean_mut_rate, sep="")
	}
	return(sites_prepared)
}


#' Ensures background megabase mutation rates use legal coordinates
#'
#' @param site_index Numeric index of site
#' @param sites_mid Numeric vector of site midpoints
#' @param sites_chr Character vector of site chromosomes
#' @param chr_max Numeric vector of max coordinate per chromosome
#' @param background_width Integer indicating background size (1e6)
#'
#' @return Data frame of legal site chromosomes and coordinates
.add_fixed_background = function(site_index, sites_mid, sites_chr, chr_max, background_width) {
	site_mid = sites_mid[site_index]
	site_chr = as.character(sites_chr[site_index])
	max_coord = chr_max[[site_chr]]

	if (site_mid > background_width & site_mid < max_coord - background_width) {
		start = site_mid - background_width
		end = site_mid + background_width
	} else if (site_mid < background_width) {
		start = 1
		end = 2 * background_width + 1
	} else if (site_mid > max_coord - background_width) {
		end = max_coord
		start = max_coord - 2 * background_width
	}
	return(data.frame(site_chr, start, end, stringsAsFactors = F))
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
