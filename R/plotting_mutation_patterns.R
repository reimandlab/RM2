#' Counts mutations at positions along flanks and sites
#'
#' @param mutations Data frame of mutations
#' @param sites Data frame of sites
#' @param window_size Half the of site width
#'
#' @return Data frame containing the relative position of each base (numeric), count (numeric) of mutations and site status (boolean)
#' @export
get_mutations_in_flanked_sites = function(mutations, sites, window_size) {

  sites_mid = floor((sites$start+sites$end)/2)
  gr_sites = GenomicRanges::GRanges(sites$chr, IRanges::IRanges(sites_mid - window_size, sites_mid + window_size - 1))

  gr_maf = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations$start, mutations$end))

  # unique mutations identifier
	GenomicRanges::mcols(gr_maf)[,1] = paste(mutations[,1], mutations[,2], mutations[,3], mutations[,4],
							mutations[,5], mutations[,7], sep = ":")

  gr_sites_flank = gr_sites
  GenomicRanges::start(gr_sites_flank) = GenomicRanges::start(gr_sites_flank) - GenomicRanges::width(gr_sites) - 1
  GenomicRanges::end(gr_sites_flank) = GenomicRanges::end(gr_sites_flank) + GenomicRanges::width(gr_sites) + 1

  gr_this_maf = gr_maf[S4Vectors::subjectHits(GenomicRanges::findOverlaps(gr_sites_flank, gr_maf))]
	gr_this_sites_flank = gr_sites_flank[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_sites_flank, gr_maf))]

	# remove duplicated mutations
	index_keep = which(!duplicated(GenomicRanges::mcols(gr_this_maf)[,1]))
	gr_this_maf = gr_this_maf[index_keep]
	gr_this_sites_flank = gr_this_sites_flank[index_keep]

  # get relative position
  mut_rel_pos = (GenomicRanges::start(gr_this_maf) -
                                  floor(0.5 * (GenomicRanges::start(gr_this_sites_flank) + GenomicRanges::end(gr_this_sites_flank)))) + 1

  dfr = data.frame(dist=as.numeric(names(table(mut_rel_pos))),
                                  count=c(table(mut_rel_pos)),
                                  is_site=FALSE, stringsAsFactors=F)
  dfr = dfr[dfr$dist >= -3*window_size & dfr$dist <= 3*window_size,]

  dfr[dfr$dist >= -window_size & dfr$dist <= window_size, "is_site"] = TRUE
  dfr
}

#' Visualizes mutations in sites and flanks
#'
#' @param dfr Data frame of mutation counts at relative positions along stacked sites/flanks
#'\describe{
#'     \item{dist}{relative position of mutations}
#'     \item{count}{number of mutations at relative position}
#'     \item{is_site}{boolean indicator of whether relative position is within the site}
#' }
#' @param window_size Half of the site width. Added to each side of site midpoint
#'
#' @return ggplot object
#' @export
plot_mutations_in_flanked_sites = function(dfr, window_size=100, y_range=NULL) {
  p = ggplot2::ggplot(dfr, ggplot2::aes(dist, count, fill=is_site)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::geom_smooth(ggplot2::aes(dist, count), fill="blue", method="loess", span=0.333) +
    ggplot2::scale_x_continuous("Relative position",
                       breaks=c(0, window_size, -window_size,
                                  2*window_size, -2*window_size,
                                  3*window_size, -3*window_size)) +
    ggplot2::scale_y_continuous("# mutations", limits=y_range) +
    ggplot2::scale_fill_manual("", labels=c("TRUE"="site", "FALSE"="flank"), values=c("TRUE"="darkred", "FALSE"="gray")) +
    ggplot2::theme_bw()
  p
}
