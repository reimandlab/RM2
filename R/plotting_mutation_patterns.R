#' Mutation counts along relative positions
#'
#' Prepares data frame of mutations at relative positions along flanks and sites
#' @param mutations Data frame of mutations
#' @param this_sites Data frame of sites
#' @param window_size Integer indicating the half-width of sites and flanking regions (added to both sides of midpoint)
#' @param n_patients Integer indicating number of patients
#' 
#' @return Data frame containing the relative position of each base (numeric), mut_count (numeric) of mutations, 
#' mut_freq (numeric) calulated as mut_count / n_patients / total bp * 1e6, and site status (boolean)
#' @export
get_mutations_in_flanked_sites = function(mutations, this_sites, window_size, n_patients) {
  
  sites_mid = floor((this_sites$start + this_sites$end) / 2)
  gr_sites_mid = GenomicRanges::GRanges(this_sites$chr, IRanges::IRanges(sites_mid, sites_mid))
  
  mutations_mid = floor((mutations$start + mutations$end) / 2)
  gr_maf_mid = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations_mid, mutations_mid))
  
  ov_nearest = GenomicRanges::distanceToNearest(gr_maf_mid, gr_sites_mid)	
  
  nearest_muts_mid = mutations_mid[S4Vectors::queryHits(ov_nearest)]
  nearest_muts_chr = mutations[S4Vectors::queryHits(ov_nearest), "chr"]
  
  nearest_sites_mid = sites_mid[S4Vectors::subjectHits(ov_nearest)]
  nearest_sites_chr = this_sites[S4Vectors::subjectHits(ov_nearest), "chr"]
  
  dfr_dist = data.frame (
    muts_chr = nearest_muts_chr, muts_mid = nearest_muts_mid,
    sites_chr = nearest_sites_chr, sites_mid = nearest_sites_mid, 
    stringsAsFactors = FALSE	
  )
  
  dfr_dist$min_dist = dfr_dist$muts_mid - dfr_dist$sites_mid
  
  # sites on the same chromosome as muts; min distance has to be within windows and flanks
  dfr_dist = dfr_dist[
    dfr_dist$muts_chr == dfr_dist$sites_chr & 
      abs(dfr_dist$min_dist) <= 3 * window_size ,, drop = FALSE]
  
  dfr_freq = reshape2::dcast(dfr_dist, min_dist ~ . , value.var = "min_dist", 
                   fun.aggregate = length)
  colnames(dfr_freq) = c("dist", "mut_count")
  
  dfr_freq$is_site = FALSE
  dfr_freq[abs(dfr_freq$dist) <= window_size, "is_site"] = TRUE
  
  # normalising constant to get mutation frequency besides count
  gr_sites_flanks = GenomicRanges::GRanges(this_sites$chr, 
                                           IRanges::IRanges (sites_mid - 3 * window_size, sites_mid + 3 * window_size - 1))
  gr_sites_flanks = GenomicRanges::reduce(gr_sites_flanks)
  total_bp = sum(IRanges::width(gr_sites_flanks))
  total_bp_per_site_position = total_bp / (6 * window_size)
  
  dfr_freq$mut_freq = dfr_freq$mut_count / n_patients / (total_bp_per_site_position / 1e6)
  
  dfr_freq
}



#' Visualizes mutations in sites and flanks
#'
#' Barplot of mutations along relative positions to site midpoints
#' @param dfr Data frame of mutation counts at relative positions along stacked sites/flanks
#'\describe{
#'     \item{dist}{relative position of mutations}
#'     \item{count}{number of mutations at relative position}
#'     \item{is_site}{boolean indicator of whether relative position is within the site}
#' }
#' @param window_size Integer indicating the half-width of sites and flanking regions (added to both sides of midpoint)
#' @param what_plot Character (default "mut_freq") indicating y-axis variable 
#' @param y_range Numeric y-axis limits
#'
#' @return ggplot object
#' @export
plot_mutations_in_flanked_sites = function(dfr, window_size=100, what_plot = "mut_freq", y_range=NULL) {
  
  dfr1 = dfr[, c(what_plot, "dist", "is_site")]
  colnames(dfr1) = c("value", "dist", "is_site")
  
  p = ggplot2::ggplot(dfr1, ggplot2::aes(dist, value, fill=is_site)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::geom_smooth(formula = y ~ x, ggplot2::aes(dist, value), fill="blue", method="loess", span=1/3) +
    ggplot2::scale_x_continuous("Distance to site midpoint (bps)",
                       breaks=c(0, window_size, -window_size,
                                  2*window_size, -2*window_size,
                                  3*window_size, -3*window_size)) +
    ggplot2::scale_y_continuous(what_plot) +
    ggplot2::coord_cartesian(ylim = y_range) +
    ggplot2::scale_fill_manual("", labels=c("TRUE"="site", "FALSE"="flank"), values=c("TRUE"="darkred", "FALSE"="gray")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 10),
          axis.text = ggplot2::element_text(colour = "black"))
  p
}
