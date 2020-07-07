#' Determines upstream coordinates of flanking regions to sites
#'
#' @param gr_sites GenomicRanges object of sites
#'
#' @return GenomicRanges object containing the coordinates of upstream flanking regions
.get_up_flank = function(gr_sites) {
	gr_sites_upstream = gr_sites
	GenomicRanges::start(gr_sites_upstream) = GenomicRanges::start(gr_sites) - GenomicRanges::width(gr_sites)
	GenomicRanges::end(gr_sites_upstream) = GenomicRanges::start(gr_sites) - 1
	gr_sites_upstream
}


#' Determines downstream coordinates of flanking regions to sites
#'
#' @param gr_sites GenomicRanges object of sites
#'
#' @return GenomicRanges object containing the coordinates of downstream flanking regions
.get_down_flank = function(gr_sites) {
	gr_sites_downstream = gr_sites
	GenomicRanges::start(gr_sites_downstream) = GenomicRanges::end(gr_sites) + 1
	GenomicRanges::end(gr_sites_downstream) = GenomicRanges::end(gr_sites) + GenomicRanges::width(gr_sites)
	gr_sites_downstream
}


#' Determines trinucleotide, pyrimidine-centered trinucleotide contexts
#'
#' @return Character vector of all possible trinucleotide mutation contexts with pyrimidines in the middle position
.get_all_trinucs = function() {
	apply(as.matrix(expand.grid(c("A", "C", "G", "T"), c("C", "T"), c("A", "C", "G", "T"))), 1, paste, collapse="")
}


#' Creates a data frame of quadnucleotide mutation contexts with counts of corresponding trinucleotides
#'
#' @param gr_sites GenomicRanges object of sites
#' @param site_ids Character vector of unique site ids
#'
#' @return Data frame of all possible quadnucleotide mutation contexts with the number of positions at which the trinucleotides occur
.get_nucs = function(gr_sites, site_ids) {
	quadnuc = .get_quadnuc_map()
	trinuc_counts = .get_sequence_trinucleotides(gr_sites, site_ids)
	trinuc_counts1 = apply(trinuc_counts[,-1], 2, sum)
	trinuc_counts1 = data.frame(trinuc=names(trinuc_counts1), posi_count=trinuc_counts1, stringsAsFactors=F)
	trinuc_counts2 = merge(trinuc_counts1, quadnuc, by="trinuc")
	trinuc_counts2
}


#' Determines all trinucleotide, pyrimidine-centered trinucleotide contexts and mutation/quadnucleotide contexts
#'
#' @return Data frame trinucleotides, alternate allele and possible quadnucleotide contexts
.get_quadnuc_map = function() {
	quadnuc = S4Vectors::expand.grid(.get_all_trinucs(), c("A", "C", "G", "T"))
	quadnuc$tag = paste(quadnuc[,1], quadnuc[,2], sep="_")
	quadnuc = quadnuc[gsub("^.(.).$", "\\1", quadnuc[,1]) != quadnuc[,2],]
	quadnuc = data.frame(as.matrix(quadnuc), stringsAsFactors=F)
	quadnuc = rbind(quadnuc, c("indel", "indel", "indel"))
	colnames(quadnuc) = c("trinuc", "alt", "quadnuc")
	quadnuc
}


#' Selects median value with priority given to larger value
#'
#' @param x Numeric vector of p-values
#'
#' @return Position of median value
#' @export
which_median = function(x) {
	mid_pos = ceiling(length(x)/2)
	which(x == sort(x)[mid_pos])[1]
}
