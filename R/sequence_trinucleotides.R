#' Counts the number of each trinucleotide within each site
#'
#' @param gr_sites GenomicRanges object of sites
#' @param side_ids Character vector of site ids
#'
#' @return Data frame containing each site id (character) and its trinucleotide counts (numeric)
.get_sequence_trinucleotides = function(gr_sites, site_ids) {

	all_trinucs = .get_all_trinucs()

	# extend by one nucleotide to get the neighbouring context for each nucleotide of the sequence
	GenomicRanges::start(gr_sites) = GenomicRanges::start(gr_sites) - 1
	GenomicRanges::end(gr_sites) = GenomicRanges::end(gr_sites) + 1

	this_seq = strsplit(as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, gr_sites)), '')

	# extract trinucleotides from sequence
	get_sq_trinucs = function(x) sapply(1:(length(x)-2), function(i) paste0(x[i:(i+2)], collapse=""))
	tri_nucleotide = lapply(this_seq, get_sq_trinucs)

	# convert to tall format, map A/G-centred trinucleotides to complementary strand
	tri_nucleotide = do.call(rbind, lapply(1:length(gr_sites), function(i) cbind(ID=site_ids[i], trinuc=tri_nucleotide[[i]])))
	where_reverse = substr(tri_nucleotide[, "trinuc"], 2, 2) %in% c("A", "G")
	tri_nucleotide[where_reverse, "trinuc"] =
		as.character(Biostrings::complement(Biostrings::DNAStringSet(tri_nucleotide[where_reverse, "trinuc"])))

	# add up trinucleotides and convert into a wide format
	tri_nucleotide = data.frame(tri_nucleotide[,1:2], count=1, stringsAsFactors=F)
	trinuc_dfr = reshape2::dcast(ID~trinuc, data=tri_nucleotide, value.var="count", fun.aggregate=sum)

	# add missing quad-nucleotide classes as zero columns
	missing_trinucs = setdiff(all_trinucs, colnames(trinuc_dfr))
	dfr_tail = matrix(0, nrow(trinuc_dfr), length(missing_trinucs))
	colnames(dfr_tail) = missing_trinucs
	trinuc_dfr = cbind(trinuc_dfr, dfr_tail)

	# indel can occur at any given position, add up trinucleotide positions and make new column for indel
	trinuc_dfr$indel = apply(trinuc_dfr[,all_trinucs], 1, sum)
	trinuc_dfr[,c("ID", all_trinucs, "indel")]
}


#' Derive quadnucleotide contexts, strandedness and single-base substitution
#'
#' get_mut_trinuc_strand() takes a data frame of mutations and determines the quadnucleotide context (trinucleotide context + alternate allele) and strand.
#' @param maf Data frame of mutations with the following columns: chr, start, end, ref and alt
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 1 coordinates}
#'     \item{end}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#' }
#'
#' @return Dataframe containing quadnucleotide context (character), the strand, classified as w for Watson and c for Crick (character), and the single-base substitution (character)
#' @export
get_mut_trinuc_strand = function(maf) {

  gr_trinuc = GenomicRanges::GRanges(maf$chr, IRanges::IRanges(start=maf$start-1, end=maf$end+1), strand="*")
  triples = as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, gr_trinuc))

  legal_nucl = c("A", "C", "G", "T")
  complement_nucl = c("A", "G") ## these nucleotides will be replaced by complementary nucleotide

  mut_strand = mut_class = rep("indel", length(triples))

  mut_class[maf$ref %in% legal_nucl & maf$alt %in% legal_nucl] = "SNV"
  mut_strand[mut_class=="SNV"] = "w"

  mid_nucl = gsub("(.)(.)(.)", "\\2", triples)

  which_to_complement = which(mut_class=="SNV" & mid_nucl %in% complement_nucl)

  new_triples = triples
  new_alt = maf$alt

  new_triples[which_to_complement] = as.character(Biostrings::complement(Biostrings::DNAStringSet(new_triples[which_to_complement])))
  new_alt[which_to_complement] = as.character(Biostrings::complement(Biostrings::DNAStringSet(new_alt[which_to_complement])))

  mut_strand[which_to_complement] = "c"

  mut_trinuc = paste0(new_triples,"_",new_alt)
  mut_trinuc[mut_class=="indel"] = "indel"

  ref_alt = paste0(gsub("^.(.).$", "\\1", new_triples), "_", new_alt)
	ref_alt[mut_class=="indel"] = "indel"

	data.frame(mut_trinuc, mut_strand, ref_alt, stringsAsFactors=F)
}
