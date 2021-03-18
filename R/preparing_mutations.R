#' Derive quadnucleotide contexts, strandedness and single-base substitutions
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
#' @return Data frame containing additional mutation context columns
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 1 coordinates}
#'     \item{end}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{mut_trinuc}{pyrimidine-centered quadnucleotide mutation context and indel}
#'     \item{mut_strand}{the mutation strand classified as w for Watson and c for Crick}
#'     \item{ref_alt}{pyrimidine single-base substitution}
#' }
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