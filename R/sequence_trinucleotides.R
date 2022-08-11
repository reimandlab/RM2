#' Return possible pyrimidine-centered trinucleotide contexts
#'
#' @return Character vector of all possible trinucleotide mutation contexts with pyrimidines in the middle position
.get_all_trinucs = function() {
  nuc_trios = as.matrix(S4Vectors::expand.grid(
    c("A", "C", "G", "T"), 
    c("C", "T"), 
    c("A", "C", "G", "T")))
  apply(nuc_trios, 1, paste, collapse = "")
}



#' Extract trinucleotides from sequence
#'
#' @param x List from strsplit containing site sequences
#'
#' @return Character vector of sequence trinucleotides
.split_seq_trinucs = function(x) {
  sapply(1:(length(x)-2), function(i) paste0( x[i:(i+2)], collapse = ""))
}



#' Counts the number of each trinucleotide within each site
#'
#' @param gr_seq GenomicRanges object of sites
#'
#' @return Data frame containing each site id (character) and its trinucleotide counts (numeric)
.get_trinucs2 = function(gr_seq) {
  
  # get all potential NNN>N combos
  quadnuc = .get_quadnuc_map()
  
  # we need to extend by one nucleotide to get the neighbouring context 
  # for each nucleotide of the sequence 
  GenomicRanges::start(gr_seq) = GenomicRanges::start(gr_seq) - 1
  GenomicRanges::end(gr_seq) = GenomicRanges::end(gr_seq) + 1
  
  # obtain sequences and split these into overlapping trinucleotides
  this_seq = strsplit(as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, gr_seq)), '')
  tri_nucleotide = unlist(lapply(this_seq, .split_seq_trinucs))
  
  # map A/G-centred trinucleotides to complementary strand	
  where_reverse = substr(tri_nucleotide, 2, 2) %in% c("A", "G")
  tri_nucleotide[where_reverse] = 
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(tri_nucleotide[where_reverse])))
  
  # get trinuc mutation frequency, remove trinucs including N's 
  trinuc_dfr = data.frame(table(tri_nucleotide), stringsAsFactors = FALSE)
  trinuc_dfr = trinuc_dfr[trinuc_dfr$tri_nucleotide %in% quadnuc$trinuc,, drop = FALSE]
  
  # merge trinucleotide freq.table with all potential trinucs and quad.nucs
  # assign zeroes to trinucleotides that were not counted
  trinuc_dfr = merge(trinuc_dfr, quadnuc, 
                     by.x = "tri_nucleotide", by.y = "trinuc", all = TRUE)
  trinuc_dfr[is.na(trinuc_dfr$Freq), "Freq"] = 0			
  
  trinuc_dfr = trinuc_dfr[,c("tri_nucleotide", "Freq", "quadnuc")]
  colnames(trinuc_dfr) = c("trinuc", "posi_count", "quadnuc")
  
  # indels have all the space to occur regardless of trinucleotide context
  trinuc_dfr[trinuc_dfr$trinuc == "indel", "posi_count"] = length(tri_nucleotide)
  
  trinuc_dfr
}




#' Determines all trinucleotide, pyrimidine-centered trinucleotide contexts and mutation/quadnucleotide contexts
#'
#' @return Data frame trinucleotides, alternate allele and possible quadnucleotide contexts
.get_quadnuc_map = function() {
  quadnuc = S4Vectors::expand.grid(.get_all_trinucs(), c("A", "C", "G", "T"))
  quadnuc$tag = paste(quadnuc[,1], quadnuc[,2], sep="_")
  quadnuc = data.frame(as.matrix(quadnuc), stringsAsFactors=F)
  colnames(quadnuc) = c("trinuc", "alt", "quadnuc")
  
  # remove where ref and alt are the same
  quadnuc = quadnuc[gsub("^.(.).$", "\\1", quadnuc[,1]) != quadnuc[,2],]
  
  # indel as pseudo quadnuc
  quadnuc = rbind(quadnuc, c("indel", "indel", "indel"))
  quadnuc
}
