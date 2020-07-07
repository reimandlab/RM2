#' PCAWG Liver-HCC chr3 mutations
#'
#' Liver-HCC chr3 simple somatic mutations from PCAWG
#'
#' @name mutations_chr3
#'
#' @docType data
#'
#' @usage data(mutations_chr3)
#'
#' @format A data frame containing the following columns: chr, pos1, pos2, ref, alt, patient.
#' \describe{
#'     \item{chr}{chr3 only}
#'     \item{start}{the start position of the mutation in base 1 coordinates}
#'     \item{end}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#' }
#'
#' @keywords datasets
#'
#' @references The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium. "Pan-cancer analysis of whole genomes". Nature 578, pages82–93(2020).
#'
#' @source \href{https://dcc.icgc.org/releases/PCAWG}
#'
#' @examples
#' data(mutations_chr3)
#' \donttest{
#' data(ctcf_chr3)
#' muts = cbind(mutations_chr3, get_mut_trinuc_strand(mutations_chr3))
#' window_size = 25
#' RM2(muts, ctcf_chr3, window_size=window_size
#' }