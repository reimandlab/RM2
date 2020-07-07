#' CTCF binding sites
#'
#' High confidence CTCF binding sites on chr3 and chr4 from ENCODE
#'
#' @name ctcf_chr3_4
#'
#' @docType data
#'
#' @usage data(ctcf_chr3)
#'
#' @format A data frame containing the following columns: chr, start, end
#' \describe{
#'     \item{chr}{chr3 and chr4 only}
#'     \item{start}{the start position of the mutation in base 0 coordinates}
#'     \item{end}{the end position of the mutation in base 0 coordinates}
#' }
#'
#' @keywords datasets
#'
#' @references ENCODE Project Consortium. "The ENCODE (ENCyclopedia Of DNA Elements) Project" Science 306(5696):636-40 (2004).
#'
#' @source \href{http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/}{FTP Server}
#'
#' @examples
#' data(ctcf_chr3_4)
#' \donttest{
#' data(mutations_chr3_4)
#' muts = cbind(mutations_chr3_4, get_mut_trinuc_strand(mutations_chr3_4))
#' window_size = 25
#' RM2(muts, ctcf_chr3_4, window_size=window_size)
#' }
