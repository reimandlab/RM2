% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preparing_mutations.R
\name{get_mut_trinuc_strand}
\alias{get_mut_trinuc_strand}
\title{Derive quadnucleotide contexts, strandedness and single-base substitutions}
\usage{
get_mut_trinuc_strand(maf)
}
\arguments{
\item{maf}{Data frame of mutations with the following columns: chr, start, end, ref and alt
\describe{
\item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
\item{start}{the start position of the mutation in base 1 coordinates}
\item{end}{the end position of the mutation in base 1 coordinates}
\item{ref}{the reference allele as a string containing the bases A, T, C or G}
\item{alt}{the alternate allele as a string containing the bases A, T, C or G}
}}
}
\value{
Data frame containing additional mutation context columns
\describe{
\item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
\item{start}{the start position of the mutation in base 1 coordinates}
\item{end}{the end position of the mutation in base 1 coordinates}
\item{ref}{the reference allele as a string containing the bases A, T, C or G}
\item{alt}{the alternate allele as a string containing the bases A, T, C or G}
\item{mut_trinuc}{pyrimidine-centered quadnucleotide mutation context and indel}
\item{mut_strand}{the mutation strand classified as w for Watson and c for Crick}
\item{ref_alt}{pyrimidine single-base substitution}
}
}
\description{
get_mut_trinuc_strand() takes a data frame of mutations and determines the quadnucleotide context (trinucleotide context + alternate allele) and strand.
}
