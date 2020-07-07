Regression Models for Local Mutations (RM2) is a tool for evaluating differential mutation rates and processes across classes of functional sites. RM2 uses negative binomial regression models to evaluate whether sites of interest show an enrichment or depletion of total mutations, or specific signatures, strandedness, transcription direction, and other features, compared to their aggregated flanking regions of the same length. This allows for quick and systematic characterization of mutational processes and can be easily extended to site-based studies like pathway analysis.   

As input, RM2 requires: 1) mutation file with optional annotation columns, and 2) set of genomic regions containing the chr, start and end position of the sites.

Installation
devtools:
Using the R package devtools, run devtools::install_github('https://github.com/reimandlab/RM2')

From source
Clone the repository: https://github.com/reimandlab/RM2.git Open R in the directory you cloned the package in and run install.packages('RM2', repos=NULL)

Example
library(RM2final)

# load sites
sites = read.delim("ENCODE_CTCF_sites_chr3.txt.gz", stringsAsFactors=F, header=F)
colnames(sites) = c("chr", "start", "end")

# load mutations
muts = read.delim("PCAWG_Liver-HCC_mutations_chr3.txt.gz", stringsAsFactors=F, header=F)
colnames(muts) = c("chr", "start", "end", "ref", "alt")

# add mutation annotations
muts = cbind(muts, get_mut_trinuc_strand(muts))

# run regression
window_size = 25
results = RM2(muts, sites, window_size=window_size)

# visualize results
dfr = get_mutations_in_flanked_sites(muts, sites, window_size)
plot_mutations_in_flanked_sites(dfr, window_size)
