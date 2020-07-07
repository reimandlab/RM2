# Regression Models for Localized Mutations
Regression Models for Localized Mutations (RM2) is a tool for evaluating differential mutation rates and processes across classes of functional sites. RM2 uses negative binomial regression models to assess whether sites of interest show an enrichment or depletion of total mutations compared to flanking control regions of the same length. Further, RM2 was designed to test whether subclasses of mutations, like those of specific signatures, strandedness, transcription direction, and other features show differential patterns. This allows for quick and systematic characterization of mutational processes and can be easily extended to site-based studies like pathway analysis.   

As input, RM2 requires: 1) mutation file with optional annotation columns, and 2) set of genomic regions containing the chr, start and end position of the sites.

## Installation
#### devtools:
Using the R package devtools, run devtools::install_github('https://github.com/reimandlab/RM2')

#### From source
Clone the repository: https://github.com/reimandlab/RM2.git Open R in the directory you cloned the package in and run install.packages('RM2', repos=NULL)

## Usage
library(RM2)

#### load sites
data("ctcf_chr3_4")

#### load mutations
data("mutations_chr3_4")

#### add mutation annotations
muts = cbind(mutations_chr3_4, get_mut_trinuc_strand(mutations_chr3_4))

#### run regression
window_size = 50 <br />
results = RM2(muts, ctcf_chr3_4, window_size=window_size)
results

#### visualize results
dfr = get_mutations_in_flanked_sites(muts, ctcf_chr3_4, window_size) <br />
plot_mutations_in_flanked_sites(dfr, window_size)
