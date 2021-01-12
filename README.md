# Regression Models for Localised Mutations
Regression Models for Localised Mutations (RM2) is a tool for evaluating differential mutation rates and processes across classes of functional sites. RM2 uses negative binomial regression models to assess whether sites of interest show an enrichment or depletion of total mutations compared to flanking control regions of the same length. Further, RM2 was designed to test whether subclasses of mutations, like those of specific signatures, strandedness, transcription direction, and other features show differential patterns. This allows for quick and systematic characterization of mutational processes and can be easily extended to site-based studies like pathway analysis.   

As input, RM2 requires:
1. Mutation file with optional annotation columns
2. Set of genomic regions containing the chr, start and end position of the sites

## Installation
#### devtools:
Using the R package devtools, run: devtools::install_github('https://github.com/reimandlab/RM2')

#### From source
Clone the repository: https://github.com/reimandlab/RM2.git <br />
Open R in the directory you cloned the package in and run install.packages('RM2', repos=NULL)

## Usage
library(RM2)

#### Load sites
data("ctcf_chr3_4")

#### Load mutations
data("mutations_chr3_4")

#### Add mutation annotations
muts = cbind(mutations_chr3_4, get_mut_trinuc_strand(mutations_chr3_4))

#### Run regression
window_size = 50 <br />
results = RM2(muts, ctcf_chr3_4, window_size = window_size) <br />
results

#### Visualize results
dfr = get_mutations_in_flanked_sites(muts, ctcf_chr3_4, window_size) <br />
plot_mutations_in_flanked_sites(dfr, window_size)

#### Run regression with additional mutation subclasses
mut_class_columns = c("mut_strand", "ref_alt") <br />
results = RM2(muts, ctcf_chr3_4, window_size = window_size, mut_class_columns = mut_class_columns) <br />
results

#### Run regression with additional mutation co-factor
muts$cofactor_col = sample(c(0,1), nrow(muts), replace=T) <br />
results = RM2(muts, sites = ctcf_chr3_4, window_size = window_size, cofactor_column = "cofactor_col") <br />
results
