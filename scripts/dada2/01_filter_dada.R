#!/usr/bin/Rscript
#install.packages("devtools",repos='http://cran.us.r-project.org')
#library(devtools)
#devtools::install_github("benjjneb/dada2")
suppressPackageStartupMessages(library(dada2))
sessionInfo()
track.filt <- filterAndTrim(snakemake@input[['r1']],snakemake@output[['r1']], 
                            snakemake@input[['r2']],snakemake@output[['r2']], 
                            maxEE=snakemake@config[["maxEE"]], 
                            compress=TRUE,
                            verbose=TRUE,
                            multithread=snakemake@threads)
row.names(track.filt) <- snakemake@params[["samples"]]
colnames(track.filt) = c('raw','filtered')
write.table(track.filt,snakemake@output[['nreads']],  sep='\t')