#!/usr/bin/Rscript

# documented in supplementary materials

# load the packages. Get QTLseqr using:
# library(devtools)
# install_github("bmansfeld/QTLseqr")
library(QTLseqr)
library(dplyr)
library(ggplot2)

# possible bulks:
# - group-1-EF - 11 accessions
# - group-2-IF - 8 accessions
# - group-3-LF - 11 accessions
# - group-5-NF - 13 accessions

# Set sample and file names
HighBulk <- "group-5-NF"
LowBulk <- "group-3-LF"
file <- "SNPs_from_GATK-LF-NF.table"

# Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- paste0(rep("C", 9), 1:9)

# Import SNP data from file
df <- importFromGATK(
	file = file,
	highBulk = HighBulk,
	lowBulk = LowBulk,
	chromList = Chroms
)

# Optional: plot the read depth distribution
# XXX the coverage *range* looks very similar to that in the vignette, although they have more
# SNPs with low coverage, i.e. more skew towards the origin. Our data is actually a bit better.
# ggplot(data = df) + 
#  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
#  xlim(0,1000)

# Optional: plot total reference allele frequency
# XXX this will look vastly different from the QTLseqr vignette because this is an
# F1, so each reference allele will be present at least 50% (from the one, homozygous, T1000 parent) +
# the other parent (Jersey kale), which might also have the allele. The average would be closer to 75%.
# This would be of interest to come up with a value for refAlleleFreq, however, the default of 0.2 looks
# reasonable. In fact, somewhat conservative due to the skewness in our distribution.
# ggplot(data = df) +
#  geom_histogram(aes(x = REF_FRQ))

# Optional: plot per-bulk SNP-index
# XXX as per the vignette, we do indeed see "two small peaks on each end", but the rest of the distribution
# is not approximately normal because this is an F1
# ggplot(data = df) +
#  geom_histogram(aes(x = SNPindex.LOW))

# Filter SNPs based on some criteria
df_filt <- filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.20, # default, is conservative because F1 (see above)
        minTotalDepth = 100, # default, at least 100x coverage, i.e. well supported SNPs
        maxTotalDepth = 400, # default, no more than 400x coverage, i.e. avoid stacked repeats
        minSampleDepth = 40, # default, per sample
        minGQ = 99, # GATK's Genotype Quality score
        verbose = TRUE
)

# Run G' analysis
df_filt <- runGprimeAnalysis(
	SNPset = df_filt,
	windowSize = 1e6,
	outlierFilter = "deltaSNP"
)

# get threshold (cribbed out of plotQTLStats)
fdrT <- getFDRThreshold(df_filt$pvalue, alpha = 0.01)
logFdrT <- -log10(fdrT)
GprimeT <- df_filt[which(df_filt$pvalue == fdrT), "Gprime"]

# write high G' snps to file
snps <- filter(df_filt, Gprime > 2.5)
write.csv(snps, file="SNPs-gprime2.5-EF-IF.csv")

# Plot G' results
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)

# Run QTLseq analysis
# XXX popStruc ought to be RIL
df_filt <- runQTLseqAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        popStruc = "RIL",
        bulkSize = c(10), # XXX see top of this file for bulk sizes. Average for RIL.
        replications = 10000,
        intervals = c(95, 99)
)

# export summary CSV
getQTLTable(
	SNPset = df_filt, 
	alpha = 0.01, 
	export = TRUE, 
	fileName = "QTL-regions-EF-IF.csv"
)
