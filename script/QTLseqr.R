#!/usr/bin/Rscript

# load the packages
library("QTLseqr")
library("dplyr")

# possible bulks:
# - group-1-EF
# - group-2-IF
# - group-3-LF
# - group-5-NF

# Set sample and file names
HighBulk <- "group-2-IF"
LowBulk <- "group-1-EF"
file <- "SNPs_from_GATK-EF-IF.table"

# Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- paste0(rep("C", 9), 1:9)

# Import SNP data from file
df <- importFromGATK(
	file = file,
	highBulk = HighBulk,
	lowBulk = LowBulk,
	chromList = Chroms
)

# Filter SNPs based on some criteria
df_filt <- filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.20,
        minTotalDepth = 100,
        maxTotalDepth = 400,
        minSampleDepth = 40,
        minGQ = 99
)

# Run G' analysis
df_filt <- runGprimeAnalysis(
	SNPset = df_filt,
	windowSize = 1e6,
	outlierFilter = "deltaSNP"
)

# write high G' snps to file
snps <- filter(df_filt, Gprime > 2.5)
write.csv(snps, file="SNPs-gprime2.5-EF-LF.csv")

# Plot G' results
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)

# Run QTLseq analysis
# XXX popStruc ought to be RIL
df_filt <- runQTLseqAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        popStruc = "RIL",
        bulkSize = c(11, 11),
        replications = 10000,
        intervals = c(95, 99)
)

# Plot deltaSNP results
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

# export summary CSV
getQTLTable(
	SNPset = df_filt, 
	alpha = 0.01, 
	export = TRUE, 
	fileName = "SNPs-gprime-filtered-EF-IF.csv"
)
