#!/usr/bin/Rscript

# load the package
library("QTLseqr")

# Set sample and file names
HighBulk <- "group-3-LF"
LowBulk <- "group-1-EF"
file <- "BSA/SNPs_from_GATK.table"

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

# Plot G' results
pdf("BSA/Gprime.pdf")
print(plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01))
dev.off()

# Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        popStruc = "F2",
        bulkSize = c(11, 11),
        replications = 10000,
        intervals = c(95, 99)
)

# Plot deltaSNP results
pdf("BSA/deltaSNP.pdf")
print(plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE))
dev.off()

# export summary CSV
getQTLTable(
	SNPset = df_filt, 
	alpha = 0.01, 
	export = TRUE, 
	fileName = "BSA/my_BSA_QTL.csv"
)
