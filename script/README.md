# Contains pipeline and data exploration scripts

1. [bwa.sh](bwa.sh) - map samples against reference genome
2. [snp.sh](snp.sh) - compute all variants for each sample
3. [genotype.sh](genotype.sh) - joint genotyping over two merged samples
4. [qtlseqr.sh](qtlseqr.sh) - prepare joint genotyping results for input in QTLSeqr.R
5. [QTLseqr.R](QTLseqr.R) - compute QTL regions and G' values for SNPs
6. [genes_in_qtl_regions.pl](genes_in_qtl_regions.pl) - get gene IDs for QTL regions
7. [snps_in_cds.pl](snps_in_cds.pl) - get SNPs for gene IDs
8. [biomart.pl](biomart.pl) - get UniProtKB IDs for gene IDs
9. [go_filter.pl](go_filter.pl) - filter enrichment test results
10. [qtl_regions_for_gene.pl](qtl_regions_for_gene.pl) - get QTL regions in which the gene ID occurs
11. [get_intervals.pl](get_intervals.pl) - get genomic coordinate(s) for the stable feature ID
