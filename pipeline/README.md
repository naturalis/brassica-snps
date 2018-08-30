This directory contains the shell scripts that were executed in the course of the overall
pipeline. In order, these are:

1. bwa.sh - do a paired-end mapping of the different samples agains the reference genome,
   then convert the SAM output to BAM, sort the reads in the BAM file, then index it. 
2. snp.sh - run the gatk HaplotypeCaller to compute all variants for each sample
3. genotype.sh - merge two samples, and then do a joint genotyping over the merge
4. qtlseqr.sh - prepare the joint genotyping results for input in QTLSeqr.R
