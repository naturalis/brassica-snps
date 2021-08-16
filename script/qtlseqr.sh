#!/bin/bash

# documented in supplementary materials

# all pairwise comparisons between phenotypes
CONTRASTS="EF-IF EF-LF EF-NF IF-LF IF-NF LF-NF"

# iterate over pairwise comparisons
for C in ${CONTRASTS}; do

	# TSV table to be used by QTLseqr.R
	OUT="${DATA}/contrasts/${C}/SNPs_from_GATK.table"

	# previously computed, joint genotypes
	VCF="${DATA}/contrasts/${C}/joint-genotypes.vcf"

	# create table
	if [ ! -e "${OUT}" ]; then
		gatk VariantsToTable \
			-F CHROM \
			-F POS \
			-F REF \
			-F ALT \
			-GF AD \
			-GF DP \
			-GF PL \
			-GF GQ \
			--variant "${VCF}" \
			--output "${OUT}"
  fi
done
