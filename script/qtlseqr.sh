#!/bin/bash

# all pairwise comparisons between phenotypes
CONTRASTS="EF-IF EF-LF EF-NF IF-LF IF-NF LF-NF"

# iterate over pairwise comparisons
for CONTRAST in $CONTRASTS; do 

	# table to be used by QTLseqr.R
	OUT=$DATA/contrasts/$CONTRAST/SNPs_from_GATK.table

	# computed, joint genotypes
	VCF=$DATA/contrasts/$CONTRAST/joint-genotypes.vcf

	# create table
	if [ ! -e ${OUT} ]; then
		gatk VariantsToTable \
			-F CHROM \
			-F POS \
			-F REF \
			-F ALT \
			-GF AD \
			-GF DP \
			-GF PL \
			-GF GQ \
			--variant ${VCF} \
			--output ${OUT}
fi

done
