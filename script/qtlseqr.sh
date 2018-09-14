#!/bin/bash

# table to be used by QTLseqr.R
OUT=BSA/SNPs_from_GATK.table

# computed, joint genotypes
VCF=BSA/joint-genotypes.vcf

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

# now execute the R script
./QTLseqr.R
