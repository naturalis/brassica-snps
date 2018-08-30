#!/bin/bash

# reference genome
REF=./reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa

# SNPs from bulk of early flowering plants
RG1=./BSA/group-1-EF/group-1-EF_pe.sorted.bam.RG.vcf

# SNPs from bulk of late flowering plants
RG2=./BSA/group-3-LF/group-3-LF_pe.sorted.bam.RG.vcf

# location for the output file with the combined SNPs
OUT1=./BSA/combined-snps.vcf

# location for the computed genotypes
OUT2=./BSA/joint-genotypes.vcf

# combine the separate SNP files
gatk CombineGVCFs \
        --reference ${REF} \
        --variant ${RG1} \
        --variant ${RG2} \
        --output ${OUT1}

# do the genotyping
gatk GenotypeGVCFs \
	--reference ${REF} \
	--variant ${OUT1} \
	--output ${OUT2}

