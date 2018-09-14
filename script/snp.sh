#!/bin/bash

# reference genome
REF=./reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa

# mapped genomes as sorted BAM files
BAMS="./BSA/group-1-EF/group-1-EF_pe.sorted.bam ./BSA/group-2-IF/group-2-IF_pe.sorted.bam ./BSA/group-3-LF/group-3-LF_pe.sorted.bam ./BSA/group-5-NF/group-5-NF_pe.sorted.bam ./gDNA/group-4/group-4_pe.sorted.bam"


# iterate over BAM files
for BAM in $BAMS; do

	echo "operating on " ${BAM}

	# index the BAM file
	if [ ! -f ${BAM}.bai ]; then
		samtools index ${BAM}
	fi

	# do the haplotype calling
	if [ ! -f ${BAM}.RG.vcf ]; then
		gatk HaplotypeCaller \
			--reference ${REF} \
			--input ${BAM} \
			--output ${BAM}.RG.vcf \
			--emit-ref-confidence GVCF
	fi
done

