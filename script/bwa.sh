#!/bin/bash

# reference genome, only mapping against chromosomes
REF=./reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa

# these are prefixes, add the suffixes below in an inner loop
SAMPLES="./BSA/group-1-EF/group-1-EF ./BSA/group-2-IF/group-2-IF ./BSA/group-3-LF/group-3-LF ./BSA/group-5-NF/group-5-NF ./gDNA/group-4/group-4"

# threads=ncores
THREADS=16

# index reference
bwa index $REF

# do the mapping
for SAMPLE in $SAMPLES; do

	# parse sample ID out of location string
	SM=`echo ${SAMPLE} | cut -f 3 -d '/'`

	# do the mapping, include the @RG tag to identify samples when merging
	# XXX this tag is important for the subsequent GATK operations and the R script
	bwa mem \
		-R "@RG\tID:NA\tSM:${SM}\tPL:ILLUMINA\tPI:NA" \
		-t ${THREADS} \
		$REF \
		${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz \
		> ${SAMPLE}_pe.sam

	# convert to BAM
	samtools view -S -b ${SAMPLE}_pe.sam > ${SAMPLE}_pe.bam
	rm ${SAMPLE}_pe.sam

	# sort the reads in the BAM
	samtools sort ${SAMPLE}_pe.bam -o ${SAMPLE}_pe.sorted.bam
	rm ${SAMPLE}_pe.bam

	# index the BAM
	samtools index ${SAMPLE}_pe.sorted.bam
done
