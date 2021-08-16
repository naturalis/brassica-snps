#!/bin/bash

# documented in supplementary materials

# reference genome, only mapping against chromosomes
REF=${DATA}/reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa

# these are sample IDs, which we will turn into file paths below
SAMPLES="group-1-EF group-2-IF group-3-LF group-5-NF group-4"

# threads=ncores
THREADS=16

# index reference
bwa index "${REF}"

# do the mapping
for SAMPLE in $SAMPLES; do
	
	# create file stem
	STEM=${DATA}/BSA/${SAMPLE}/${SAMPLE}

	# do the mapping, include the @RG tag to identify samples when merging
	# XXX this tag is important for the subsequent GATK operations and the R script
	bwa mem \
		-R "@RG\tID:NA\tSM:${SAMPLE}\tPL:ILLUMINA\tPI:NA" \
		-t ${THREADS} \
		"${REF}" \
		"${STEM}"_R1.fastq.gz "${STEM}"_R2.fastq.gz \
		> "${STEM}"_pe.sam

	# convert to BAM
	samtools view -S -b "${STEM}"_pe.sam > "${STEM}"_pe.bam
	rm "${STEM}"_pe.sam

	# sort the reads in the BAM
	samtools sort "${STEM}"_pe.bam -o "${STEM}"_pe.sorted.bam
	rm "${STEM}"_pe.bam
done
