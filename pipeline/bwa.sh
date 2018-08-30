#!/bin/bash

# reference genome, only mapping against chromosomes
# REF=./reference/Brassica_oleracea.v2.1.dna.toplevel.fa.gz
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

	# deduplicate the reads: unzip, dedup, rezip
#	gunzip ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
#	fastuniq -i ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq -t q -o ${SAMPLE}_R1.dedup.fastq -p ${SAMPLE}_R2.dedup.fastq
#	gzip -9 ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq ${SAMPLE}_R1.dedup.fastq ${SAMPLE}_R2.dedup.fastq
#	ParDRe \
#		-i ${SAMPLE}_R1.fastq.gz \
#		-p ${SAMPLE}_R2.fastq.gz \
#		-o ${SAMPLE}_R1.fastq.dedup.gz \
#		-r ${SAMPLE}_R2.fastq.dedup.gz \
#		-t 16 \
#		-z 9

	# do the mapping
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
