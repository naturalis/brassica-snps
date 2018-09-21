#!/bin/bash
INTERVALS=$1
DATA=/home/ubuntu/data
REF=$DATA/reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa
BSA=$DATA/BSA
BULKS="group-1-EF group-2-IF group-3-LF group-5-NF"
SUFFIX=_pe.sorted.bam.RG.vcf.gz

for BULK in $BULKS; do
  java -jar /usr/local/src/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
    -R ${REF} \
    -o ${BULK}.fasta \
    -L ${INTERVALS} \
    -V ${BSA}/${BULK}/${BULK}${SUFFIX}
done
