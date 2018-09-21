#!/bin/bash
INTERVALS=$1
DATA=/home/ubuntu/data
REF=$DATA/reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa
BSA=$DATA/BSA
BULKS="group-1-EF group-2-IF group-3-LF group-5-NF"
SUFFIX=_pe.sorted.bam.RG.vcf
OUTFILE=merged.fasta
ALIGNED=merged.aligned.fasta

# iterate over bulk assays
for BULK in $BULKS; do

  # create focal bulk's consensus sequence from SNPs and reference
  java -jar /usr/local/src/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
    -R ${REF} \
    -o ${BULK}.fasta \
    -L ${INTERVALS} \
    -V ${BSA}/${BULK}/${BULK}${SUFFIX}
  
  # append focal bulk's name as FASTA header
  echo ">${BULK}" >> ${OUTFILE}
  
  # append all CDSs concatenated below the header
  grep -v '>' ${BULK}.fasta >> ${OUTFILE}
done

# align
muscle -in ${OUTFILE} -out ${ALIGNED}
