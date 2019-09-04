#!/bin/bash

CDSS=$(ls ../results/genes/*/combined.fasta)
for CDS in $CDSS; do

  # do normal alignment
  OUTFILE=${CDS//.fasta/-aligned.fasta}
  muscle -in "$CDS" > "$OUTFILE"

  # calculate proportion of invariant sites
  PINVAR=$(perl -MBio::Phylo::IO=parse -e 'print parse->[0]->calc_prop_invar' type dna format fasta file "$OUTFILE")
  echo "$PINVAR $CDS"

  # if less than 95% we probably should reverse complement
  TEST=$(perl -e "print $PINVAR < 0.95 ? 'TRUE' : 'FALSE'")
  if [ "$TEST" == "TRUE" ]; then

    # do revcom and align again
    perl ./revcom.pl -r Brassica_oleracea.v2.1.dna.toplevel -f "$CDS" -o "$CDS".revcom
    muscle -in "$CDS".revcom > "$OUTFILE"
  fi
done