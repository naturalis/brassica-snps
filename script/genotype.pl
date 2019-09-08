#!/usr/bin/perl
use strict;
use warnings;

# make reference genome location and list of samples
my $DATA = $ENV{'DATA'};
my $REF  = $DATA . '/reference/Brassica_oleracea.v2.1.dna.toplevel.chromosomes.fa';
my @SAM  = qw[group-1-EF group-2-IF group-3-LF group-5-NF];

# make all permutations
for my $i ( 0 .. $#SAM-1 ) {
  for my $j ( $i+1 .. $#SAM ) {

    # dereference focal samples
    my ( $S1, $S2 ) = @SAM[$i,$j];

    # create input file paths
    my $RG1 = "${DATA}/BSA/${S1}/${S1}_pe.sorted.bam.RG.vcf";
    my $RG2 = "${DATA}/BSA/${S2}/${S2}_pe.sorted.bam.RG.vcf";

    # create phenotypes
    my ( $P1, $P2 ) = ( $S1, $S2 );
    $P1 =~ s/group-\d-//;
    $P2 =~ s/group-\d-//;

    # output files
    my $OUT1 = "${DATA}/contrasts/${P1}-${P2}/combined-snps-${P1}-${P2}.vcf";
    my $OUT2 = "${DATA}/contrasts/${P1}-${P2}/joint-genotypes-${P1}-${P2}.vcf";

    # combine the read groups
    system("gatk CombineGVCFs --reference ${REF} --variant ${RG1} --variant ${RG2} --output ${OUT1}");

    # do the joint genotyping
    system("gatk GenotypeGVCFs --reference ${REF} --variant ${OUT1} --output ${OUT2}");
  }
}