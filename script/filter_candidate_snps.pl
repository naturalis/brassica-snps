#!/usr/bin/perl
use strict;
use warnings;
use List::Util 'sum';
my @header;
LINE: while(<>) {
  chomp;
  my @line = split /\t/, $_;
  
  # process header
  if ( not @header ) {
    for ( @line ) {
      s/^group-4\.//;
      push @header, $_;
    }
    print join("\t", @header), "\n";
  }
  
  # process record
  else {
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    
    # phred-scaled genotype quality
    next LINE if $record{'GQ'} < 99;
    
    # must be true *S*NP: REF and ALT can only be length 1
    next LINE if length($record{'REF'}) != 1 or length($record{'ALT'}) != 1;
    
    # avoid low coverage (<100) or repeats (>400)
    my @cover = split /,/, $record{AD};
    next LINE if sum(@cover) < 100 or sum(@cover) > 400;
    
    # middle genotype, 0/1, has highest score, i.e. heterozygous SNP
    my @genotypes = split /,/, $record{PL};
    next LINE if $genotypes[1] != 0; 
    
    # print remaining SNPs
    print join("\t", @record{@header}), "\n";
  }
}