#!/usr/bin/perl
use strict;
use warnings;
use List::Util 'sum';
my @header;
my $snp_id = 1;
LINE: while(<>) {
  chomp;
  my @line = split /\t/, $_;
  
  # process header
  if ( not @header ) {
    for ( @line ) {
      s/^group-4\.//;
      push @header, $_;
    }
    
    # no header, so we can import into SQLite
    # print join("\t", @header), "\n";
  }
  
  # process record
  else {
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    
    # phred-scaled genotype quality
    next LINE if $record{'GQ'} < 99;
    
    # must be true *S*NP: REF and ALT can only be length 1
    next LINE if length($record{'REF'}) != 1 or length($record{'ALT'}) != 1;
    
    # avoid low coverage (<100) or repeats (>400)
    my @cover = split /,/, $record{'AD'};
    next LINE if sum(@cover) < 100 or sum(@cover) > 400;
    
    # rightmost genotype, 1/1, has highest score, i.e. homozygous alternative SNP
    my @genotypes = split /,/, $record{'PL'};
    next LINE if $genotypes[2] != 0; 
    
    # change chromosome name to ID
    $record{'CHROM'} =~ s/^C//;
    
    # print remaining SNPs
    print join("\t", $snp_id++, @record{qw(CHROM POS REF ALT)}), "\n";
  }
}