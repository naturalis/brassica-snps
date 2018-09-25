#!/usr/bin/perl
use strict;
use warnings;
use Vcf;
use Getopt::Long;
use Data::Dumper;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

# process command line arguments
my $vcffile;
my $intervals;
my $readgroup;
my $mincover = 0.5;
GetOptions(
  'vcffile=s'   => \$vcffile,
  'intervals=s' => \$intervals,
  'readgroup=s' => \$readgroup,
  'mincover=f'  => \$mincover,
);

# open intervals file, iterate over lines
open my $fh, '<', $intervals or die $!;
while(<$fh>) {
  
  # parse focal interval
  chomp;
  my $line = $_;
  
  # instantiate VCF reader
  my $vcf = Vcf->new(
    'file'   => $vcffile,
    'region' => $line,
  );
  $vcf->parse_header();
  
  # do some simple parsing. Most thorough but slowest way how to get the data.
  while( my $x = $vcf->next_data_hash() ) {

    # only continue if there is at least one alternative allele
    my @alts = grep { $_ ne '<NON_REF>' } @{ $x->{ALT} };
    if ( @alts ) {

      # get total depth and position
      my $dp  = $x->{INFO}->{DP};
      my $pos = $x->{POS};
      next unless $dp;

      # lookup allelic depths
      my @ads = split /,/, $x->{gtypes}->{$readgroup}->{AD};
      my $string;
      for my $i ( 0 .. $#alts ) {
        if ( ( $ads[ $i + 1 ] / $dp ) > $mincover ) {
          $string .= $alts[ $i ] . ':' . $ads[ $i + 1 ] . '/' . $dp . ' ';
        }
      }
      if ( $string ) {
        print $line, "\t", $pos, "\t", $string, "\n";
      }
    }
  }
}
