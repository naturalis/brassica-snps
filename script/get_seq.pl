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
GetOptions(
  'vcffile=s'   => \$vcffile,
  'intervals=s' => \$intervals,
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
    if ( scalar( @{ $x->{ALT} } ) > 1 ) {
      DEBUG Dumper($x);
    }
  }
}
