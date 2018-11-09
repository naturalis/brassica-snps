#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use My::Brassica;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

# process command line files
my $infile = '/home/ubuntu/brassica-snps/results/linkages/linkages_to_abs.tsv';
my $db     = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $range  = 100;
GetOptions(
  'infile=s' => \$infile,
  'db=s'     => \$db,
  'range=i'  => \$range,
);

# connect to database
my $schema = My::Brassica->connect("dbi:SQLite:$db");

# start reading the infile
my @header;
open my $fh, '<', $infile or die $!;
while(<$fh>) {
  chomp;
  my @line = split /\t/, $_;
  
  # process header
  if ( not @header ) {
    @header = @line;
    print join("\t", @header, 'snp_pos', 'snp_ref', 'snp_alt'), "\n";
  }
  
  # process record 
  else {
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    my $chromosome_id = substr $record{'linkage_group'}, 1;
    my @coordinates = sort { $a <=> $b } grep { $_ } @record{qw(fw_primer_start fw_primer_end rev_primer_start rev_primer_end)};
    my $snps = $schema->resultset('KaleSnp')->search({
      'chromosome_id' => $chromosome_id,
      'position'      => {
        '>=' => $coordinates[0 ] - $range,
        '<=' => $coordinates[-1] + $range,
      }
    });
    while( my $snp = $snps->next ) {
      no warnings 'uninitialized';
      print join("\t", @record{@header}, $snp->position, $snp->ref, $snp->alt), "\n";
    }
  }
}