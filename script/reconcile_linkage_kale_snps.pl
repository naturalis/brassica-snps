#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use My::Brassica;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

# process command line files
my $infile   = '/home/ubuntu/brassica-snps/results/linkages/linkages_to_abs.tsv';
my $db       = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $range    = 0;
my $offset   = 100;
GetOptions(
  'infile=s' => \$infile,
  'db=s'     => \$db,
  'range=i'  => \$range,
  'offset=i' => \$offset,
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
    print join("\t", @header, 'snp_pos', 'snp_ref', 'snp_alt', 'window_size'), "\n";
  }
  
  # process record 
  else {
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    my $chromosome_id = substr $record{'linkage_group'}, 1;
    my @coordinates = sort { $a <=> $b } grep { $_ } @record{qw(fw_primer_start fw_primer_end rev_primer_start rev_primer_end)};
    my ( $start, $stop ) = ( $coordinates[0], $coordinates[-1] );    
    
    # look for @filtered snps while growing the window
    my $window = $range;
    my @filtered;
    while( not @filtered ) {
      @filtered = get_snps(
        'chrom' => $chromosome_id,
        'start' => $start,
        'stop'  => $stop,
        'range' => $window,
      );
      $window += 100;
    }
    
    # write output
    for my $f ( @filtered ) {
      no warnings 'uninitialized';
      print join("\t", @record{@header}, $f->position, $f->ref, $f->alt, $window), "\n";
    }
  }
}

sub get_snps {
  my %args = @_;
  
  # get all the SNPs within the coordinate range
  my @snps = $schema->resultset('KaleSnp')->search(
  {
    'chromosome_id' => $args{'chrom'},
    'position'      => {
      '>=' => $args{'start'} - $args{'range'},
      '<=' => $args{'stop'}  + $args{'range'},
    }
  },
  { 'order_by' => { '-asc' => 'position' } })->all;  
  
  # compute the distance to next and previous
  my @intermediates;
  for my $i ( 1 .. $#snps - 1 ) {
    push @intermediates, {
      'snp'      => $snps[$i],
      'previous' => ( $snps[$i  ]->position - $snps[$i-1]->position ),
      'next'     => ( $snps[$i+1]->position - $snps[$i  ]->position ),
    }
  }
  
  # retain snps whose distance to next and previous is >= $offset
  my @filtered = map { $_->{'snp'} } grep { $_->{'previous'} >= $offset && $_->{'next'} >= $offset } @intermediates;  
  return @filtered;
}