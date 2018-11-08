#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw(tempfile);
use Bio::Tools::Run::StandAloneBlastPlus;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($ERROR);

# process command line arguments
my $infile    = '/home/ubuntu/data/reference/12864_2012_4560_MOESM1_ESM.tsv';
my $reference = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes.fa';
GetOptions(
  'infile=s'    => \$infile,
  'reference=s' => \$reference,
);

# instantiate BLAST object
my $bp = Bio::Tools::Run::StandAloneBlastPlus->new( '-db_name' => $reference );

# start reading the file
my @header;
open my $fh, '<', $infile or die $!;
while(<$fh>) {
  chomp;
  my @line = split /\t/, $_;
  
  # process header
  if ( not @header ) {
    @header = @line;
  }
  
  # process record
  else {
    
    # populate the record
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    
    # run blast
    my @results;
    for ( qw(fw_primer rev_primer) ) {
      
      # do logging
      my $template = 'going to blast the %s for marker number %d (name: %s, type: %s) on chromosome %s';
      DEBUG sprintf($template, @record{($_, qw(position marker_name marker_type linkage_group))});
      push @results, do_blast(%record, 'primer_type' => $_);
    }
    
  }
}

sub do_blast {
  my %args = @_;
  
  # write input file
  my ( $fh, $filename ) = tempfile();
  print $fh ">query\n", $args{$args{'primer_type'}};
  close $fh;
  
  # run query
  my $result = $bp->blastn( '-query' => $filename );
  my @hits;
	while( my $hit = $result->next_hit ) {
		if ( $hit->accession eq $args{'linkage_group'} ) {
		  DEBUG "identity: " . $hit->frac_identical('total');
		  print Dumper($hit);
		}
	}
  
  unlink $filename;
}