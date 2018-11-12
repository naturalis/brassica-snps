#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw(tempfile);
use Bio::Tools::Run::StandAloneBlastPlus;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

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
LINE: while(<$fh>) {
  chomp;
  my @line = split /\t/, $_;
  
  # process header
  if ( not @header ) {
    @header = @line;
    push @header, qw(fw_primer_start fw_primer_end rev_primer_start rev_primer_end);
    print join("\t", @header), "\n";
  }
  
  # process record
  else {
    
    # populate the record
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    # next LINE if $record{'marker_type'} ne 'SNP';
    

    # run blast
    my $have_hits = 0;
    for ( qw(fw_primer rev_primer) ) {
      
      # do logging
      my $template = 'going to blast the %s for marker number %d (name: %s, type: %s) on chromosome %s';
      DEBUG sprintf($template, $_, @record{qw(position marker_name marker_type linkage_group)});
      my @results = do_blast(%record, 'primer_type' => $_);
      if ( @results == 1 ) {
        $record{$_ . '_start'} = $results[0]->{'start'};
        $record{$_ . '_end'}   = $results[0]->{'end'};
        $have_hits++;
      }
    }
    if ( $have_hits ) {
      no warnings 'uninitialized';
      print join("\t", @record{@header}),"\n";
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
  my $result = $bp->blastn( 
    '-query'     => $filename,
    '-method_args' => [
      '-evalue'    => 1000,
      '-word_size' => 7,
      '-gapopen'   => 5,
      '-gapextend' => 2,
      '-reward'    => 1,
      '-penalty'   => -3,
    ]
  );
  my @hits;
  HIT: while( my $hit = $result->next_hit ) {

    # only want hits on the same chromosome
    if ( $hit->accession eq $args{'linkage_group'} ) {
      my @hsps;
      HSP: while( my $hsp = $hit->next_hsp ) {
        next HSP unless $hsp->{'HSP_LENGTH'} == $hsp->{'QUERY_LENGTH'};
        next HSP unless ( $hsp->{'HSP_LENGTH'} - $hsp->{'IDENTICAL'} ) <= 1;
        push @hsps, $hsp;
      }
      next HIT if scalar(@hsps) > 1;
      push @hits, {
        'start' => $hsp->{'HIT_START'},
        'end'   => $hsp->{'HIT_END'},
      }
    }
  }
  unlink $filename;
  return @hits;
}
