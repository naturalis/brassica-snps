#!/usr/bin/perl
use strict; 
use warnings; 
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Getopt::Long; 
use My::Brassica;
use Data::Dumper;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

# process command line arguments
my $id;
my $contrast;
GetOptions(
  'id=s'       => \$id,
  'contrast=s' => \$contrast,
);

# initialize variables and databases
my $db     = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $ref    = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes';
my $schema = My::Brassica->connect("dbi:SQLite:$db");
my $ctable = Bio::Tools::CodonTable->new();

# select * from features where attributes like '%ID=CDS:Bo6g103730%';
my $cdss = $schema->resultset("Feature")->search({ attributes => { LIKE => "%ID=CDS:$id.%" } });

# store reference and alternative sequence segments here
my ( @ref, @alt );

# iterate over CDS features
while( my $cds = $cdss->next ) {

  # get coordinates
  my $chr    = $cds->chromosome_id;
  my $start  = $cds->feat_start;
  my $end    = $cds->feat_end;
  my $phase  = $cds->phase;
  my $strand = $cds->strand;
  
  # get coding, in-frame, reference sequence
  my $seq = get_refseq(
  	'chr'    => $chr,
  	'start'  => $start,
  	'stop'   => $end,
  	'phase'  => $phase,
  	'strand' => $strand,
  )->seq;
  push @ref, $seq;
  
  # search SNPs
  my $snps = $schema->resultset("Snp")->search({
  	chromosome_id => $chr,
    	contrast      => $contrast,
  	position      => { '>=' => $start, '<=' => $end },
  });
  
  # iterate over SNPs
  while( my $snp = $snps->next ) {
  
    # lookup variables
    my $pos = $snp->position;
    my $ref = $snp->ref;
    my $alt = $snp->alt;
    DEBUG "$chr:$pos $ref -> $alt";
  
    # adjust SNP coordinate
    my $snp_coord;
    if ( $strand eq '-' ) {
      
      # CDS is on '-' strand, count backward from 3' location
      $snp_coord = $end - $pos;
    }
    else {
      
      # CDS is on '+' strand, count forward relative to CDS start
      $snp_coord = $pos - $start;
    }
    $snp_coord -= $phase; # subtract, might be 0, 1 or 2
  					
    # revcom $ref & $alt if on '-' strand
    my ( $cref, $calt ) = ( $ref, $alt );
    if ( $strand eq '-' ) {
      $cref =~ tr/ACGT/TGCA/;
      $calt =~ tr/ACGT/TGCA/;
      $cref = reverse($cref);
      $calt = reverse($calt);
    }
  
    # splice alternative allele					
    my $retval = splice_snp(
      'seq' => $seq,
      'ref' => $cref,
      'alt' => $calt,
      'pos' => $snp_coord,
      'str' => $strand,
    );
  					
    # replace $seq with $retval if any
    $seq = $retval if $retval;
  }
  
  # store spliced seq
  push @alt, $seq;
}

# print output
print ">REF\n", join('',@ref), "\n";
print ">$contrast\n", join('',@alt), "\n";

sub splice_snp {
	my %args = @_;
	my $raw = $args{seq};
	
	# Compute starting position differently depending on forward or reverse 
	# strand: we need to go upstream for the the reverse strand.
	my $pos = $args{pos};
	if ( $args{str} eq '-' ) {
	
		# on the '-' strand the start coordinate for alleles needs to
		# be adjusted, but using 0-based indexing. I.e. this has no
		# effect for alleles of 1bp length, which is nearly all of them.
		my $l = length($args{ref}) - 1;
		$pos -= $l;
	}
	
	# allele cannot be nonsynonymous if we've moved upstream outside of the phased CDS
	return if $pos < 0;
	
	# Compute length so that we don't go outside of the phased CDS
	my $length = length($args{ref});
	if ( ( $length + $pos ) > length($raw) ) {
		$length = length($raw) - $pos;
	}
	
	# extract observed allele
	my $obs_ref = substr( $raw, $pos, $length ); 
	my $exp_ref = substr( $args{ref}, 0, $length );
	
	# allele cannot be nonsynonymous if we've moved downstream outside of the phased CDS
	return unless $obs_ref;
	
	# replace the reference allele with the alternative
	substr( $raw, $pos, $length, $args{alt} );
	return $raw;
}

sub get_refseq {
	my %args  = @_;
	my ( $chr, $start, $stop ) = @args{ qw(chr start stop) };
	$chr = "C$chr" if $chr !~ /^C/; # translate chromosome foreign key to FASTA ID
	my $fasta = `fastacmd -d $ref -s $chr -L $start,$stop`;
	DEBUG $fasta;
	
	# parse raw FASTA
	my $seq = Bio::SeqIO->new( 
		'-string' => $fasta, 
		'-format' => 'fasta',
	)->next_seq;
	
	# reverse complement if CDS is on '-' strand
	$seq = $seq->revcom if $args{'strand'} eq '-';
	
	# truncate sequence if there is a phase offset
	$seq = $seq->trunc( $args{'phase'} + 1, $args{'stop'} - $args{'start'} ) if $args{'phase'};
	return $seq;
}
