#!/usr/bin/perl
use strict;
use warnings;
use Vcf;
use Bio::SeqIO;
use Getopt::Long;
use Data::Dumper;
use My::Brassica;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

# process command line arguments
my $vcffile;
my $gene;
my $readgroup;
my $mincover = 0.5;
GetOptions(
  'vcffile=s'   => \$vcffile,
  'gene=s'      => \$gene,
  'readgroup=s' => \$readgroup,
  'mincover=f'  => \$mincover,
);

# initialize variables and databases
my $db     = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $ref    = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes';
my $schema = My::Brassica->connect("dbi:SQLite:$db");

# lookup CDSs
my $cdss = $schema->resultset("Feature")->search({ attributes => { LIKE => "%ID=CDS:${gene}.%" } });
CDS: while( my $cds = $cdss->next ) {
  
	# get coordinates
	my $chr    = $cds->chromosome_id;
	my $start  = $cds->feat_start;
	my $end    = $cds->feat_end;
	my $phase  = $cds->phase;
	my $strand = $cds->strand;  
  my $region = "C${chr}:${start}-${end}";
  
	# get coding, in-frame, reference sequence
	my $seq;
	eval {
		$seq = get_refseq(
			'chr'    => $chr,
			'start'  => $start,
			'stop'   => $end,
			'phase'  => $phase,
			'strand' => $strand,
		);
	};
	if ( $@ ) {
		ERROR "Problem extracting ${chr}[${strand}]:${start}..${end} (phase: ${phase})";
		ERROR $@;
		next CDS;
	}  
	my $raw = $seq->seq;
  
  # instantiate VCF reader
  my $vcf = Vcf->new(
    'file'   => $vcffile,
    'region' => $region,
  );
  $vcf->parse_header();
  my $offset = 0;
  
  # do some simple parsing. Most thorough but slowest way how to get the data.
  while( my $x = $vcf->next_data_hash() ) {

    # only continue if there is at least one alternative allele
    my @alts = grep { $_ ne '<NON_REF>' } @{ $x->{ALT} };
    if ( @alts ) {

      # get total depth and position
      my $dp  = $x->{INFO}->{DP};
      my $pos = $x->{POS};
      my $ref = $x->{REF};
      next unless $dp;
      die Dumper($x) unless $ref;

      # lookup alternative whose depth exceeds $mincover
      my @ads = split /,/, $x->{gtypes}->{$readgroup}->{AD};
      my $alt;
      for my $i ( 0 .. $#alts ) {
        if ( ( $ads[ $i + 1 ] / $dp ) > $mincover ) {
          $alt = $alts[ $i ];
        }
      }
      
      # we have an alternative allele
      if ( $alt ) {
        WARN "$alt $strand" if length($alt) > 1;
        INFO "$gene $region $pos: $alt";
        
				# adjust SNP coordinate
				my $snp_coord;
				if ( $strand eq '-' ) {
					# CDS is on '-' strand, count 
					# backward from 3' location
					$snp_coord = $end - $pos;
				}
				else {
					# CDS is on '+' strand, count
					# forward relative to CDS start
					$snp_coord = ( $pos - $start ) + $offset;
				}
				$snp_coord -= $phase; # is 0, 1 or 2      
				
				# revcom $ref & $alt if on '-' strand
				my ( $cref, $calt ) = ( $ref, $alt );
				if ( $strand eq '-' ) {
					$cref =~ tr/ACGT/TGCA/;
					$calt =~ tr/ACGT/TGCA/;
					$cref = reverse($cref);
					$calt = reverse($calt);
				}	
				
				# do sanity checks, splice snp
				my $retval = splice_snp(
					'raw' => $raw,
					'ref' => $cref,
					'alt' => $calt,
					'pos' => $snp_coord,
					'str' => $strand,
				);
				if ( $retval ) {
				  $raw = $retval;
				  $offset += ( length($calt) - 1 );
				}
      }
    }
  }
  print ">${readgroup} ${gene} ${region}\n${raw}\n";
}

sub splice_snp {
	my %args = @_;
	my $raw = $args{raw};
	
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
	
	# check if observed matches expected
	if ( $obs_ref ne $exp_ref ) {
		ERROR "Error: $obs_ref != $exp_ref (relative position: $pos, strand: ".$args{str}.")";
	}
	
	# replace the reference allele with the alternative
	substr( $raw, $pos, $length, $args{alt} );
	return $raw;
}

sub get_refseq {
	my %args  = @_;
	my ( $chr, $start, $stop ) = @args{ qw(chr start stop) };
	$chr = "C$chr" if $chr !~ /^C/; # translate chromosome foreign key to FASTA ID
	my $fasta = `fastacmd -d $ref -s $chr -L $start,$stop`;

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