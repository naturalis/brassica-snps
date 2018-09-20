#!/usr/bin/perl
use strict; 
use warnings; 
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Getopt::Long; 
use My::Brassica;
use Data::Dumper;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($ERROR);

=pod

Given a list of EnsEMBL gene IDs, provided on STDIN line by line, looks for SNPs within
CDSs of these genes. The output is a tab-separated table with the following columns:

- gene ID
- chromosome number
- start coordinate of the CDS
- stop coordinate of the CDS
- codon phase
- strand
- SNP position
- syn/nonsyn
- reference allele
- alternative allele
- contrasts in which the SNP occurs

=cut

# process command line arguments
my @genes;
while(<>) {
	chomp;
	push @genes, $_ if /\S/;
}

# initialize variables and databases
my $db     = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $ref    = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes';
my $schema = My::Brassica->connect("dbi:SQLite:$db");
my $ctable = Bio::Tools::CodonTable->new();

# iterate over gene IDs
for my $id (@genes) {

	# select * from features where attributes like '%ID=CDS:Bo6g103730%';
	my $cdss = $schema->resultset("Feature")->search({ attributes => { LIKE => "%ID=CDS:$id.%" } });

	# iterate over CDS features
	CDS: while( my $cds = $cdss->next ) {

		# get coordinates
		my $chr    = $cds->chromosome_id;
		my $start  = $cds->feat_start;
		my $end    = $cds->feat_end;
		my $phase  = $cds->phase;
		my $strand = $cds->strand;

		# search SNPs
		DEBUG "Fetching SNPs for ${id}, CDS C${chr}:${start}..${end}, ${strand}-strand, offset ${phase}";
		my $snps = $schema->resultset("Snp")->search({
			chromosome_id => $chr,
			position => { '>=' => $start, '<=' => $end },
		});

		# iterate over SNPs, merge contrasts
		my %merged;
		while( my $snp = $snps->next ) {

			# lookup variables to use in hash keys
			my $pos = $snp->position;
			my $ref = $snp->ref;
			my $alt = $snp->alt;

			# populate data structure
			$merged{$pos} = {} if not $merged{$pos};
			$merged{$pos}->{$ref} = {} if not $merged{$pos}->{$ref};
			$merged{$pos}->{$ref}->{$alt} = [] if not $merged{$pos}->{$ref}->{$alt};

			# store contrast
			push @{ $merged{$pos}->{$ref}->{$alt} }, $snp->contrast;

		}

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

		# navigate data structure
		POS: for my $pos ( sort { $a <=> $b } keys %merged ) {
			REF: for my $ref ( sort { $a cmp $b } keys %{ $merged{$pos} } ) {
				ALT: for my $alt ( sort { $a cmp $b } keys %{ $merged{$pos}->{$ref} } ) {				
					
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
						$snp_coord = $pos - $start;
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
					
					# compute whether synonymous and whether observed
					# reference allele matches the expectation
					my ( $is_nonsyn, $exp ) = is_nonsyn(
						'seq' => $seq,
						'ref' => $cref,
						'alt' => $calt,
						'pos' => $snp_coord,
						'str' => $strand,
					);
					
					# if there's no return values there's nothing to report
					next ALT unless $is_nonsyn;

					# prepare and print result
					my $contrast = join ',', @{ $merged{$pos}->{$ref}->{$alt} };
					my @result = ( $id, $chr, $start, $end, $phase, $strand, $pos, $is_nonsyn, $ref, $alt, $exp, $contrast );
					print join("\t", @result), "\n";
				}
			}
		}

	}
}

sub is_nonsyn {
	my %args = @_;
	my $raw = $args{seq}->seq;
	
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
	
	# the $retval is a switch where '=' means the observed allele extracted from
	# the reference genome matches the expected from GATK; '!' means there is
	# a mismatch.
	my $retval = '=';
	if ( $obs_ref ne $exp_ref ) {
		ERROR "Error: $obs_ref != $exp_ref (relative position: $pos)";
		$retval = '!';
	}
	
	# replace the reference allele with the alternative and test for synonymy
	substr( $raw, $pos, $length, $args{alt} );
	return $ctable->translate($raw) eq $ctable->translate($args{seq}->seq) ? ( 'syn', $retval ) : ( 'nonsyn', $retval );
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
