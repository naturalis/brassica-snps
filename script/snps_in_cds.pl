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
	while( my $cds = $cdss->next ) {

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

		# navigate data structure
		for my $pos ( sort { $a <=> $b } keys %merged ) {
			for my $ref ( sort { $a cmp $b } keys %{ $merged{$pos} } ) {
				for my $alt ( sort { $a cmp $b } keys %{ $merged{$pos}->{$ref} } ) {
				
					# get coding, in-frame, reference sequence
					my $seq = get_refseq(
						'chr'    => $chr,
						'start'  => $start,
						'stop'   => $end,
						'phase'  => $phase,
						'strand' => $strand,
					);
					
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
					
					# compute whether synonymous
					my $is_nonsyn = is_nonsyn(
						'seq' => $seq,
						'ref' => $ref,
						'alt' => $alt,
						'pos' => $snp_coord,						
					);

					# prepare and print result
					DEBUG Dumper(\%merged);
					my $contrast = join ',', @{ $merged{$pos}->{$ref}->{$alt} };
					my @result = ( $id, $chr, $start, $end, $phase, $strand, $pos, $is_nonsyn, $ref, $alt, $contrast );
					print join("\t", @result), "\n";
				}
			}
		}

	}
}

sub is_nonsyn {
	my %args = @_;
	my $raw = $args{seq}->seq;
	my $exp_ref = substr( $raw, $args{pos}, length($args{ref}) );
	if ( $exp_ref ne $args{ref} ) {
		WARN $args{pos}, "\t", $raw;
		WARN "Error: $exp_ref != " . $args{ref};
	}
	substr( $raw, $args{pos}, length($args{ref}), $args{alt} );
	return $ctable->translate($raw) eq $ctable->translate($args{seq}->seq) ? 'syn' : 'nonsyn';
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
