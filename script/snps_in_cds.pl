#!/usr/bin/perl
use strict; 
use warnings; 
use Bio::SeqIO;
use Getopt::Long; 
use My::Brassica;

# process command line arguments
my @genes;
while(<>) {
	chomp;
	push @genes, $_ if /\S/;
}

# connect to database
my $db     = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $ref    = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes';
my $schema = My::Brassica->connect("dbi:SQLite:$db");

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

					# prepare and print result
					my $contrast = join ',', @{ $merged{$pos}->{$ref}->{$alt} };
					my @result = ( $id, $chr, $start, $end, $phase, $strand, $pos, $ref, $alt, $contrast );
					print join("\t", @result), "\n";
				}
			}
		}

	}
}

sub get_seq {
	my %args  = @_;
	my ( $chr, $start, $stop ) = @args{ qw(chr start stop) };
	my $fasta = `fastacmd -d $ref -s $chr -L $start,$stop`;
	my $seq = Bio::SeqIO->new( -string => $fasta, -format => 'fasta' )->next_seq;
	$seq = $seq->revcom if $args{'strand'} eq '-';
	$seq = $seq->trunc( $args{'phase'} + 1, $args{'stop'} - $args{'start'} ) if $args{'phase'};
	return $seq;
}
