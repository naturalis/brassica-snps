#!/usr/bin/perl
use strict;
use warnings;
use My::Brassica;
use Getopt::Long;

=pod

Given a gene ID, looks up its coding sequences and the QTL
regions within which it falls.

=cut

# process command line arguments
my $db = '/home/ubuntu/data/reference/sqlite/snps.db';
my $gene;
GetOptions(
	'db=s'   => \$db,
	'gene=s' => \$gene,
);
die "Usage: $0 -g <gene> [-d <db>]" if not $gene;

# connect to database
my $schema = My::Brassica->connect("dbi:SQLite:$db");

# select * from features where attributes like '%ID=CDS:Bo6g103730%';
my $cdss = $schema->resultset('Feature')->search({ attributes => { LIKE => "%ID=CDS:${gene}.1%" } });

# iterate over coding sequences
while( my $cds = $cdss->next ) {

	# fetch coordinates
	my $chr   = $cds->chromosome_id;
	my $start = $cds->feat_start;
	my $end   = $cds->feat_end;

	# search for QTL regions
	my $qtls = $schema->resultset('QtlRegion')->search({
		chromosome_id => $chr,
		start => { '<=' => $start },
		end   => { '>=' => $end }
	});

	# iterate over QTL regions
	while( my $qtl = $qtls->next ) {

		# fetch ID and contrast
		my $id = $qtl->qtl_region_id;
		my $bc = $qtl->bsa_contrast;

		# print output
		print join("\t", $gene, $chr, $start, $end, $id, $bc), "\n";
	}
}
