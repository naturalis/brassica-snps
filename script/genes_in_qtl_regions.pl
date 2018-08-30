#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use My::Brassica;

# process command line arguments
my $db = '/home/ubuntu/data/reference/sqlite/snps.db';
my $contrast = 'EF-LF';
GetOptions(
	'db=s'       => \$db,
	'contrast=s' => \$contrast,
);

# connect to database
my $schema = My::Brassica->connect("dbi:SQLite:$db");

# query regions
my $regions = $schema->resultset("QtlRegion")->search({ bsa_contrast => $contrast });
while( my $r = $regions->next ) {

	# get coordinates
	my $chr   = $r->chromosome->id;
	my $start = $r->start;
	my $end   = $r->end;

	# query features
	my $features = $schema->resultset("Feature")->search({ 
		'chromosome_id' => $chr,
		'feat_start'    => { '>' => $start },
		'feat_end'      => { '<' => $end },
		'feature_type'  => "gene"
	});
	while( my $f = $features->next ) {
		my $att = $f->attributes;
		if ( $att =~ m/ID=gene:([^;]+)/ ) {
			my $id = $1;
			print join("\t", $chr, $start, $end, $id), "\n";
		}
	}		
}
