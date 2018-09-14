#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use GO::OntologyProvider::OboParser;

# process command line arguments
my $infile;
my $term = 'GO:0003006'; # reproductive developmental process
my $ontology = 'go-basic.obo';
my $pval = 0.05;
my $section = 'P';
GetOptions(
	'infile=s'   => \$infile,
	'term=s'     => \$term,
	'ontology=s' => \$ontology,
	'pval=f'     => \$pval,
	'section=s'  => \$section,
);

# parse ontology
my $go = GO::OntologyProvider::OboParser->new(
	'ontologyFile' => $ontology,
	'aspect'       => 'P'
);

# fetch ancestor
my $anc = $go->nodeFromId($term);


# iterate over file
open my $fh, '<', $infile or die $!;
my @header;
LINE: while(<$fh>) {
	my @line = split /\t/;
	
	# store header line
	if ( not @header ) {
		@header = @line;
		print $_;
		next LINE;
	}
	
	# process record
	my %record = map { $header[$_] => $line[$_] } 0 .. $#line;
	
	# apply section and pval filter
	if ( $record{term_type} eq $section and $record{pvalue} <= $pval ) {
		if ( my $node = $go->nodeFromId($record{GO_acc}) ) {
		
			# print line if $node descends from provided node
			if ( $node->isADescendantOf($anc) ) {
				print $_;
			}
		}
		else {
			warn "NO NODE: $_";
		}
	}
}