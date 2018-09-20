#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use GO::OntologyProvider::OboParser;

=pod

This script filters the records in the $infile to retain only those lines where the
listed GO term's enrichment: 

- has a significance level of <= $pval (default 0.05)
- belongs to the defined GO $section (default biological_process, P)
- is subtended by the provided $term (default GO:0003006 reproductive developmental process)

In addition, the script needs to full GO, in OBO format.

The input file is produced by the enrichment test at: http://bioinfo.cau.edu.cn/agriGO

=cut

# process command line arguments
my $infile;
my $term = 'GO:0003006'; # reproductive developmental process
my $ontology = 'go-basic.obo';
my $pval = 0.05;
my $fdrt = 0.05; # false discovery rate threshold
my $section = 'P';
GetOptions(
	'infile=s'   => \$infile,
	'term=s'     => \$term,
	'ontology=s' => \$ontology,
	'pval=f'     => \$pval,
	'fdrt=f'     => \$fdrt,
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
	if ( $record{term_type} eq $section and $record{pvalue} <= $pval and $record{FDR} <= $fdrt ) {
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