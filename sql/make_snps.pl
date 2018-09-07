#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;

my $id = 1; 
my @files = @ARGV; 

# use CSV reader because some fields are quoted, comma-separated triples themselves
my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: ".Text::CSV->error_diag ();	
for my $file (@files) {

	# switch to jump over the first (header) row
	my $header;
	open my $fh, "<:encoding(utf8)", $file or die $!;
	while( my $r = $csv->getline($fh) ) {
		next unless $header++;
		my @fields = @{ $r };

		# remove all quotes: will emit tab-separated anyway
		$_ =~ s/"//g for @fields;

		# input files each have IDs starting at 1, hence use counter across files
		$fields[0] = $id++;

		# turn column 2 into foreign key
		$fields[1] =~ s/C//;

		print join("\t", @fields), "\n";
	}
}
