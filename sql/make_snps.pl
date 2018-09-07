#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;

my $id = 1; 
my @files = @ARGV; 
my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: ".Text::CSV->error_diag ();	
for my $file (@files) {
	my $header;
	open my $fh, "<:encoding(utf8)", $file or die $!;
	while( my $r = $csv->getline($fh) ) {
		next unless $header++;
		my @fields = @{ $r };
		$fields[0] = $id++;
		print join("\t", @fields), "\n";
	}
}
