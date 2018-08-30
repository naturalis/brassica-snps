#!/usr/bin/perl
use strict;
use warnings;

# define globals
my @files = @ARGV;
my $region_id = 1;

# iterate over files
for my $file (@files) {

	# prepare constraint identifier
	my $constraint = $file;
	$constraint =~ s/^.+?([A-Z][A-Z]-[A-Z][A-Z])\.csv$/$1/;

	# iterate over records in file
	open my $fh, '<', $file or die $!;
	while(<$fh>) {

		# skip header line
		next if /"CHROM"/;

		# pre-process record line
		chomp;
		my @fields = split /,/, $_;
		$fields[0] =~ s/"C(\d+)"/$1/; # turn chromosome name into FK

		# prefix primary key, postfix constraint identifier
		unshift @fields, $region_id++;
		push @fields, $constraint;

		# print output
		print join("\t", @fields), "\n";
	}
}
