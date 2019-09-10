#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV 'csv';
use Getopt::Long;
use Data::Dumper;
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $infile; # results/refinement/gProfiler_boleracea_9-9-2019 21-58-35__intersections.csv
my $verbosity = WARN;
GetOptions(
    'infile=s' => \$infile,
    'verbose+' => \$verbosity,
);

# read the file
Bio::Phylo::Util::Logger->new( '-level' => $verbosity, '-class' => 'main' );
my $aoh = csv( 'file' => $infile, 'headers' => 'auto' );

# collect the gene ids, update the intersections as we go along
my ( %gene, @header );
for my $term ( @$aoh ) {
    if ( my $i = $term->{intersections} ) {

        # parse the string with genes
        my @genes = split ',', $i;
        INFO "Parsed ${\scalar(@genes)} from record";

        @header = keys %$term;         # should stay the same for all records
        $term->{$_}++ for @genes;      # will be a boolean 0/1 for intersection or no
        $gene{$_}++ for @genes;        # will be total intersections, maybe conditional formatting
        delete $term->{intersections}; # clean up
    }
    else {
        ERROR "Missing record field 'intersections': " . Dumper($term);
    }
}

# start printing the output. genes are now sorted by chromosomal position
my @genes = sort { $a cmp $b } keys %gene;
print join("\t", @header, @genes ), "\n";

# print the number of intersections for the focal gene
print join("\t", map { $gene{$_} || '' } @header, @genes ), "\n";
for my $term ( @$aoh ) {

    # here we complete the boolean, i.e. '0' if no recorded intersection
    print join("\t", map { $term->{$_} || 0 } @header, @genes ), "\n";
}