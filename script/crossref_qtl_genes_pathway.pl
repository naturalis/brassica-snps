#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# e.g. results/genes.txt and results/targeted/oleracea_flowering_time_genes_frontiers.tsv
my ( $genelist, $table ) = @ARGV;

# read QTL gene
my %qtlgene;
{
    open my $in, '<', $genelist or die $!;
    while (<$in>) {
        chomp;
        if (/(\d+) (Bo.+)/) {
            my ($count, $gene) = ($1, $2);
            $qtlgene{$gene} = $count;
        }
    }
}

# read table
my @header;
open my $in, '<', $table or die $!;
while(<$in>) {
    chomp;
    my @line = split /\t/, $_;
    if (not @header) {
        @header = @line;
        push @header, 'N_CONTRASTS';
        print join("\t", @header), "\n";
    }
    else {
        my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
        my $gene = $record{'OLERACEA_ID'};
        if ( my $count = $qtlgene{$gene} ) {
            $record{'N_CONTRASTS'} = $count;
            my @values = @record{@header};
            print join("\t", @values), "\n";
        }
    }
}