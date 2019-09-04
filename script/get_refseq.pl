#!/usr/bin/perl
use v5.10;
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

# process command line arguments
my $intervals;
my $seqfile;
GetOptions(
    'intervals=s' => \$intervals,
    'seqfile=s'   => \$seqfile,
);

# read intervals
my %interval;
open my $in, '<', $intervals or die $!;
while(<$in>) {
    chomp;
    if ( /^(C\d+):(\d+?)-(\d+)$/ ) {
        my ( $chromo, $start, $stop ) = ( $1, $2, $3 );
        $interval{$chromo} = [] if not $interval{$chromo};
        push @{ $interval{$chromo} }, [ $start, $stop ];
    }
}

# sort intervals
for my $chromo ( keys %interval ) {
    my @sorted = sort { $a->[1] <=> $b->[1] } @{ $interval{$chromo} };
    $interval{$chromo} = \@sorted;
}

# setup reader args, optionally decompress GZ on the fly
my %args = ( '-format' => 'fasta' );
if ( $seqfile =~ /\.gz$/ ) {
    open my $zcat, "gunzip -c $seqfile |" or die $!;
    $args{'-fh'} = $zcat;
}
else {
    $args{'-file'} = $seqfile;
}

# create fasta line
my $line = $seqfile;
$line =~ s/.+\//>/;
$line =~ s/\.gz$//;
$line =~ s/\.fa$//;

# start reading
my $io = Bio::SeqIO->new(%args);
while( my $seq = $io->next_seq ) {
    my $chromo = $seq->display_id;
    if ( $interval{$chromo} ) {
        my @segments = map { $seq->subseq(@$_) } @{ $interval{$chromo} };
        say join "\n", $line, @segments;
    }
}