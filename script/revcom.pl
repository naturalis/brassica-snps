#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments
my $refseq;
my $file;
my $outfile;
GetOptions(
    'refseq=s'  => \$refseq,
    'file=s'    => \$file,
    'outfile=s' => \$outfile,
);

my ( @deflines, %data );
open my $fh, '<', $file or die $!;
while(<$fh>) {
    chomp;
    if ( />/ ) {
        my $def = $_;
        push @deflines, $def;
        $data{$def} = [];
    }
    else {
        my $def = $deflines[-1];
        push @{ $data{$def} }, $_;
    }
}

open my $outfh, '>', $outfile or die $!;
for my $def ( @deflines ) {
    my @seq;
    if ( $def =~ />\Q$refseq\E/ ) {
        for my $line ( @{ $data{$def} } ) {
            $line =~ tr/aAcCgGtT/tTgGcCaA/;
            $line = reverse $line;
            push @seq, $line;
        }
    }
    else {
        @seq = @{ $data{$def} };
    }
    @seq = reverse @seq;
    print $outfh join( "\n", $def, @seq ), "\n";
}

