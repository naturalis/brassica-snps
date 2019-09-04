#!/usr/bin/perl
use strict;
use warnings;

my @files = @ARGV;
my ( %data, @global_header );
for my $file ( @files ) {
    my @header;
    open my $fh, '<', $file or die $!;
    while(<$fh>) {
        chomp;
        my @line = split /\t/, $_;
        if ( not @header ) {
            @header = @line;
        }
        else {
            my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
            my $id = $line[0];
            if ( not $data{$id} ) {
                $record{'Dataset'} = [ $file ];
                $data{$id} = \%record;
            }
            else {
                warn;
                push @{ $data{$id} }, $file;
            }
        }
    }
    @global_header = @header;
}

push @global_header, 'Dataset';
print join("\t", @global_header), "\n";
for my $id ( sort { $a cmp $b } keys %data ) {
    my %record = %{ $data{$id} };
    my @values = map { ref $_ ? join(', ', @$_) : $_ } @record{@global_header};
    print join("\t", @values), "\n";
}