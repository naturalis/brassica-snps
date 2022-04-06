#!/usr/bin/perl
use v5.14;
use strict;
use warnings;
use Getopt::Long;
use Convert::Color;
use Text::CSV 'csv';
use Text::Wrap 'wrap';
use List::Util qw'min max';
use GO::OntologyProvider::OboParser;
use Bio::Phylo::Util::Logger ':simple';

# documented in supplementary materials

# Usage:
# go_treebuilder.pl --obo=go-basic.obo | dot -Tsvg -o go_subgraph.svg

# process command line arguments
my $obo       = 'go-basic.obo';
my $base      = 'GO:0003006';
my $root      = '../results';
my @contrasts = qw(EF-IF EF-LF EF-NF IF-LF IF-NF LF-NF);
my $verbosity = WARN;
GetOptions(
    'obo=s'       => \$obo,
    'base=s'      => \$base,
    'root=s'      => \$root,
    'verbose+'    => \$verbosity,
    'contrasts=s' => \@contrasts,
);

# configure objects and services
Bio::Phylo::Util::Logger->new( '-level' => $verbosity, '-class' => 'main' );
my $go = GO::OntologyProvider::OboParser->new( 'ontologyFile' => $obo, 'aspect' => 'P' );
local $Text::Wrap::columns = 30;
local $Text::Wrap::separator='<BR/>';

# initialize the data structure that maps between safe IDs (GO_\d+), term objects, and gene counts per contrast
my $safe_base = $base =~ s/:/_/r;
my %terms = ( $safe_base => { obj => $go->nodeFromId($base) } );
my @counts; # for sorting by value, to compute palette

# scan the files, which should be ${root}/${contrast}/enriched_${base}.tsv, where $base =~ s/:/_/;
for my $contrast ( @contrasts ) {
    my $file = "${root}/${contrast}/enriched_${safe_base}.tsv";
    if ( -e $file ) {
        my $tsv = csv(
            'in'       => $file,
            'headers'  => 'auto',
            'sep_char' => "\t",
        );
        for my $t ( @$tsv ) {
            my $term  = $t->{GO_acc} =~ s/:/_/r;
            # my $count = $t->{queryitem};
            $terms{$term} = { obj => $go->nodeFromId( $t->{GO_acc} ) } if not $terms{$term};
            $terms{$term}->{$contrast} = {
                'count' => $t->{queryitem},
                'pval'  => log($t->{pvalue}),
            };
            push @counts, log($t->{pvalue});
        }
    }
    else {
        WARN "$file not found";
    }
}

# build the initial graph topology
my @nodes = keys %terms;
my %ancestors;
for my $i ( 0 .. $#nodes ) {
    for my $j ( 0 .. $#nodes ) {
        next if $i == $j;
        my ( $id1, $id2 ) = @nodes[ $i, $j ];
        if ( $terms{$id1}->{obj}->isADescendantOf($terms{$id2}->{obj}) ) {
            $ancestors{$id1} = {} if not $ancestors{$id1};
            $ancestors{$id1}->{$id2}++;
        }
    }
}

# prune graph to remove transitive ancestors, keeping only direct ones
for my $node ( @nodes ) {
    my @anc = keys %{ $ancestors{$node} };
    if ( @anc > 1 ) {
        for my $i ( 0 .. $#anc ) {
            for my $j ( 0 .. $#anc ) {
                next if $i == $j;
                my ( $id1, $id2 ) = @anc[ $i, $j ];
                if ( $terms{$id1}->{obj}->isADescendantOf($terms{$id2}->{obj}) ) {
                    delete $ancestors{$node}->{$anc[$j]};
                }
            }
        }
    }
}

# print dot language header
print "graph go_sub {\n";
print "\tnode [fontname = \"arial\", width = 1.8, height = 1.6, fontsize = 12, labelloc = \"t\"];\n";

# print the node statements with wrapped term definitions and go ID
for my $node ( @nodes ) {

    # wrapped term, GO ID in gray
    my $label  = wrap('', '', $terms{$node}->{obj}->term );
       $label .= '<BR/><FONT COLOR="gray">' . $terms{$node}->{obj}->goid . '</FONT>';

    # build HTML-ish table
    my @cells;
    for my $c ( @contrasts ) {
        my $count = $terms{$node}->{$c}->{count} || 0;
        my $pval  = $terms{$node}->{$c}->{pval}  || 1;
        my $hsv;
        if ( $count ) {
            $hsv = make_color($pval, @counts);
        }
        else {
            $hsv = '0 0 0.7'; # white
        }
        push @cells, sprintf( '<TD WIDTH="25" HEIGHT="25" FIXEDSIZE="TRUE" BGCOLOR="%s">%01d</TD>', $hsv, $count );
    }
    my $table = sprintf(
        '<TABLE BORDER="0"><TR><TD COLSPAN="6">&nbsp;</TD></TR><TR>%s</TR><TR><TD COLSPAN="%i">%s</TD></TR></TABLE>',
        join( '', @cells ),
        scalar(@contrasts),
        $label
    );


    printf "%s [shape=box, label=<%s>];\n", $node, $table;
}

# print the topology
for my $node ( @nodes ) {
    for my $anc ( keys %{ $ancestors{$node} } ) {
        if ( $anc eq $safe_base ) {
            print "\t$anc -- $node [style=dotted];\n";
        }
        else {
            print "\t$anc -- $node;\n";
        }
    }
}

# print dot language footer
print "}\n";

# creates hexadecimal RGB code
sub make_color {
    my ( $value, @set ) = @_;
    my $min = min @set;
    my $max = max @set;
    my $ratio = ( $value - $min ) / ( $max - $min ); # i.e. normalized b/w 0 & 1

    # let's say that we want the highest value to be 'red' (0° HSV), i.e. hot, and the lowest 'yellow' (60°)
    my $H = ( $ratio ) * 90;
    return '#' . Convert::Color->new( "hsv:$H,1,1" )->as_rgb8->hex;
}
