#!/usr/bin/perl
use strict;
use warnings;
use My::GeneDB;
use Getopt::Long;
use Bio::DB::GenBank;
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $verbosity = WARN;
my $infile;
GetOptions(
    'infile=s' => \$infile,
    'verbose+' => \$verbosity,
);

# configure services
Bio::Phylo::Util::Logger->new( '-level' => $verbosity );
my $db = My::GeneDB->new($infile);
my $gb = Bio::DB::GenBank->new;

for my $r ( $db->records ) {

    # fetch the ref seq identifier from the record
    my $acc = $r->link_name('refseq');
    INFO "Going to fetch sequence $acc from GenBank";

    # do the lookup
    my $seq = $gb->get_Seq_by_acc($acc);
    DEBUG "Returned by GenBank: $seq";

    # iterate over features
    for my $f ( $seq->get_SeqFeatures ) {

        # there is a top-level gene annotation
        if ( $f->primary_tag eq 'gene' ) {

            # top-level gene does not have to have a sub-tag
            if ( $f->has_tag('gene') ) {
                my ($name) = $f->get_tag_values('gene');
                $r->gene_name($name);
                INFO "Found gene name '$name'";
            }
            else {
                DEBUG "Top-level gene tag had no gene sub tag"
            }
        }
    }
}

print $db->to_markdown;