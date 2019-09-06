#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use My::GeneDB;
use Getopt::Long;

# process command line arguments
my $infile;
my $rootfolder;
GetOptions(
    'infile=s'     => \$infile,
    'rootfolder=s' => \$rootfolder,
);

# create db
my $db = My::GeneDB->new($infile);

# iterate over genes
for my $gene ( $db->records ) {
    my $id = $gene->gene_id;

    # pick the file to parse. Might be reverse complemented, if it exists
    my $seqfile = "${rootfolder}/${id}/combined-aligned.fasta";
#    if ( -e "${rootfolder}/${id}/combined.fasta.revcom" ) {
#        $seqfile = "${rootfolder}/${id}/combined.fasta.revcom";
#    }
#    else {
#        $seqfile = "${rootfolder}/${id}/combined.fasta";
#    }

    # instantiate the FASTA reader
    my $seqio = Bio::SeqIO->new(
        '-format' => 'fasta',
        '-file'   => $seqfile,
    );

    # start reading
    while( my $seq = $seqio->next_seq ) {
        my $aaseq = $seq->translate->seq;
        if ( $aaseq =~ /\*.+/ ) {
            warn $seqfile;
        }
    }
}