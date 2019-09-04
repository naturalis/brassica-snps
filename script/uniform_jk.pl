#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use My::GeneDB;
use Getopt::Long;

# process command line arguments
my $infile;
my $rootfolder;
my $modelseq = 'Brassica_oleracea.v2.1.dna.toplevel';
GetOptions(
    'infile=s'     => \$infile,
    'rootfolder=s' => \$rootfolder,
    'modelseq=s'   => \$modelseq,
);

# configure services
my $db = My::GeneDB->new($infile);

# iterate over genes
GENE: for my $gene ( $db->records ) {
    my $id = $gene->gene_id;

    # pick the file to parse. Might be reverse complemented, if it exists
    my $seqfile = "${rootfolder}/${id}/combined-aligned.fasta";

    # instantiate the FASTA reader
    my $seqio = Bio::SeqIO->new(
        '-format' => 'fasta',
        '-file'   => $seqfile,
    );

    # start reading
    my ( %matrix, $length );
    while( my $seq = $seqio->next_seq ) {
        my $aaseq = $seq->translate->seq;
        my $id = $seq->id;
        $matrix{$id} = [ split //, $aaseq ];
        $length = length($aaseq);
    }

    # start checking
    my $mismatches = 0;
    for my $i ( 0 .. $length - 1 ) {
        my %seen;
        SEQ: for my $id ( keys %matrix ) {
            next SEQ if $id eq $modelseq;
            my $aa = $matrix{$id}->[$i];
            $seen{$aa}++;
        }

        # check if model AA is different from seen
        my $maa = $matrix{$modelseq}->[$i];
        $mismatches++ if not $seen{$maa};
    }
    $gene->uniform_jk_snp($mismatches);
}

print $db->to_markdown;
