#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::GenBank;
use Log::Log4perl qw(:easy);
use Bio::Tools::Run::StandAloneBlastPlus;
Log::Log4perl->easy_init($DEBUG);

# process command line arguments
my $interval;
my $taxon = 3702;
GetOptions(
	'interval=s' => \$interval,
	'taxon=i'    => \$taxon,
);

# location of reference genome
my $ref = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes';

# read intervals file, concatenate CDSs
my $fasta = ">query\n";
open my $fh, '<', $interval or die $!;
while(<$fh>) {
	chomp;
	if ( /^(C\d):(\d+)-(\d+)/ ) {
		my ( $chr, $start, $stop ) = ( $1, $2, $3 );
		my $cds = `fastacmd -d $ref -s $chr -L $start,$stop`;
		my @lines = split /\n/, $cds;
		shift @lines;
		$fasta .= join '', @lines;
	}
}
DEBUG $fasta;

# create seq object
my $seq = Bio::SeqIO->new(
	'-string' => $fasta,
	'-format' => 'fasta',
)->next_seq;

# instantiate online blast factory 
my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
	'-db_name' => 'nr',
 	'-remote'  => 1
);

# instantiate genbank lookup client
my $gb = Bio::DB::GenBank->new;

# perform the blast query
my $result = $fac->tblastn( '-query' => $seq );

# iterate over results
my @filtered;
while( my $hit = $result->next_hit ) {
	my $acc = $hit->accession;
	DEBUG $acc;

	# versioned refseq
	if ( $acc =~ /^N(?:M|R|P)_.+?\.\d/ ) {
		my $seq = $gb->get_Seq_by_version($acc);
		if ( $seq->species->ncbi_taxid == $taxon ) {
			print $acc, "\n";
		}
	}
}

