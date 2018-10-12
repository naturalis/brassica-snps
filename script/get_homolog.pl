#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::GenBank;
use Log::Log4perl qw(:easy);
use Bio::Tools::Run::RemoteBlast;
Log::Log4perl->easy_init($DEBUG);

# process command line arguments
my $interval;
my $taxon = 3702;
my $ref = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes';
# or: /Users/rutger.vos/Dropbox/documents/projects/dropbox-projects/brassica/Brassica_oleracea_chromosomes
GetOptions(
	'interval=s' => \$interval,
	'taxon=i'    => \$taxon,
	'ref=s'      => \$ref,
);

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
my $fac = Bio::Tools::Run::RemoteBlast->new(
	'-data' => 'nr',
 	'-prog' => 'blastn'
);

# instantiate genbank lookup client
my $gb = Bio::DB::GenBank->new;

# perform the blast query
my @hits;
$fac->submit_blast( $seq );
while ( my @rids = $fac->each_rid ) {
  for my $rid ( @rids ) {
    my $rc = $fac->retrieve_blast($rid);
    
    # result is not an object, still waiting
    if ( !ref($rc) ) {
      $fac->remove_rid($rid) if $rc < 0;
      DEBUG "waiting...";
      sleep 5;
    }
    else {
      my $result = $rc->next_result();
      $fac->remove_rid($rid);
      while ( my $hit = $result->next_hit ) {
        push @hits, $hit;
      }
    }
  }
}

# iterate over results, filter
my @filtered;
for my $hit ( @hits ) {
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

