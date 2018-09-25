#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use My::Brassica;

=pod

Given a Brassica oleracea stable ID (gene/CDS/exon) and a specification of the
feature type, print the genomic coordinate(s) of the feature in a GATK-style
.intervals or .list format, see: https://software.broadinstitute.org/gatk/documentation/article.php?id=1319

=cut

# process command line arguments
my $db = '/home/ubuntu/data/reference/sqlite/snps.db';
my $id;
my $type = 'gene'; # can be 'gene', 'CDS', 'exon', etc. I.e. terms from the SO
my $bed;
GetOptions(
	'db=s'   => \$db,
	'id=s'   => \$id,
	'type=s' => \$type,
	'bed'    => \$bed,
);

# instantiate database
my $schema = My::Brassica->connect("dbi:SQLite:$db");

# run query
my $features = $schema->resultset('Feature')->search({ attributes => { LIKE => "%ID=${type}:${id}%" } });

# iterate over resultset
while( my $f = $features->next ) {
  if ( $bed ) {
    print 'C', $f->chromosome_id, "\t", $f->feat_start - 1, "\t", $f->feat_end - 1, "\n";
  }
  else {
	  print 'C', $f->chromosome_id, ':', $f->feat_start, '-', $f->feat_end, "\n";
  }
}
