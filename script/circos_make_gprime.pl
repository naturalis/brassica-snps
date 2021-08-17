#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use My::Brassica;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':simple';

# documented in supplementary materials

# process command line arguments
my $dbfile = 'snps.db';
my $window = 1_000_000;
my $contrast = 'EF-IF';
my $verbosity = WARN;
GetOptions(
    'dbfile=s'   => \$dbfile,
    'window=i'   => \$window,
    'contrast=s' => \$contrast,
    'verbose+'   => \$verbosity,
);

# instantiate services
Bio::Phylo::Util::Logger->new( '-level' => $verbosity, '-class' => 'main' );
my $db = My::Brassica->connect("dbi:SQLite:$dbfile");

# iterate over chromosomes
my $crs = $db->resultset('Feature')->search({ feature_type => 'chromosome' });
while( my $chr = $crs->next ) {
    my $id = $chr->chromosome_id;
    my $length = $chr->feat_end;
    for ( my $i = 1; $i <= $length; $i += $window ) {
        my $clause =  {
            'chromosome_id' => $id,
            'contrast'      => $contrast,
            'position'      => {
                '>=' => $i,
                '<=' => $i + $window,
            }
        };
        my $sum = 0;
        my $count = 0;
        my $snps = $db->resultset('Snp')->search( $clause );
        while( my $snp = $snps->next ) {
            $sum += $snp->g_prime;
            $count++;
        }
        my $average = ( $count ? ($sum/$count) : 2.5 );
        #warn "$sum $count $average";
        printf("bol%i %i %i %f\n", $id, $i, ( ($i+$window)>$length ? $length : $i+$window ), $average );
    }
}

