use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Tools::Run::StandAloneBlastPlus;

# process command line arguments
my $dbname;
my $tmpfile;
my $linkages;
GetOptions(
	'dbname=s'   => \$dbname,
	'tmpfile=s'  => \$tmpfile,
	'linkages=s' => \$linkages,
);

# prepare globals
my $pk = 1;
my $bp = Bio::Tools::Run::StandAloneBlastPlus->new( -db_name => $dbname );
open my $in, '<', $linkages or die $!;

while(<$in>) {
	# pre-process input
	next if /^linkage_group/;
	chomp;
	my @fields = split /\t/;

	$fields[0] =~ s/^C//;	# update fields

	
	# do search
	for my $i ( -2, -1 ) {

		# write tmp file
		my $seq = $fields[$i];
		open my $fh, '>', $tmpfile or die $!;
		print $fh ">query\n$seq\n";
		close $fh;

		# do query
		my $result = $bp->blastn( -query => $tmpfile );
		print $result->num_hits, "\n";
#		while( my $hit = $result->next_hit ) {
#			print Dumper($hit);
#		}
	}


	# update keys
	unshift @fields, $pk++;
	push @fields, ( undef, undef );
#	print join( "\t", @fields ), "\n";
}
