#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use My::Brassica;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use File::Temp qw(tempfile);
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);

# process command line arguments
my $infile = '/home/ubuntu/brassica-snps/results/targeted/oleracea_flowering_time_genes_frontiers.tsv';
my $bam    = '/home/ubuntu/data/gDNA/group-4/group-4_pe.sorted.bam';
my $db     = '/home/ubuntu/data/reference/sqlite/snps.db'; 
my $ref    = '/home/ubuntu/data/reference/Brassica_oleracea_chromosomes.fa';
my $jar    = '/usr/local/src/GenomeAnalysisTK.jar';
GetOptions(
  'infile=s' => \$infile,
  'db=s'     => \$db,
  'ref=s'    => \$ref,
  'bam=s'    => \$bam,
  'jar=s'    => \$jar,
);

# connect to database, initialized codon table
my $schema = My::Brassica->connect("dbi:SQLite:$db");
my $ctable = Bio::Tools::CodonTable->new();

# read infile
my %genes; # gene_name => [ Bo1, Bo2, ... BoN ]
my @header;
open my $fh, '<', $infile or die $!;
while(<$fh>) {
  chomp;
  my @line = split /\t/, $_;
  
  # read header
  if ( not @header ) {
    @header = @line;
  }
  
  # process record
  else {
    
    # it's pointless to try to investigate the copies on scaffolds:
    # we didn't map the reads to them, didn't do the QTLSeqR analysis, etc.
    if ( $line[0] =~ /^Bo\dg\d+/ ) {
      my %record    = map { $header[$_] => $line[$_] } 0 .. $#header;
      my $gene_name = $record{'AT_ID'};
      my $bo_id     = $record{'OLERACEA_ID'};
      $genes{$gene_name} = [] if not $genes{$gene_name};
      push @{ $genes{$gene_name} }, $bo_id;
    }
  }
}

# attributes of gene locations detected in PastedGraphic-4.pdf
# - 'deletion'    - will detect this by computing coverage (should be very low)
# - 'duplication' - will detect this by computing coverage (should be high)
# - 'no variant'  - if there are no SNPs?
# - 'silent SNP'  - probably just a SNP outside of a CDS
# - 'synonymous SNP'     - splice alternative allele in reference sequence, compute AA translation, compare
# - 'non-synonymous SNP' - splice alternative allele in reference sequence, compute AA translation, compare
# - 'stop codon variant' - won't do, not seen in interesting genes anyway
# - 'splice variant'     - won't do, not seen in interesting genes anyway
# - 'InDel'              - i.e. length(ref allele) != length(alt allele)
print join( "\t", qw(gene_name copy_id coverage silent synonymous non_synonymous stop_codon indel_coding indel_non_coding) ), "\n";
for my $gene_name ( sort { $a cmp $b } keys %genes ) {
  INFO "gene: ${gene_name}";
  for my $bo_id ( sort { $a cmp $b } @{ $genes{$gene_name} } ) {
    INFO "\tcopy: ${bo_id}";

    # get gene and CDS coordinates, get SNPs
    my $cdss = $schema->resultset('Feature')->search({ 'attributes' => { 'LIKE' => "%ID=CDS:$bo_id.%" } });
    my $gene = $schema->resultset('Feature')->single({ 'attributes' => { 'LIKE' => "%ID=gene:$bo_id%" } });
    my $snps = $schema->resultset("Snp")->search({
		  'chromosome_id' => $gene->chromosome_id,
			'position'      => {
			  '>=' => $gene->feat_start, 
			  '<=' => $gene->feat_end 
			 },
		});
		
		# calculate average coverage
		my $coverage = get_coverage(
		  'chr'   => $gene->chromosome_id,
		  'start' => $gene->feat_start,
		  'stop'  => $gene->feat_end,
		);
		  
		# make lookup of coordinates
		my ( %lookup, @cds );
		while( my $cds = $cdss->next ) {
		  my $fs = $cds->feat_start;
		  $lookup{$fs} = $cds->feat_end;
		  push @cds, $cds;
		}
		  
		# iterate over SNPs
		my ($silent,$syn,$nonsyn,$indel_coding,$indel_noncoding,$stop_codon,%seen) = (0,0,0,0,0,0);
		SNP: while( my $snp = $snps->next ) {
		  
		  # lookup column values
		  my $pos = $snp->position;
		  my $ref = $snp->ref;
		  my $alt = $snp->alt;
		  
		  # populate seen hash
		  $seen{$pos} = {} if not $seen{$pos};
		  $seen{$pos}->{$ref} = {} if not $seen{$pos}->{$ref};
		  
		  # have not seen this ref/alt combination yet
		  if ( not $seen{$pos}->{$ref}->{$alt}++ ) {
  		  my @starts = grep { $pos >= $_ } keys %lookup;
  		  @starts = grep { $pos <= $lookup{$_} } @starts;
  		  
  		  # is a coding snp
  		  if ( scalar(@starts) == 1 ) {
  		    
  		    # the SNP is an indel if the length of the ref allele is not the same as the length of the alt allele
  		    $indel_coding += ( length($ref) != length($alt) ? 1 : 0 );
  		    my ($cds) = grep { $pos >= $_->feat_start && $pos <= $_->feat_end } @cds;
  		    
      		# get coding, in-frame, reference sequence
      		my $seq;
      		eval {
      			$seq = get_refseq(
      				'chr'    => $cds->chromosome_id,
      				'start'  => $cds->feat_start,
      				'stop'   => $cds->feat_end,
      				'phase'  => $cds->phase,
      				'strand' => $cds->strand,
      			);
      		};
      		if ( $@ ) {
      		  my $template = 'Problem extracting %s[%s]:%d..%d (phase: %d)';
      		  my $msg = sprintf($template, $cds->chromosome_id, $cds->strand, $cds->feat_start, $cds->feat_end, $cds->phase);
      			ERROR $msg;
      			ERROR $@;
      			next SNP;
      		}  	
      		
					# adjust SNP coordinate
					my $snp_coord;
					if ( $cds->strand eq '-' ) {
						# CDS is on '-' strand, count 
						# backward from 3' location
						$snp_coord = $cds->feat_end - $pos;
					}
					else {
						# CDS is on '+' strand, count
						# forward relative to CDS start
						$snp_coord = $pos - $cds->feat_start;
					}
					$snp_coord -= $cds->phase; # is 0, 1 or 2  
					
					# revcom $ref & $alt if on '-' strand
					my ( $cref, $calt ) = ( $ref, $alt );
					if ( $cds->strand eq '-' ) {
						$cref =~ tr/ACGT/TGCA/;
						$calt =~ tr/ACGT/TGCA/;
						$cref = reverse($cref);
						$calt = reverse($calt);
					}	
					
					# compute whether synonymous and where, relative to the
					# local seq, the SNP is (i.e. in stop codon?)
					my ( $is_nonsyn, $pos ) = is_nonsyn(
						'seq' => $seq,
						'ref' => $cref,
						'alt' => $calt,
						'pos' => $snp_coord,
						'str' => $cds->strand,
					);	
					$is_nonsyn ? $nonsyn++ : $syn++;
					$stop_codon++ if $pos >= (length($seq->seq) - 2);
  		  }
  		  elsif ( scalar(@starts) > 1 ) {
  		    ERROR "This should be impossible";
  		  }
  		  else {
  		    
  		    # the SNP is an indel if the length of the ref allele is not the same as the length of the alt allele
  		    $indel_noncoding += ( length($ref) != length($alt) ? 1 : 0 );
  		    $silent++;
  		  }
		  }
		}
		
		# print output
		print join( "\t", $gene_name,$bo_id,$coverage,$silent,$syn,$nonsyn,$stop_codon,$indel_coding,$indel_noncoding ), "\n";
  }
}

sub is_nonsyn {
	my %args = @_;
	my $raw = $args{seq}->seq;
	
	# Compute starting position differently depending on forward or reverse 
	# strand: we need to go upstream for the the reverse strand.
	my $pos = $args{pos};
	if ( $args{str} eq '-' ) {
	
		# on the '-' strand the start coordinate for alleles needs to
		# be adjusted, but using 0-based indexing. I.e. this has no
		# effect for alleles of 1bp length, which is nearly all of them.
		my $l = length($args{ref}) - 1;
		$pos -= $l;
	}
	
	# allele cannot be nonsynonymous if we've moved upstream outside of the phased CDS
	return if $pos < 0;
	
	# Compute length so that we don't go outside of the phased CDS
	my $length = length($args{ref});
	if ( ( $length + $pos ) > length($raw) ) {
		$length = length($raw) - $pos;
	}
	
	# extract observed allele
	my $obs_ref = substr( $raw, $pos, $length ); 

	# allele cannot be nonsynonymous if we've moved downstream outside of the phased CDS
	return unless $obs_ref;
	
	# the $retval is a switch where '=' means the observed allele extracted from
	# the reference genome matches the expected from GATK; '!' means there is
	# a mismatch.
	my $retval = $pos;
	
	# replace the reference allele with the alternative and test for synonymy
	substr( $raw, $pos, $length, $args{alt} );
	return $ctable->translate($raw) eq $ctable->translate($args{seq}->seq) ? ( 'syn', $retval ) : ( 'nonsyn', $retval );
}

sub get_refseq {
	my %args  = @_;
	my ( $chr, $start, $stop ) = @args{ qw(chr start stop) };
	$chr = "C$chr" if $chr !~ /^C/; # translate chromosome foreign key to FASTA ID
	my $command = "fastacmd -d $ref -s $chr -L $start,$stop";
	DEBUG $command;
	my $fasta = `$command`;
	DEBUG $fasta;
	
	# parse raw FASTA
	my $seq = Bio::SeqIO->new( 
		'-string' => $fasta, 
		'-format' => 'fasta',
	)->next_seq;
	
	# reverse complement if CDS is on '-' strand
	$seq = $seq->revcom if $args{'strand'} eq '-';
	
	# truncate sequence if there is a phase offset
	$seq = $seq->trunc( $args{'phase'} + 1, $args{'stop'} - $args{'start'} ) if $args{'phase'};
	return $seq;
}

sub get_coverage {
  my %args = @_;
  
  # make bed file
  my ($fh,$filename) = tempfile( 'tmpXXXX', SUFFIX => '.intervals' );
  print $fh 'C' . $args{'chr'} . ':' . $args{'start'} . '-' . $args{'stop'};
  close $fh;
  
  # run the command
  open my $stdout, "java -jar $jar -T DepthOfCoverage -R $ref -I $bam -L $filename -l fatal |" or die $!;
  my ( $seen_header, $average_coverage );
  LINE: while(<$stdout>) {
    if ( /^Target/ ) {
      $seen_header++;
      next LINE;
    }
    if ( $seen_header ) {
      my @line = split /\t/, $_;
      $average_coverage = $line[2];
      last LINE;
    }
  }
  
  # delete bed file
  unlink $filename;
  return $average_coverage;
}