#!/usr/bin/perl
# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use lib '/usr/local/src/biomart-perl/lib';
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

# documented in supplementary materials

=pod

Given lines of gene IDs (from the EnsEMBL Brassica oleracea genome annotation) provided on 
STDIN, do a lookup in BioMart. The result is a tab separated table with the following
columns:

- the input gene ID
- the UniProtKB identifier that corresponds to the gene product
- a GO term ID with which the UniProtKB record is annotated
- the GO term name
- the GO term definition
- the upper level category in the GO (e.g. biological process)

=cut

my $confFile = "/usr/local/src/biomart-perl/conf/martURLLocation.xml";

# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
my $action='cached';

my $initializer = BioMart::Initializer->new( "registryFile" => $confFile, "action" => $action );
my $registry = $initializer->getRegistry;
my $query = BioMart::Query->new( "registry" => $registry, "virtualSchemaName" => "default" );

my @gene_ids;
while(<>) {
	chomp;
	push @gene_ids, $_;
}

$query->setDataset("boleracea_eg_gene");
$query->addFilter("ensembl_gene_id", \@gene_ids);
$query->addAttribute("ensembl_gene_id");
$query->addAttribute("uniprotsptrembl"); # we want UniProtKB identifiers to use in the AgriGO enrichment test
$query->addAttribute("go_id");
$query->addAttribute("name_1006");
$query->addAttribute("definition_1006");
$query->addAttribute("namespace_1003");

$query->formatter("TSV");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################
