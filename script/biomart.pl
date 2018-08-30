# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use lib '/usr/local/src/biomart-perl/lib';
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;


# apiExampleRegistry.xml  log4perl.conf  mart_0_4_0_5.xsl  mart_0_5_0_6.xsl  martDBLocation.xml  martURLLocation.xml  mime.types  projectRegistry.xml  registryDBPointer.xml  registryURLPointer.xml  settings.conf  templates

my $confFile = "/usr/local/src/biomart-perl/conf/martURLLocation.xml";

#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

		
	$query->setDataset("boleracea_eg_gene");
	$query->addFilter("ensembl_peptide_id", ["Bo9g186210.1"]);
	$query->addAttribute("ensembl_gene_id");
	$query->addAttribute("ensembl_transcript_id");
	$query->addAttribute("go_id");
	$query->addAttribute("name_1006");
	$query->addAttribute("definition_1006");

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
