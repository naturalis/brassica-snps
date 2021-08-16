BioMart configuration
=====================

The [XML](martURLLocation.xml) file in this folder is used in conjunction with
[../../script/biomart.pl](../../script/biomart.pl) to translate gene IDs from 
the _B. oleracea_ reference genome onto UniprotKB identifiers that can be used
as input for GO term enrichment tests. Scripts and configuration files such as
this can be generated (more or less) automatically from queries created in the
BioMart graphical web interface so that they can be rerun locally, provided the
dependencies (BioPerl and the BioMart API) are in place.
