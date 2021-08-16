SQL database creation and code generation
=========================================

This directory contains files to setup a relational database in SQLite and then export the schema 
to [DBIx::Class](https://metacpan.org/pod/DBIx::Class) modules for handy object-relational mapping. 
The overall workflow is as follows:

1. the helper scripts `make_*.pl` are run to pre-process various reference data tables in a 
   format that is acceptable to the schema. These pre-processing steps are minimal, boiling
   down to stripping off header rows, pre-fixing auto-incrementing primary keys to the 
   records, and processing the names of chromosomes such that they become foreign keys
2. the snps.sql script is piped into the sqlite3 executable, creating the schema and loading
   the tables produced in the previous step
3. the `make_dbix_api.sh` script is run to read the schema from sqlite and generate object-
   relational mappings from this, which are written as perl5 modules for later usage in scripts


