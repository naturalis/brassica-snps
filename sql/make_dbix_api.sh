# documented in supplementary materials
perl -MDBIx::Class::Schema::Loader=make_schema_at,dump_to_dir:./lib -e 'make_schema_at("My::Brassica", { debug => 1 }, [ "dbi:SQLite:dbname=snps.db","sqlite" ])'

