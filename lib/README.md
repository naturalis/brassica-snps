Library code
============

This folder and the folders below it contain library code for the Perl5 programming language.
The library code is used by some of the scripts in the [../script](../script) folder, so those
scripts need to 'know' the location of this folder by way of an environment variable. Assuming
the location of a local clone of this repository as a whole is `$BRASSICA_REPO`, the location
of the library code is then specified as `export PERL5LIB=$PERL5LIB:$BRASSICA_REPO/lib`.

Most of the code here was generated automatically from a SQLite database so that Perl scripts
can interact with the contents of that database more easily. For more about that database and
the code generation tooling, refer to the [../sql](../sql) folder.
