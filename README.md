brassica-snps
=============
Scripts and library code for Gene Ontology-mediated refinement of Bulk Segregant Analyses
-----------------------------------------------------------------------------------------

### Description

This repository contains the code for performing the GO-mediated refinement of BSA results
as described in [10.1101/2021.08.11.455982](https://doi.org/10.1101/2021.08.11.455982).
The repository's organization is inspired by the guidelines of 
[10.1371/journal.pcbi.1000424](https://doi.org/10.1371/journal.pcbi.1000424) and is therefore
intended to be understandable and reusable. Nevertheless, additional documentation has been
added to the subfolders in this repository to further explain the various moving parts.

### Contents

- [conf](conf) - configuration files for BioMart mappings and Circos plots
- [doc](doc) - documentation
- [lib](lib) - library code for object-relational mapping
- [results](results) - intermediate and final output files
- [script](script) - analysis scripts in Shell, Perl and R
- [sql](sql) - scripts for creating and managing database in Shell, Perl and SQL

### License

MIT License

Copyright (c) 2018 Naturalis Biodiversity Center

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
