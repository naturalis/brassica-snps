[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5211375.svg)](https://doi.org/10.5281/zenodo.5211375)

Refining Bulk Segregant Analyses: Ontology-Mediated Discovery of Flowering Time Genes in _Brassica oleracea_
------------------------------------------------------------------------------------------------------------

### Description

This repository contains the code for performing the GO-mediated refinement of BSA results
as described in [10.1101/2021.08.11.455982](https://doi.org/10.1101/2021.08.11.455982).
The repository's organization is inspired by the guidelines of 
[10.1371/journal.pcbi.1000424](https://doi.org/10.1371/journal.pcbi.1000424) and the 
preferred project layouts promulgated by the Perl5 community and is therefore
intended to be understandable and reusable. Nevertheless, additional documentation has been
added to the subfolders in this repository to further explain the various moving parts.

### Contents

- [conf](conf) - configuration files for BioMart mappings and Circos plots
- [doc](doc) - documentation
- [lib](lib) - library code for object-relational mapping
- [results](results) - intermediate and final output files
- [script](script) - analysis scripts in Shell, Perl and R
- [sql](sql) - scripts for creating and managing database in Shell, Perl and SQL

### How to use

The contents of this repository (scripts, library code) are meant to be used in
conjunction with a fixed folder structure of input and reference data. The
[supplementary data files](https://doi.org/10.5281/zenodo.3402201) of the 
[paper](https://doi.org/10.1101/2021.08.11.455982) that describes our methods may 
serve as a template. Alternatively, the following structure may be replicated:

    root ($DATA)
    |-- reference (contains reference genome in FASTA, annotations in GFF3)
    |-- BSA 
    |   |-- $SAMPLE1 (contains $SAMPLE1_R1.fastq.gz and $SAMPLE1_R2.fastq.gz)
    |   |-- $SAMPLE2
    |   |-- $SAMPLE3
    |   `-- $SAMPLE4
    |-- contrasts
    |   |-- $CONTRAST1 (will contain outputs)
    |   |-- $CONTRAST2
    |   |-- $CONTRAST3
    |   |-- $CONTRAST4
    |   |-- $CONTRAST5
    |   `-- $CONTRAST6
    `-- sqlite (contains chromosomes.tsv)

> **Note** in this structure, the root folder can have any name, as long
> as this is specified as the value of the $DATA environment variable.
> Likewise, the samples and contrasts can have any name. In our analysis, they 
> were labelled SAMPLE1=group-1-EF, SAMPLE2=group-2-IF, SAMPLE3=group-3-LF and 
> SAMPLE4=group-5-NF but this can be anything reasonable. For the contrasts, we
> use CONTRAST1=EF-IF, CONTRAST2=EF-LF, CONTRAST3=EF-NF, CONTRAST4=IF-LF,
> CONTRAST5=IF-NF and CONTRAST6=LF-NF. As always, it is advisable to avoid 
> special characters (such as directory separator symbols) and spaces in folder 
> names.      

With the above structure in place, the analysis steps as described in the 
[supplementary methods PDF](doc/SUPPLEMENTARY_METHODS.pdf) can be performed. 
Where applicable, the PDF indicates how various variables need to be updated 
for this to work with other file naming schemes than we followed.

We provide a Docker [container](https://hub.docker.com/r/naturalis/brassica-snps) 
inside of which are all the tools and their dependencies for running the analysis 
steps. That container should be run in interactive mode (i.e. the user 'boots into' 
the container to run a shell session). Additionally, the container must be made
aware of the location of the above folder structure on the host machine such that
it mounts it under `/home/ubuntu/data` inside the container. The command to start
the container should therefore be something like:

    docker run -it -v $DATA:/home/ubuntu/data naturalis/brassica-snps

Where $DATA is the absolute location of the root of the data folder on the host
machine.

### License

MIT License

Copyright (c) 2018-2022 Naturalis Biodiversity Center

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
