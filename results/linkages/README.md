SNP selection, mapping linkages to absolute locations
=====================================================

In order for Klaas Vrieling's team to design SNPs it would be useful if we had a sense
of how linkage (centiMorgan) distance relates to absolute locations in the genome:
the ideal set of SNPs should be fairly evenly spread out in linkage distance. To get
a sense of this we used the [supplementary table](12864_2012_4560_MOESM1_ESM.xls) from
[this paper](https://doi.org/10.1186/1471-2164-13-523) as an input.

Using [this script](../../script/blast_gene_map.pl) we BLASTed each forward and reverse
primer from the spreadsheet for a SNP (not an SSR). The results are [here](linkages_to_abs.tsv).
This table should be interpreted as follows:

- except for the last four columns, all are identical between the spreadsheet and our table,
  but our table has headings that are a bit more amenable to importing into databases and
  statistical software (because there are no spaces or capital letters).
- from the input spreadsheet, only records whose *Marker Type* is SNP were BLASTed. Of these,
  only those hits that were on the expected *LG* (i.e. chromosome) were retained and that were
  exact matches.
- *Position* means the rank order along the chromosome, the first marker being 1, and so on.
- *cM* means the distance in centiMorgan from the first marker.
- the last four columns in our table are the absolute locations in the reference genome of
  the hit

Given there are 322 hits that are fairly evenly dispersed across the genome, we have the
option of omitting two out of three to come up with a list of 100.