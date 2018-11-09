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

SNP selection, finding heterozygous SNPs in the Jersey Kale
-----------------------------------------------------------

First, we transform and filter the VCF file to a tab-separated table:

```bash
gatk VariantsToTable \
  -F CHROM \
  -F POS \
  -F REF \
  -F ALT \
  -GF AD \
  -GF PL \
  -GF GQ \
  --variant group-4_pe.sorted.bam.genotypes.vcf \
  --output group-4_pe.sorted.bam.genotypes.tsv 
```

This gives a table like this:

| CHROM | POS   | REF             | ALT | group-4.AD | group-4.PL  | group-4.GQ |
|-------|-------|-----------------|-----|------------|-------------|------------|
| C1    | 33524 | C               | A   | 128,16     | 242,0,5162  | 99         |
| C1    | 33753 | CTT             | C   | 76,126     | 5122,0,2797 | 99         |
| C1    | 33757 | ATGTCCAAGACAGAT | A   | 86,128     | 5088,0,3255 | 99         |
| C1    | 33816 | G               | C   | 0,236      | 10420,710,0 | 99         |
| C1    | 33835 | A               | G   | 0,233      | 10023,699,0 | 99         |
| C1    | 33863 | CA              | C   | 0,217      | 6557,652,0  | 99         |

We are going to filter these according to the following criteria:

- the total depth, which is the sum of the two numbers in the AD column, must be between
  100 and 400, by the same logic we applied [here](https://github.com/naturalis/brassica-snps/blob/master/script/QTLseqr.R#L56-L57)
- the quality, in the GQ column, must be 99, just like [here](https://github.com/naturalis/brassica-snps/blob/master/script/QTLseqr.R#L59)
- the SNP must be heterozygous, which means the middle of the three numbers in the PL
  column must be 0
- the SNP must be a truly single nucleotide, so REF and ALT must both be length 1
