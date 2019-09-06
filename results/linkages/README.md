SNP selection, mapping linkages to absolute locations
=====================================================

In order for Klaas Vrieling's team to design SNPs it would be useful if we had a sense
of how linkage (centiMorgan) distance relates to absolute locations in the genome:
the ideal set of SNPs should be fairly evenly spread out in linkage distance. To get
a sense of this we used the [supplementary table](12864_2012_4560_MOESM1_ESM.xls) from
[this paper](https://doi.org/10.1186/1471-2164-13-523) as an input. (The values in this
table are based on the _capitata_ cultivar.)

Using [this script](../../script/blast_gene_map.pl) we BLASTed each forward and reverse
primer from the spreadsheet for a SNP (not an SSR) against the _Brassica oleracea_
TO1000 reference genome. The results are [here](linkages_to_abs.tsv). This table should 
be interpreted as follows:

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

Given there are 322 hits that are fairly evenly dispersed across the genome, we can look for
SNPs in their vicinity. For example, for each of ±1/3 of these hits, find a good SNP nearby,
for a total of 100 SNPs.

SNP selection, finding heterozygous SNPs in the Jersey Kale
-----------------------------------------------------------

Now we have to find candidate SNPs in the Jersey Kale mapped genome, on which we've done
variant calling, so we have a VCF file to work with.

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

This gives a table like this (but very large):

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

This is implemented in [this script](../../script/filter_candidate_snps.pl). Besides the filtering,
it takes the following steps to prepare the output for ingesting into SQLite:

- strips the 'C' prefix from the CHROM columns, making it a foreign key to the chromosomes table
- pre-pend an autoincrementing integer that becomes the primary key
- only emit primary key, foreign key, POS, REF, ALT

We execute as:

```bash
filter_candidate_snps.pl group-4_pe.sorted.bam.genotypes.tsv > group-4_pe.sorted.bam.genotypes.filtered.tsv
```

Which yields (sans header):

| snp_id | chromosome_id | position | ref | alt |
|--------|---------------|----------|-----|-----|
| 1      | 1             | 1517     | T   | C   |
| 2      | 1             | 21618    | G   | A   |
| 3      | 1             | 21738    | G   | A   |

We [add a table definition](https://github.com/naturalis/brassica-snps/commit/cdf480c8ff40e0b188058f63cc8000bf6b09e2e2),
and, using [this script](../../sql/make_dbix_api.sh), we then
[update the api](https://github.com/naturalis/brassica-snps/commit/46578876f1e4f8914cbcc2ac0f1b14bdd094bbad).
Finally, we import the data:

```sql
.separator "\t"
.import group-4_pe.sorted.bam.genotypes.filtered.tsv kale_snps
create index kale_snp_location_idx on kale_snps(chromosome_id,position);
```

SNP selection, finding SNPs near linkage markers
------------------------------------------------

Now that we have mapped markers from relative to absolute distances (first section), and have filtered
SNPs (second section), we now want to find SNPs near the markers. Using this
[script](../../script/reconcile_linkage_kale_snps.pl) we iteratively build a growing window around
each marker, until we find a SNP. The results of this are [here](kale_snps_near_markers.tsv).
