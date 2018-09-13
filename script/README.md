# Contains data exploration scripts and results

## genes_in_qtlregions.pl

Usage:

```shell
genes_in_qtlregions.pl -db <sqlite.db> -c <contrast>
```

For any given contrast (such as `EF-LF`), emits all gene IDs that lie with the QTL regions of that
contrast. To obtain from this [the listing of the genes over all 6 contrasts, sorted by number of 
contrasts in which they occur](genes.txt), one might do this:

```shell
CONTRASTS=`ls /home/ubuntu/data/BSA/contrasts`

# create one big table
for C in $CONTRASTS; do 
   genes_in_qtl_regions.pl -c $C >> genes.tsv
done

# collapse to file with counts
cut -f4 genes.tsv | sort | uniq -c | sort --reverse > genes.txt
```

## biomart.pl

Usage:

```shell
biomart.pl $ID1 $ID2 .. $IDn
```

For any given list of transcript IDs such as may be encountered in the _Brassica oleraceae_ genome annotation 
to a GO term/name/definition using biomart. To get [the table of terms](go_terms.tsv) associated with the genes 
that occur four or more times, one might do this:

```shell
egrep '^\s*[456]' genes.txt | sed -e 's/      . //' | biomart.pl > go_terms.tsv
```

And the [table of terms that match 'flower' anywhere](flowering.tsv) is generated like this:

```shell
grep flower go_terms.tsv > flowering.tsv
```
## snps_in_cds.pl

Usage:

```shell
snps_in_cds.pl $ID1 $ID2 .. $IDn
```

For any given list of transcript IDs such as may be encountered in the _Brassica oleraceae_ genome annotation,
look up which SNPs occur in coding sequences. For those SNPs, look up whether these are synonymous or not. To 
get the [table of SNPs](snps.tsv), one might do this:

```shell
egrep '^\s*[456]' genes.txt | sed -e 's/      . //' | snps_in_cds.pl > snps.tsv
```
To get a list of UniProtKB identifiers (used for [GO enrichment tests](http://bioinfo.cau.edu.cn/agriGO)), one 
might do this:

```shell
for C in $CONTRASTS; do
   # 1. only get the SNPs associated with the focal contrast $C
   # 2. only get the nonsynonymous SNPs
   # 3. retain the column with transcript IDs (first column)
   # 4. sort | uniq => retain distinct IDs
   # 5. lookup the UniProtKB ids and go terms via biomart
   # 6. retain the column with UniProtKB ids (second column)
   # 7. sort | uniq => retain distinct IDs
   # 8. filter out the header and write to file 
   grep $C snps.tsv | grep nonsyn | cut -f1 | sort | uniq | biomart.pl | cut -f2 | sort | uniq | grep -v UniProtKB > ${C}.txt
done
