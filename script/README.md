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
