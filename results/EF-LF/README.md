## [uniprot.txt](uniprot.txt)

![G' prime](gprime.png)

A list of UniProtKB identifiers of genes that were in the QTL regions 
of the contrast shown here that contained nonsynonymous SNPs. To be precise, the SNPs
are filtered on G' >= 2.5, as opposed to the benchmark value determined for the focal
contrast.

The list of UniProtKB identifiers is used for the online GO term enrichment test here: 
http://bioinfo.cau.edu.cn/agriGO (session: 283547849)

## [enriched.tsv](enriched.tsv)

The list of enriched GO terms, and the UniProtKB identifiers that are annotated with 
these terms. This table corresponds to the subgraph below.

![SEA](enriched.png)

The table [enriched_GO_0003006.tsv](enriched_GO_0003006.tsv) is produced from this table
by only retaining the terms whose p-value <= 0.05 that are subtended by GO:0003006,
`reproductive developmental process`
