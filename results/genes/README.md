The folders below this one contain the consensus sequences for the CDSs of interest as identified previously. Since the bulks consist of 
F1 offspring from a cross between (putatively) homozygous parent TO1000 and Jersey Kale, any variants with respect to the reference
(i.e. TO1000) can have an allele frequency of at most 50% if the Kale is homozygous for the alternative, or 25% if it's heterozygous.

Here we compute the consensus for those SNPs where the Kale is homozygous:

```bash
for gene in $genes; do 
	for group in $groups; do 
		get_seq.pl \
			-gene ${gene} \
			-readgroup ${group} \
			-mincover 0.5 \
			-vcffile /home/ubuntu/data/BSA/${group}/${group}_pe.sorted.bam.RG.vcf.gz \
				> ${gene}/${group}.fasta
	done 
done
```
