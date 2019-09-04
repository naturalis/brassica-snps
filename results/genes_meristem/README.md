The folders below this one contain the consensus sequences for the CDSs of interest as identified in 
[intersection.tsv](intersection.tsv). Since the bulks consist of F1 offspring from a cross between the (putatively) 
homozygous parent TO1000 and Jersey Kale, any variants with respect to the reference (i.e. TO1000) can have an allele 
frequency of at most 50% if the Kale is homozygous for the alternative, or 25% if it's heterozygous.

# Creating the codon alignments

[Here](https://github.com/naturalis/brassica-snps/commit/9978494239cab129923e6205d9112c22ef44b8a0) we computed the 
consensus for those SNPs where the Kale is homozygous:

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

[Here](https://github.com/naturalis/brassica-snps/commit/d02f8abee7d083ee62e430c6caf50c810b1088e0) we 
concatenated the CDS records so that each group has a single FASTA record, in a `combined.fasta` file:

```bash
for gene in $genes; do 
	cd $gene
	for group in $groups; do 
		echo ">$group" >> combined.fasta
		grep -v '>' ${group}.fasta >> combined.fasta
	done
	cd ../
done
```

[Here](https://github.com/naturalis/brassica-snps/commit/d2483df99a9ee9691584129180eb7c4c40de26ee) we appended the spliced sequence of the reference genome:

```bash
REFSEQ=/home/ubuntu/data/reference/Brassica_oleracea.v2.1.dna.toplevel.fa.gz
for gene in $genes; do 
	cd $gene; 
	get_refseq.pl -i cds.intervals -s $REFSEQ >> combined.fasta; 
	cd ../; 
done
```

For genes on the `-` strand we did the reverse complement 
[here](https://github.com/naturalis/brassica-snps/commit/11bb7bd4825e49b6767b78c0eb35bd89921e1f31), and
then aligned all genes (with muscle, default settings) [here](https://github.com/naturalis/brassica-snps/commit/32f3b64233584d1986c63fbb938f78ccfff74961)
