#!/bin/bash
genes=`cut -f1 intersection.tsv | grep -v 'Gene'`
groups="group-1-EF group-2-IF group-3-LF group-5-NF"
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

