---
layout: post
title:  "Assignment"
date:   2024-10-16 12:48:07 -0700
categories: jekyll update
---

### The genetic basis of a flesh colour polymorphism in Chinook Salmon


**The slides from the lecture that go over the Assignment are at the end of the slide deck for [Tutorial 9](/pages/topic_9/topic_9.pdf)**


As the final part of the workshop, you are asked to complete an analysis that weaves together the various concepts from across the week. 

Your colaborators have sequenced the genomes of 22 individual Chinook Salmon (11 with red flesh, 11 with white flesh). Each individual has been sequenced to approximately 10x coverage using Illumina HiSeq paired-end reads. The FASTQs have been trimmed and screened for quality. It is your job to build a pipeline to analyse these genomic data and to help identify candidate loci underlying the flesh colour polymorphism.

________________

The assignment is to build a pipeline to use sequence data to identify the causal gene.
* 75% of your mark comes from building the pipeline - using comments to show that you understand the different steps
* You get the remaining 25% of your mark if the pipeline works on the server you are assigned to
* Extra credit will be given if you correctly identify the causal gene underlying the flesh colour polymorphism


### Assignment Data

*The data for the assignment has been deposited at the following location on the servers:*
 ```/mnt/data/Assignment/fastq/```

The file names of the FASTQ files indicate whether the fish had white or red flesh.


In addition, you'll need the reference genome and the annotations for the Salmon Genome. Those are located here:
* ```/mnt/data/Assignment/SalmonReference.fasta```
* ```/mnt/data/Assignment/SalmonAnnotations.gff```


### Calulating F~ST~ from VCF Files

Recall that we are aiming to build a script to identify the genetic basis of the colour polymorphism. You can take any approach you wish, but I'd suggest sticking with an F~ST~ based analysis.

There are a variety of different software packages implementing  F~ST~, we have installed ```vcftools``` onto the server for you to use. Once you have your final, analysis ready VCF you will need a plain text file containing a list of the samples belonging to each population (in this case white or red-fleshed fish). Assuming you have the files:

* VCF: ```my_lovely_SNPs.vcf.gz```
* Red fish sample list: ```red_fish.txt```
* White fish sample list ```white_fish.txt```

You can run VCF tools to calculate F~ST~ as follows:

```bash
vcftools \
	--gzvcf my_lovely_SNPs.vcf.gz \
	--weir-fst-pop red_fish.txt \
	--weir-fst-pop white_fish.txt \
	--out red_white_comparison_10kb \
	--fst-window-size 10000 \
	--fst-window-step 10000
```

Inspect the output from that file, particularly the file called: ```red_white_comparison_10kb.windowed.weir.fst```. Using that output construct a figure or find another way to examine the output and identify regions with the highest FST. Once you're satisfied with your analysis, compare the regions you identify with the annotations held in the ```SalmonAnnotations.gff```. Can you identify any genes that look like good candidates for controlling the coulour polymorphism?

### Hints

* You do not need to filter SNPs for linkage disequilibrium for this analysis
* [The structure of a GFF file is described here](https://useast.ensembl.org/info/website/upload/gff.html)
* The ```WEIGHTED_FST``` column is what you're after - not ```MEAN_FST```.