---
layout: post
title:  "Extra Challenges"
date:   2024-10-16 12:48:07 -0700
categories: jekyll update
---

For those of you who are working through the workshop and find yourself completing the daily tutorials in plenty of time, we have a set of exercises for you to work on to hone your skills. 



# Exercise 1: Convert a GFF file to BED format

A routine task in bioinformatics is to convert files from one form to another. In this exercise we'll use the Unix command line to take a file from one format and make the necessary changes to put it in another.

For this exercise, pretend that a collaborator of yours has written a piece of analysis software. This analysis examines the total length of protiein-coding genes and does a sophisticated evolutionary analysis with that information. However, your colleague was not very helpful and did not know that the data is currently in GFF format. The program they wrote only takes BED files. Here's the specification of the file format from your collaborator:

_ _ _ _ _ _ _ _ _ _
*My new and improved algorithm for analysing selection in the genome requires, as input, a files with the following characteristics:*
* *One file per chromosome*
* *Chromosomes are referred to using a single integer name*
* *The input **MUST** be in BED format*
* *The program can detect when there is duplicate sequences, but it is **much** faster if it does not have to identify them*
_ _ _ _ _ _ _ _ _ _

For this exercise can you generate BED files that your collaborator's software would accept?

The annotation file for the entire genome is on the server - we suggest you make a local copy and store it in a convenient location. Here's a command to grab the annotations:

```
cp /mnt/data/anno/SalmonAnnotations.gff ./
```

 _________________


# Exercise 2: Genome-Wide Association Study

**Only try after completing all tasks**

In the SNP analysis tutorial we used bcftools to filter snps, and plink & admixture to investigate population structure. We will plot these results to get a better visual sense of this data, but not quite yet! Before we move to R, we are going to run some analyses to infer potential patterns of selection and adaptation in our dataset. 

Remember, we want to test hypothesis about how our populations of Chinook salmon might have adapted to different temperature environments. We could do this a bunch of different ways, but lets try two. The approach we'll take is to see if we can identify loci that show particularly extreme differentiation between environments of varying temperatures. We'll use the Fst statistic to find these loci that diverge in allele frequency between extreme temperatures. We've already identified a polymorphic phenotype across this temperature gradient, so it might be interesting to see how well environmentally associated allelic variation lines up with alleles that associate with our phenotype. We can find these phenotypically-associated alleles with a GWA. 



## Genome-wide Association ##

A nice thing about plink, the program we used in the previous tutorial, is that alot of programs take the .bed/.fam format for input, including the GWA program GEMMA. We're going to use the same VCF we used to infer patterns of population structure. It's amazing how easy (and fast!!) it is to run a GWA, but we have to be careful about the statistical design and interpretation of these type of analyses.

By default, GEMMA knows to take the 6th column of the plink .fam file as the dependent variable. Take a look at our fam file right now using 'less -S vcf/Chinook_GWAS_filtered_fixedsamps.fam'. Darn, it doesn't have any phenotype info! Let's replace the current column with phenotype info in phenos.txt

```bash
#to do this, lets first make sure we the phenotype information we have is in the same order as our samples
comm --help #make sure you know how to interpret the output!
comm  <( cut -d" " -f1 vcf/Chinook_GWAS_filtered_fixedsamps.fam  ) <( cut -d"," -f-1 /mnt/data/vcf/phenos.txt ) #notice only column 2 (rows in common between datasets) prints
#note the subprocessing

#make modified fam file, by pasting the first 5 columns of our fam with the 2nd column of the phenotype file.
#this will make plink happy as now our 6 column = phenotype information

paste -d " "  <( cut -d" " -f1-5 vcf/Chinook_GWAS_filtered_fixedsamps.fam) <( cut -d"," -f2 /mnt/data/vcf/phenos.txt) > vcf/Chinook_GWAS_filtered_fixedsamps.fammod

```
For operations like above, we're lucky that our information we're interested in merging is in the same order (alphanumeric sample order). Otherwise we'd want to incorporate a couple sort commands in there. This is a super important thing to make sure of before doing any merging

```bash
#plink expects the phenotype to be in the -bfile <prefix>.fam, but right now its in <prefix>.fammod. lets do some quick renaming
mv vcf/Chinook_GWAS_filtered_fixedsamps.fam vcf/Chinook_GWAS_filtered_fixedsamps.famnophenos
mv vcf/Chinook_GWAS_filtered_fixedsamps.fammod vcf/Chinook_GWAS_filtered_fixedsamps.fam

#actually running the GWA is as simple as:
/mnt/software/gemma-0.98.1-linux-static \
	-bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	-lm \
	-o Chinook_GWAS \
	-miss .10 \
	-maf 0.01

#gemme provides options for specifying whether we want to exclude sites based on missing data (-miss) or minor allele frequency (-maf).
```

We talked aobut in lecture how population structure can have effects on GWA that are very important to control for. Lets compare the above GWA to one that controls for population structure via the relatedness matrix in a linear mixed model framework

```r
 #first use the option -gk to generate the relatedness matrix 
/mnt/software/gemma-0.98.1-linux-static \
	-bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	-gk \
	-o Chinook_GWAS_filtered_fixedsamps

#run the GWA controlling for relatedness, with the output matrix (*.cXX.txt)
/mnt/software/gemma-0.98.1-linux-static \
	-bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	-k output/Chinook_GWAS_filtered_fixedsamps.cXX.txt \
	-lmm 4 \
	-o Chinook_GWAS_relatedness
```

GEMMA writes results from these analyses to a folder it makes called output/. Lets rename this folder and move it inside analysis/, and because we are impatient we'll just take a quick look at the the min p-values from these two different GWA runs, before we read it into R.

```bash

mv output/ analysis/gwas/
head analysis/gwas/Chinook_GWAS.assoc.txt #lets look at the p_wald values, column 11
head analysis/gwas/Chinook_GWAS_relatedness.assoc.txt #note that the p_wald column for the linear mixed effect GWA is 13

#a quick way to figure out what SNP is most signficant would just be resorting by those columns, taking just the first row (and dropping all the other columns we don't want to see)

sort -g -k11,11 analysis/gwas/Chinook_GWAS.assoc.txt | head -n 2 | cut -d$'\t' -f11  #g tells sort to interpret scientific notation, -k to sort on column 11
sort -g -k13,13 analysis/gwas/Chinook_GWAS_relatedness.assoc.txt | head -n 2 | cut -d$'\t' -f13

#Note that our minimum p-value is much lower when we account for relatedness

```

Now that we've run some analyses to investigate potential patterns of selection, lets move to R to do some statistics and data visualization.
