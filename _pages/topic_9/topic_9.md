---
title: "Topic 9 - SNP Filtering and Analysis"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---


Take a look at this. We've ran the exact same pipeline from read data simulated to a higher coverage. Check out how important coverage is for identifying high quality SNPs.

![](stats by coverage.jpeg)


### Coding challenge
* Use command line tools to extract a list of all the samples in your VCF file, from the vcf file itself. They should be one name per line.
* Take the original vcf file produced and create a vcf of only biallelic SNPs for P1 samples. 
* Use bcftools to filter your vcf file and select for sites with alternate allele frequencies > 0.01, including multi-allelic sites. 
