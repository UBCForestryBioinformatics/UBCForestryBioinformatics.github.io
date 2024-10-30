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

### If you don't know the commands, call one of the instructors over to walk through the steps you would take to complete the exercise



<details>
<summary markdown="span">**Click for the solution!**
</summary>
```bash
   > cat /mnt/data/anno/SalmonAnnotations.gff | awk 'BEGIN {OFS = "\t"};{if ($3=="gene") print  $1,$4-1,$5}'

# The program `bedtools` would come in handy to deal with duplicate entries

```
</details>
