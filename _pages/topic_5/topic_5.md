---
title: "Topic 5 - Genome Assembly"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---


## Accompanying material
* Slides 2020: [UBC - De novo Assembly 2021](https://github.com/UBC-biol525D/UBC-biol525D.github.io/blob/master/Topic_4/Assembly2021Julia%20%20-%20JK.pdf)
* Background reading: The present and future of de novo whole-genome assembly [Paper](https://academic.oup.com/bib/article/19/1/23/2339783?login=true#119542667). 

Programming Resources (in this tutorial or from lecture)
* [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) and [Flye Paper](https://doi.org/10.1038/s41587-019-0072-8).
* [Hslr Manual](https://github.com/vpc-ccg/haslr) and [Hslr Paper](https://www.cell.com/iscience/fulltext/S2589-0042(20)30577-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004220305770%3Fshowall%3Dtrue).
* [Spades Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) and [Spades Paper](https://cab.spbu.ru/files/release3.15.3/manual.html).
* [HiFiASM manual](https://github.com/chhylp123/hifiasm) and [HiFiASM paper](https://www.nature.com/articles/s41592-020-01056-5)



## Tutorial 

Your goal for today is to see how well we can assemble a portion of the Salmon reference genome with PacBio HiFi reads and the '''hifiasm''' assembler.

HiFi reads represent the state-of-the art for genome assembly. However, they are fairly expensive and are not necessarily what is going to be available to depend on many factors (i.e. project budget, sample availability etc.). For the sake of today's tutorial, we're just going to focus on HiFi reads, but make sure to look around for the newest and most promising methods for other data types. The methods adopted will also be highly dependent on the architecture of the genome you are trying to assemble (i.e. depending on repeat content, ploidy, etc.). Beyond just running the commands to run an assembler, we're going to focus on the things you should be considering once those commands have finished running. 

You've already spent some time getting familiar with the data we're working with (e.g. fastq files). While our overall aims may be to understand important evolutionary processes, high quality reference genomes will be indispensible for examining patterns of polymorphism across the genome and across populations.

For today's tutorial, we're going tp be working with two sets of HiFi reads, one high coverage and one lower coverage set. 

## First, lets check out our sequencing data!

We just got back our PacBio HiFi reads from the sequencing facility! Let's see how it turned out. First lets specify variables for the paths to our data since we'll be working with them alot...

```bash
highCov="/mnt/data/fastq/hifi/"
lowCov="/mnt/data/fastq/hifi/"
```

Lets break this down into steps.

1) Subset only lines containing read information from our fastqs
```bash
cat $longreads/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz | head # why doesn't this work?
zcat $longreads/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz | grep "chr" -A1 #what  does the A option do?
# don't forget about man to read the manual of the command grep!

# we still have some extra information being printed.
# what could you add to this pipe to only keep the sequence? hint: grep -v -E "match1|match2"
```

2) Get the length of every line (read)
```bash
#pipe the output above (every line = seperate read) to the following awk command. 
awk -F "" '{print NR "\t" NF}'  #from the output, can you figure out what the awk variables NR and NF are? what does the -F "" option do?
#save the output to longread_lengths.txt
```

Instead of just investigating by eye, it would be nice to get some quick stats of the read lengths of these datasets _i.e. what is the average read length for our long read datasets? how much variation is there around the mean?_ \

While R is typically the go to place for statistical analyses, bash can also handle doing some intuitve stats, which will save us the headache of importing data into R.

For example, awk is a super useful tool for working with column and row based analyses. A one-liner in awk can help us calculate the mean read length.

```bash
awk '{ total += $2 } END { print total/NR }' longread_lengths.txt #this gets the mean of column two.
```

**total += $2** <= set the variable "total" equal to the sum of all items in column 2 \

**print total/NR** <= once finished (END), print the column sum after dividing by NR (number of rows) \


We can get the quartiles of read length pretty easily with bash and some simple arithmetic.

```bash
LINECOUNT=`wc -l longread_lengths.txt | cut -d" " -f1`
FIRSTQUART=`echo "$LINECOUNT * 1 / 4" | bc` 
THIRDQUART=`echo "$LINECOUNT * 3 / 4" | bc` 
cut -f2 longread_lengths.txt | sort -n | sed -n "$FIRSTQUART p" 
cut -f2 longread_lengths.txt | sort -n | sed -n "$THIRDQUART p" 
```
**wc -l** <= number of lines \

**cut -d" " -f1** <= keeps only the first column (based on a space delimeter) \

**sed -n "N p"** <= prints (p) the line at value N \


Nice. While theres some variance in our long-read lengths, its nice to see its actually quite consistent. 

Another thing we might want to know is how much coverage our reads give us for each locus in our genome. We are working with a genome size of 10Mb. What is the expected estimate of coverage/bp in our long reads?


```bash
READLENGTH=`cat longread_lengths.txt | awk '{ total += $2 } END { print total/NR }' `
NUMREADS=`wc -l longread_lengths.txt | cut -d" " -f1`
MEANCOV=`echo "$READLENGTH * $NUMREADS / 10000000" | bc`
echo "mean coverage = $MEANCOV" 

```

Both the average read length and expected coverage are good things to know for setting expectations for your genome assembly!

## Genome Assembly

Genome assembly can take a long time. Because our course is short, we won't have time to run these commands in tutorial. Instead, today we'll focus on comparing the quality of the assemblies produced through a few different programs, using short, short + long, or long read only data. 


#### Short read assembly: SPADES - *don't run*
#using a single core this takes ~14 minutes
```bash
#mkdir spades
#/ohta1/julia.kreiner/software/SPAdes-3.
#15.3-Linux/bin/spades.py \
#	--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
#	--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#	-o ./spades/
````


#### Hybrid (short + long) assembly - SPADES & HASLR - *don't run*
#haslr takes ~15 minutes
```bash
#make new dir
#mkdir hybridspades
#run
#/ohta1/julia.kreiner/software/SPAdes-3.15.3-Linux/bin/spades.py \
#	--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
#	--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#	--pacbio ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
#	-t 20 \
#	-o /hybridspades/ 

#mkdir hybridhaslr
#install
#conda install -c bioconda haslr
#run
#haslr.py \
#	-t 20 \
#	-s ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#	-l ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
#	-x pacbio \
#	-g 10m \
#	-o hybridhaslr \
#	--cov-lr 40 #takes 15 minutes
#
```

#### Long read assembly - FLYE - *don't run*
#this took 8 minutes with 20 threads

```bash
#install
#conda install flye

#run flye assuming lower quality pacbio reads
#flye \
#	--pacbio-raw ${longreads}/SalmonSim.Stabilising.p1.3.30k.PacBio.fastq.gz \
#	--threads 20 \
#	-o flye/ 
#	--genome-size 10m
```

Now we have four assemblies, 1 short read (spades), 2 hybrid (spades, haslr), and one long-read (flye).

## Assess quality of assemblies

Lets compare assemblies. First, copy these out of the /mnt/data/fasta directory

```bash
mkdir assemblies
cd assemblies

cp /mnt/data/fasta/spades_shortreadonly.fasta ./
cp /mnt/data/fasta/spades_hybrid.fasta ./
cp /mnt/data/fasta/haslr_hybrid.fasta ./
cp /mnt/data/fasta/flye_longread.fasta ./

#we could have done this in one line
#cp /mnt/data/fasta/*.fasta
```

bbmap is command line alignment program that has a collection of nice scripts for library and sequence quality control. We're going to use its stats.sh script to get at some basic stats related to the number and length of sequences in our assembly.

An important stat is N50: the number of contigs making up 50% of your assembly. 
Relatedly, L50 describes for those largest sequences that make up 50% of your assembly, what the minimum sequence length is.
This becomes more intutitive when you look at the lower table outputted by bbmap (unfortunately there is some inconsitency between these definitions in that L/N usage is often swapped, as you'll see later)

We have several assemblies we want to generate stats for, so lets automate this using a for loop:
```bash
mkdir assembly_stats
for i in spades_shortreadonly spades_hybrid haslr_hybrid flye_longread
do

/mnt/software/bbmap/stats.sh in=${i}.fasta > assembly_stats/${i}.stats

done
```

Question: which genome has the best *contiguity* via N50/L50 metrics? Whats an issue with relying just on contiguity?

Having a reference genome of a closely related species can really help asses how *accurate* our assemblies are. In our case we have one better - the high quality reference genome from which our reads were simulated. Lets test out Quast - a program for getting detail stats when comparing to a closely related high quality reference.

Importantly, this allows us to get not just the number and length of our contigs/scaffolds, but also _completeness_ (possible when you know the target genome size) and _correctness_.

```bash
#it would be nice to know how well our assemblies have captured known genes. Lets extract this information from the true reference genome's gff file and format it as quast expects 
head /mnt/data/gff/SalmonAnnotations_forIGV.gff

#quast expects the column order as follows: chromosome, gene id, end, start
cd /assemblies
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}' /mnt/data/gff/SalmonAnnotations_forIGV.gff | sed 's/ID=//g' > SalmonReference.genes

#what does $1 and $5 represent in the command above? Awk lets you manipulate columns in all sorts of ways - take note of "split". Here we're splitting column 9 by the delimiter ";" and save each index of the split to an array "a" (i.e. a[1], a[2], a[3], etc)

#run quast on all 4 assemblies at once
/mnt/software/quast-5.2.0/quast.py --help #check out the manual
#notice we're not going to worry about the advanced options 

/mnt/software/quast-5.2.0/quast.py flye_longread.fasta haslr_hybrid.fasta spades_hybrid.fasta spades_shortreadonly.fasta \
	-r /mnt/data/fasta/SalmonReference.fasta \
	-g SalmonReference.genes \
	-o quast_out

scp -r <username@ip.address>:/home/<usr>/assemblies/quast_out/ ./
```

Open the report.html results file in your browser and explore the outcome of our assembly efforts. Make sure to eventually click the "View in icarus contig browser".

Question: what new information does this report give us?



## Review some important command line operations we used today.

For loops - this can be a list of items seperated by spaces or you can specify of range numeric values to use a an iterator
```bash
for i in {1..100}
do
echo $i
done
```
or

```bash
for i in 1 2 3 4 5
do
echo $i
done
```

Making new folders `mkdir` \

Renaming files `mv` and moving them  `mv file path/new_file_name` \

Counting `wc -l` \

Find and replace `sed 's/find/replace/g'` \

Printing a particuar row `sed -n "10p" SalmonReference.genes` \

Column means with `awk '{ total += $2 } END { print total/NR }'` \

Assigning variables `shortreads="/mnt/data/shortreads/"` and calling them `echo ${shortreads}` \

Printing and splitting certain columns (specified with $) with awk 
```bash
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}'
#split column 9 by the delimeter ";" and save every split value into the array a
#print column 1, the first value of array a, column 5, and column 4
```

Cut also allows your to print certain columns but the order will be the same as provided in the input. For e.g., compare the following outputs.
```bash
cut -f2,1 SalmonReference.genes 
cut -f1,2 SalmonReference.genes
```
