---
title: "Topic 4 - Sequence Alignment"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---

# Accompanying material
[Lecture Slides](/pages/topic_4/Topic_4.pdf)


# 1 - k-mers  

1. Write a one liner to find all the overlaps (i.e. the beginning or the end) exactly 4bp in length between CTCTAGGCC and a list of other sequences in the file /mnt/data/codebreaks/overlaps.fa

2. Find all the unique 9mers in a fasta sequence /mnt/data/codebreaks/kmer.fa

This might be a tricky one and there are likely many ways to do this. First try out these commands that might help you complete this task. It might also be helpful to determine how many characters are in the sequence (wc -c).

Hints: test out the following commands:

```bash
cut -c1- /mnt/data/codebreaks/kmer.fa
cut -c1-4 /mnt/data/codebreaks/kmer.fa
```

```bash
wc -c /mnt/data/codebreaks/kmer.fa
```

```bash
for num in {1..10}
do
echo $num >> file.txt
done
```
What do these commands do? Can you use commands like this to find all the kmers in the sequence? 


One approach we could use involves variable assignment. Variable assignment is a super powerful part of bash coding. A variable can be used as a shortcut for long path, or it can even contain the output of a command. The notation is as follows:

```bash
shortpath="/mnt/data/fastq/shortreads/"
ls $shortpath
cmdout=`echo "test"`
echo $cmdout
```

Think about how you could incorporate some basic algebra and variable assignment to solve this problem.

```bash
#one way to do math in bash is by piping to bc (basic calculator).
echo "1+1" | bc
#another way to do arithmitic in bash is through the $(( )) syntax which tells shell to evaluate its contents
echo $((1+1))
```


<details>
<summary markdown="span">**Answer**</summary>

```bash
for i in {1..52} 
do 
k=$(($i+8))
cut -c $i-$k /mnt/data/codebreaks/kmer.fa
done
```
</details>



______________


# 2 - Nucleotide Alignment

### Exact matches

Very frequently in high-throughput bioinformatics we are interested in aligning fragments of DNA against a large reference genome. In this context, our sequencing reads (i.e. the fragmented DNA) is referred to as a *query sequence* that we want to compare with a *subject sequence* (the reference genome).

First, let's make a copy of the reference genome...


```bash

cp /mnt/data/assemblies/SalmonReference.fasta ~/ ## this will copy the reference genome to your home directory

```


Select one of the k-mers from the list you generated above.  

<details>
<summary markdown="span">**If you didn't finish the above bit click here to get a 9-mer to work with**</summary>

```ATCGCACAA```
</details>

With your 9-mer, identify regions of the reference genome that contain that sequence using ```grep```.

Hint:
  * Try looking at the options for ```grep``` - or try piping the output into ```wc```

*Can you think why your number is probably an underestimate?*

<details>
<summary markdown="span"> **Answer**</summary>

```bash
grep -c my_kmer

# or

grep my_kmer | wc -l
```
</details>



Ok, great. So we've found parts of the genome where our 9-mer is present. Is this an efficient way to build alignments? Why not? What could go wrong?


___________________

What would we do next? Perhaps we should build alignments around these matching kmers. 

### NCBI BLAST

The **B**asic **A**lignment **S**earch **T**ool (**BLAST**) is one of the cornerstones of modern biology. Having a search engine with which to query any sequence against a database of all publicly available sequences is tremendous. BLAST has had a profound effect on biology - the original paper describing the tool has more than 100,000 citations!  

First, let's format the ```SalmonReference.fasta``` as a Blast Database...

```bash
mkdir Tutorial_4
cd Tutorial_4


 /mnt/bin/ncbi-blast-2.16.0+/bin/makeblastdb -in ../SalmonReference.fasta -out SalmonBlast -dbtype nucl
# This should take a couple of seconds and have generated the following eight files: 
###   SalmonBlast.ndb  SalmonBlast.nin  SalmonBlast.not  SalmonBlast.ntf
###   SalmonBlast.nhr  SalmonBlast.njs  SalmonBlast.nsq  SalmonBlast.nto


```




Now try BLAST-ing the sequence we used to build k-mers against the BLAST database you just made. Here's how you can do that:

```bash

 /mnt/bin/ncbi-blast-2.16.0+/bin/blastn -db SalmonBlast -query /mnt/data/codebreaks/kmer.fa  -out kmerSearch.out  
## blastn is the name of the nucleotide BLAST program...

```

Check out the contents of the file called ```kmerSearch.out```. I'd suggest using the ```less``` command.

Did you find a good match?

______________

# Excercise 

### Aligning short reads using BLAST

You may have guessed it, but I took the sequence in the ```kmer.fa``` directly from the reference genome itself. That was just for illustration. What does it look like when we align reads to the reference?

Your task is to BLAST the sequences for the three reads contained in the file called ```three_reads.fq.gz``` located in ```/mnt/data/codebreaks```. 

First, you will need to convert the FASTQ file to a FASTA file. Remember the structure of these files...


<details>
<summary markdown="span">**If you're struggling with the file conversion, here's a little helper:**</summary>

```bash
zcat /mnt/data/codebreaks/three_reads.fq.gz | head -n12 | sed -n '1~4s/^@/>/p;2~4p'  > OUTFILE.fasta

## This is just one way of doing it - there are plenty of others. If you came up with a different way, share it with the group!

## dDon't worry if the above is baffling to you - it's a pretty complex line. 
## If you're struggling with it, try taking it apart piece by piece to see what it's doing...

```

</details>

Once you have a FASTA file, you can query it against the BLAST database as we did before. Do you find matches now? Do they look reliable? 

How long does it take to BLAST one of those files against the Salmon reference genome?


Hint:
  * Prepending the ```time``` program to your lines is a neat way of getting timings from Unix commands

______________

Would this is a feasible strategy for mapping short reads from a whole-genome-sequencing experiment? 

For many eukaryotes, with genomes that are >10-100 times the size of the simulated genome we're working with here, we would have that many more sequencing reads for comparable coverage as well as a larger reference genome to compare to. To sequence the human genome (approximately 3Gbp) to 30x coverage with 150bp paired end reads, how many reads would you need? How long would it take to map these using BLAST? Can you make a (very) rough estimate with your commands above?

