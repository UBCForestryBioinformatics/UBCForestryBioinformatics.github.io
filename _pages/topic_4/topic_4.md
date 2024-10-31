---
title: "Topic 4 - Sequence Alignment"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---

## Accompanying material


## 1. k-mers  

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
<summary markdown="span">**Answer**
</summary>
```bash
for i in {1..52} 
do 
k=$(($i+8))
cut -c $i-$k /mnt/data/codebreaks/kmer.fa
done
```
</details>



______________


## 2. Nucleotide BLAST

Very frequently in high-throughput bioinformatics we are interested in aligning fragments of DNA against a large reference genome. In this context, our sequencing reads (i.e. the fragmented DNA) is referred to as a *query sequence* that we want to compare with a *subject sequence* (the reference genome).

Select one of the k-mers from the list you generated above.  

<details>
<summary markdown="span">**If you didn't finish the above bit click here to get a 9-mer to work with**
</summary>
```ATCGCACAA```
</details>

With your 9-mer, identify regions of the reference genome that contain that sequence


would like identify regions in the genome where the fragmented DNA (e.g. short or long reads) exhibit good alignment with the reference genome.


If we made a database that contained all the 9-mers present in a species' reference genome, we could index the genome using hash tables (or something like that) to quickly look up where occurences of 

By indexing each sequence with respect to the 9-mers that you listed into a hash index  


Now we're going to shift gears a little and actually align some sequences. 

Quickly count the number of reads in the ```SHORT_READ_FASTQ``` file. How long does it take to BLAST one of those files against the Salmon reference genome? For many eukaryotes, with genomes that are >10-100 times the size of the simulated genome we're working with here, we would have that many more sequencing reads for comparable coverage as well as a larger reference genome to compare to. For example, in the human genome (approx 3Gbp) 30x coverage is acheived with ~300,000,000 paired-end reads. Roughly speaking how long would it take to BLAST that many reads against the 

