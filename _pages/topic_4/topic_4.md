---
title: "Topic 4 - Sequence Alignment"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---

## Recorded lecture

## Accompanying material
* Slides 2020: [UBC - De novo Assembly 2021](https://github.com/UBC-biol525D/UBC-biol525D.github.io/blob/master/Topic_4/Assembly2021Julia%20%20-%20JK.pdf)
* Background reading: The present and future of de novo whole-genome assembly [Paper](https://academic.oup.com/bib/article/19/1/23/2339783?login=true#119542667). 

Programming Resources (in this tutorial or from lecture)
* [Flye Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) and [Flye Paper](https://doi.org/10.1038/s41587-019-0072-8).
* [Hslr Manual](https://github.com/vpc-ccg/haslr) and [Hslr Paper](https://www.cell.com/iscience/fulltext/S2589-0042(20)30577-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004220305770%3Fshowall%3Dtrue).
* [Spades Manual](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md) and [Spades Paper](https://cab.spbu.ru/files/release3.15.3/manual.html).
* [HiFiASM manual](https://github.com/chhylp123/hifiasm) and [HiFiASM paper](https://www.nature.com/articles/s41592-020-01056-5)


## Code break questions

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


## Tutorial 
