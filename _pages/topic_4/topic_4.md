---
title: "Topic 4 - Genome Annotation"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---


# Accompanying material

[Lecture Slides](/pages/topic_4/Topic_4.pdf)

## Clean up your home directory...

Take this opportunity to organize your home directory, the work we've done so far has generated some bits and pieces of data that may be scattered around. Move all the intermediate files into directories named for the different processes they correspond to. 

(If you did this already and are looking at a clear directory you have earned a gold star)

The following UNIX commands will be useful
```
mkdir my_directory # Makes a new directory called "my_directory"
mv file_1 file_2 # Renames file_1 as file_2
mv file_2 my_directory/file_2 # Moves file_2 into the my_directory directory  
rm file_1 # Deletes file_1
```

Be sure to store your lovely reference genome somewhere obvious where you'll remember it. Feel free to name it something memorable. Here are some examples:

```
cp HiFi.ccs.50.p_ctg.fasta MySalmonReference.fasta

cp HiFi.ccs.50.p_ctg.fasta IOweTomTwentyDollars.fasta

cp HiFi.ccs.50.p_ctg.fasta DownstreamSalmon.fasta

```

The name is totally arbitrary, it just needs to be something that you'll remember down the line. Make sure it's sensible though.

## Install EviAnn from GitHub

You installed FastQC yesterday, but today you're going to install a different program: EviAnn (we discussed it in the lecture)

Take a look at the documentation for EviAnn here: [https://github.com/alekseyzimin/EviAnn_release]{https://github.com/alekseyzimin/EviAnn_release}. On Github, released versions of software get listed on the right-hand side of the page. Check there for the current version. 

```
## I lifted this straight from the EviAnn README - you'll need to fix the version numbers to get it to work with the current version of the software

tar xvzf EviAnn-X.X.X.tar.gz
cd EviAnn-X.X.X
./install.sh

```
## Prepare input files for EviAnn...

To run EviAnn, you need files containing RNA-seq data. These files are located at ```/mnt/data/fastq/RNA/```. **Copy a single set of paired-end reads** to your home directory:

```
mkdir RNA
cp /mnt/data/fastq/RNA/cold_sample_04_1.fq.gz ~/home/booker/RNA/
cp /mnt/data/fastq/RNA/cold_sample_04_2.fq.gz ~/home/booker/RNA/
```

EviAnn requires a list of the RNA data that you're going to use to assemble transcripts and map to the genome. The file should contain a single line per pair of files, with the absolute paths to the files included in the names of each. The two FASTQ file names should be separated by a single space character. 

**You can do this the easy (but annoying) way by writing the file name in nano OR you can do this by writing a little bit of UNIX code to make the file quickly**

**Give your file the name ```reads.txt```**

<details>
<summary>Click here for a hint</summary>

Here are some different ways to generate the configuration file (you'll need to adjust these)
```
ls *_R1.fastq.gz | sed 's/\(.*\)_R1\(.*\)/& \1_R2\2/' > reads.txt
```
</details>


## Run EviAnn on your Salmon genome assembly

With your RNA and EviAnn installed it's time to annotate that Salmon!

``` 
EviAnn-2.0.8/bin/eviann.sh -t 6  \
                        -g SalmonReference.fasta \
                        -r reads.txt \
                        -m 10000 ## Maximum intron size
```


*NOTE* EviAnn takes about an 30 mins to run on our servers using 6 threads, so this gives us a good opportunity to run TMUX or BYOBU. Set up EviAnn to run in a TMUX window and leave it running while we discuss what the program does (we'll come back and look at the output later)

### After a very long time...

If everything worked well, EviAnn should have given you a LONG list of files. The one that we are interested in is ```SalmonReference.fasta.gff``` - if you named your files the same way I have.

__________

### Challenges:
* How many CDS are there in your GFF?
* What is the average length of annotated mRNAs?
* What is the total fraction of the genome that corresponds to protein-coding sequence? 

__________ 
