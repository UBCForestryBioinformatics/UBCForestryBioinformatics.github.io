---
title: "Topic 3 - Sequence Data"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---

### Accompanying material

[Lecture Slides](/pages/topic_3/Topic_3.pdf)

# 1 - Overview of sequence data

In this exercise, let's pretend that you have been notified that your data is ready for download from the sequence centre. You submitted a short read sequencing library and a long read sequencing library generated using Illumina and PacBio technologies, respectively.

To start, make a local copy of the data files in your home directory. Run the following:
```
cp /mnt/data/Tutorial_3_data.tar.gz ./ # Make a local copy of the sequence data

tar -zxvf Tutorial_3_data.tar.gz  # Extract the data from the tarball that it was shipped in

```

A "tarball" is a way of compressing an entire directory full of files so that you only have to download one data packet rather than many. The options given to tar do the following:
* ``-z``, this tells the program that the file "Tutorial_3_data.tar.gz" is compressed or zipped
* ``-x``, this tells the program to extract the files from the archive
* ``-v``, this tells the program to print a log to screen
* ``f``, this tells the program that you are specifying the tar-ball as a file on the command line


### Checking data integrity

Often when downloading large files, we want to perform a check on the integrity of the data downloaded. Perhaps you were downloading a large file and the power went down before you were not able to check the data. More generally, when downloading many large files you may not have the time (or the patience) to manually curate each one to ensure that it was downloaded correctly.

It is good practice to check data integrity when moving files from place to place and there are useful functions for checking data integrity. The two main methods that are used are shasum and md5. shasum is perhaps a bit more common as it is available as standard on MacOS. Both methods generate what is called a “checksum”, a hexadecimal string (e.g. “d241941bac307bd853fd21945d029e62c83cea71”) that is unique to a given file.

The data that you obtained for Topic 1 also contained a file called ```SalmonData_checksums.sha```. If you inspect the contents of ```SalmonData_checksums.sha``` you’ll notice that there are two columns in the file, the left hand column contains a bunch of checksums, the righthand column contains the name of corresponding files. You can compare all the files in the directory you downloaded using:

Store your checksums with your data, where ever you host it.

Navigate to the location of the data and check that the downloaded data actually match those that were sent by the sequencing centre...

```
cd Tutorial_3_data
shasum -c SalmonData_checksums.sha
```

As long as everything worked well, you should have seen reassuring messages print to screen after invoking the `shasum` command.

______________________________

Now we've got our data, let's just make a couple of directories to tidy things up.

```
mkdir short_reads
mv Salmon.Illumina.R?.fastq.gz short_reads/ # Note the use of the single character wildcard!

mkdir long_reads
mv Salmon.PacBio.fastq.gz long_reads/

cd ../

```

Great, that's nice and tidy now.


# Exercise 1

A question that might be on top of mind would be, how many reads are we working with?

Can you use command line tools to get a count of the number of reads that we have from the two datasets?

*Hint, remember the structure of a fastq file*

Once you've done that, take a look at the contents of each file. How are paired-end reads identified? What information is used to identify read mates?

# 2 - Read Quality

### Install fastqc

We are going to use a program called `FastQC` to examine the quality of the sequencing reads that we are working with.

It is really easy to generate a report using `FastQC`, it's just a matter of telling the program where the files are that you want to analyse. For example:

```bash

shortreads="/home/tommyB/Tutorial_3_data/short_reads/" # Obviously use your own user name here
longreads="/home/tommyB/Tutorial_3_data/long_reads/"

mkdir qualityControl
cd qualityControl

#The program is installed on the server so can be accessed as follows:
fastqc ${shortreads}/*fastq.gz ${longreads}/*fastq.gz -o ./
```

*Oh no!* `FastQC` was not found!If you read the error message it looks FASTQC has not been installed on these servers. I guess this would be a good chance to try downloading a package from the internet.

The first thing you need is the address of the most recent version of `fastqc`. If you navigate the program's website, you'll find several download links when you click "Download Now":
[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

*DO NOT LEFT CLICK ON THE DOWNLOAD LINK!*
Instead, right click the link for Linux downloads and hit "Copy link address"

Now, we can use that link to download the package directory to our local directories:
```bash
wget <insert copied link here>

# This should have downloaded a zipped file called something like fastqc_v0.11.9.zip

# To unzip just run:
unzip fastqc_v0.11.9.zip ## Note that the version umber has probably changed...

# This will have made a directory, go into it:
cd FastQC

## If you get the following message, FASTQC is not executable:
Can't exec "java": No such file or directory at ./fastqc line 350.

## You can change that using:
chmod +x fastqc

```

Phew, now we should have a happy version of fastqc that we can work with.

### Now lets look at some properties of our reads

```bash

shortreads="~/Tutorial_3_data/short_reads/"
longreads="~/Tutorial_3_data/long_reads/"

mkdir qualityControl
cd qualityControl

#fastqc is a widely used program for looking at read quality
 ~/FastQC/fastqc  ${shortreads}/*fastq.gz ${longreads}/*fastq.gz -o ./

# Now, we are referring to the version of FastQC that we just downloaded, and it should play nicely with the Java version we are working with.

#some programs allow you to input multiple files at once using wild cards. this is convenient.
#otherwise we would have had to specify each file individually. Here's another way we could use wildcards here:

#fastqc ${shortreads}/*R?.fastq.gz ${longreads}/*fastq.gz -o ./ )

#summarise into a single report
multiqc ./
```

If the above looked like it ran happily, logout of your current session by writing:
```logout```


Download the output of multiqc run as follows:

```bash
#From a terminal window running on your local computer
scp <username@ip.address>:/home/<username>/assembly/multiqc_report.html <path on your computer where you want the file>
```

Open that file up in an internet browser and have a look.

Since we are using simulated data, our data does not have adapters. But taking a look at the signatures of our highly quality reads can help you get a sense of what very high quality data might look like, and additionally,by comparing to your real data, how sequencing technology influences data quality.  

**Discussion Question: our data effectively has had adapters removed from reads and already been trimmed for low quality BPs near the end of reads. what might be the cons of read trimming?**

__________________

# Exercise 2

Let's pretend that a collaborator just got in contact with you. Going through a shel midden they found some mouldy old bones that look kind of like a Chinook Salmon. They managed to extract DNA from these shells and sequenced them. Your collaborator has sequenced the DNA from this mouldy fish and have sent you the FASTQ files. 

For this exercise, take the ```mouldyFish.fastq.gz``` file and run it throught ```fastqc```, download the output an take a look at it on your browser. What differences are there between these reads and those that you looked at earlier?


