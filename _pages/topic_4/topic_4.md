---
title: "Topic 4 - Genome Assembly"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---


# Accompanying material


[Lecture Slides](/pages/topic_4/Topic_4.pdf)


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

# 1. Check Out Our Sequencing Data!

We just got back our PacBio HiFi reads from the sequencing facility! Let's see how it turned out. First lets specify variables for the paths to our data since we'll be working with them alot...

```bash
hifiDir="/mnt/data/fastq/hifi/"
```

Lets break this down into steps.

1) Subset only lines containing read information from our fastqs
```bash
cat $hifiDir/mergedHiFi.ccs.20.fastq.gz | head # why doesn't this work?
zcat $hifiDir/mergedHiFi.ccs.20.fastq.gz | grep "^@" -A1 #what  does the A option do?
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


Nice. While theres quite a lot of variance in our long-read lengths. With HiFi reads, we are typically aiming for >10,000bp 

Another thing we might want to know is how much coverage our reads give us for each locus in our genome. We are working with a genome size of 10Mb. What is the expected estimate of coverage/bp in our long reads?


```bash
READLENGTH=`cat longread_lengths.txt | awk '{ total += $2 } END { print total/NR }' `
NUMREADS=`wc -l longread_lengths.txt | cut -d" " -f1`
MEANCOV=`echo "$READLENGTH * $NUMREADS / 10000000" | bc`
echo "mean coverage = $MEANCOV" 

```

Both the average read length and expected coverage are good things to know for setting expectations for your genome assembly!

**What coverage should we expect from the ```$hifiDir/mergedHiFi.ccs.50.fastq.gz``` file?**

# 2. Genome Assembly

Genome assembly can take a long time. Because our course is short, we won't have time to explore lots of different options for assembly. So for today, let's just run  ```hifiasm``` and see how we go assembling our salmon genome from the hifi data.

We will compare this to an assembly constructed with short reads a little later on.

### Running HiFiasm 

Now let's take our HiFi data and use ```hifiasm```:

```bash
# Make a location for the genome assembly to live
mkdir hifi_genome_assembly
/softwareDirectory/software/hifiasm \  
	-t 2 \ 
	-o hifi_genome_assembly/HiFi.ccs.20 \ 
	$hifiDir/mergedHiFi.ccs.20.fastq.gz
 ```

*Did you make sure to check that the directory where hifiasm is located exists? (check the ```/mnt/bin/``` directory...)

The above step should take about two minutes or so. When it's done, take a look at your  ```genome_assembly/``` directory to see the files that were generated.

Take a look at the [documentation for hifiasm to get an idea of what is contained in the different output files.](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output)

### Extract primary assembly from hifiasm output

As you can see, ```hifiasm``` has generated a phased genome assembly. However, the files are not in FASTA format, they are in the ```*.gfa``` format - which stands for Graphical Fragment Assembly format. In this format, contigs are stored in tab delimited files with a different format from FASTA. Specifically, ```*.gfa``` files contain information on how segments overlap with each other. [The full specification of these files is contained here](https://gfa-spec.github.io/GFA-spec/GFA1.html). However, for most downstream purposes we'll want a FASTA file that contains just the contigs. *Note that ```hifiasm``` uses version 1 of the GFA format.*

The contigs present in the primary assembly constructed by ```hifiasm``` are stored in the file ```hifi_genome_assembly/HiFi.ccs.20.bp.p_ctg.gfa``` (if you altered the name of the assembly from the code above you'll need to account for that). The contigs themslves are stored on "segment" lines of gfa files, so we'll need to grab them to make a FASTA. 

In GFA files, lines beginning with ```S``` are "segment" lines and have the following structure... 


| Column | Field        | Type      | Regexp              | Description
|--------|--------------|-----------|---------------------|------------
| 1      | `RecordType` | Character | `S`                 | Record type
| 2      | `Name`       | String    | `[!-)+-<>-~][!-~]*` | Segment name
| 3      | `Sequence`   | String    | `\*\|[A-Za-z=.]+`    | Optional nucleotide sequence

So we'll need to grab the contig/segment name (in the second column) and the sequence name (the third column) to convert the GFA to a FASTA File. 

The following AWK comman can do it very efficiently:

```bash
awk '/^S/{print ">"$2;print $3}' hifi_genome_assembly/HiFi.ccs.20.bp.p_ctg.gfa > hifi_genome_assembly/HiFi.ccs.20.p_ctg.fasta  # get primary contigs in FASTA
```

Take a look at the FASTA file that you just made. How many contigs are there in the assembly you made? What is the total length of the assembly you generated?

___________________

In the ```$hifiDir``` directory you'll also find another file called  ```mergedHiFi.ccs.50.fastq.gz```. Using the steps as a template, construct a genome assembly with hifiasm using that file. What is the coverage in that file? Do you expect the genome assembly to be better or worse than the one we made using the ```mergedHiFi.ccs.50.fastq.gz``` file?

*If you choose to do this be aware that the assembly will take a while longer - about 5 minutes or so.*

# 3. Assessing quality of assemblies

Now that we've constructed a couple of genome assemblies, it would be a great idea to compare them to determine which one should be used for downstream analyses. 

As we mentioned above, we constructed an assembly with short reads from the same individual we got the hifi reads from. However, it takes too long to assemble that for this tutorial, so we'll just give you the resulting assembly. 




*DON'T RUN THIS CODE TODAY*
```bash
#	mkdir spades
#	/software/SPAdes-3.15.3-Linux/bin/spades.py \
#		--pe1-1 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R1.fastq.gz \
#		--pe1-2 ${shortreads}/SalmonSim.Stabilising.p1.1.6400000_R2.fastq.gz \
#		-o ./spades/
```

The output file for that short-read only assembly done with Spades can be found here in the /mnt/data/fasta directory...

```bash
mkdir assemblies

cp /mnt/data/assemblies/spades_shortreadonly.fasta assemblies/

mv hifi_genome_assembly/HiFi.ccs.lowCoverage.p_ctg.fasta assemblies/ 
## Double check the names you gave your assemblies and alter this text accordingly! 
mv hifi_genome_assembly/HiFi.ccs.highCoverage.p_ctg.fasta assemblies/
```
_______________

We're going to use ```Quast``` to compare our assemblies. ```Quast``` is command line alignment program that has a collection of nice scripts for examining the properties of our assemblies. It produces some basic stats related to the number and length of sequences in our assembly and can even compare different assemblies against each other.

Having a reference genome of a closely related species can really help asses how *accurate* our assemblies are. In our case we have one better - the high quality reference genome from which our reads were simulated. Importantly, this allows us to get not just the number and length of our contigs/scaffolds, but also _completeness_ (possible when you know the target genome size) and _correctness_.

```bash
#it would be nice to know how well our assemblies have captured known genes. Lets extract this information from the true reference genome's gff file and format it as quast expects 
head /mnt/data/anno/SalmonAnnotations.gff

#quast expects the column order as follows: chromosome, gene id, end, start
cd /assemblies
awk 'split($9,a,";") {print $1 "\t" a[1] "\t" $5 "\t" $4}' /mnt/data/anno/SalmonAnnotations.gff | sed 's/ID=//g' > SalmonReference.genes

#what does $1 and $5 represent in the command above? Awk lets you manipulate columns in all sorts of ways - take note of "split". Here we're splitting column 9 by the delimiter ";" and save each index of the split to an array "a" (i.e. a[1], a[2], a[3], etc)

#run quast on multiple assemblies at once
python3 /mnt/bin/quast-5.2.0/quast.py --help #check out the manual
#notice we're not going to worry about the advanced options 

python3 /mnt/bin/quast-5.2.0/quast.py spades_shortreadonly.fasta HiFi.ccs.20.p_ctg.fasta \
	-r /mnt/data/assemblies/SalmonReference.fasta \
	-g SalmonReference.genes \
	-o quast_out
 ```

**If you want to alter the above code to incorporate more of the assemblies held in the ```/mnt/data/assemblies/``` directory, go for it!

____________


Quast makes some nice HTML formatted reports that we can check out to see how the different assemblies compare. Log out of your terminal session and copy the Quast output to your local machine (i.e. laptop).  

```
## You'll need to edit the following with your own username and the appropriate IP
scp -r <username@ip.address>:/home/<usr>/assemblies/quast_out/ ./
```

Open the report.html results file in your browser and explore the outcome of our assembly efforts. Make sure to eventually click the "View in icarus contig browser".

Question: what new information does this report give us?

Which assembly is the best of the ones we have? 

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
