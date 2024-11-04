---
title: "Topic 9 - SNP Filtering and Analysis"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---


### Accompanying material
* [Slides](./VariantCalling 2021.pdf)

In this tutorial we're going to use SNPs called with GATK to analyse patterns of population structure in the Chinook genome and conduct a GWAS.

### 1. Filtering SNPs for analysis

As we talked about in the lecture, the filters you choose to apply to your data can have a big impact on downstream analysis. Last topic we called variants across the two chromosomes. If you look at the VCF, you'll notice there are a lot of sites only genotyped in a small subset of the samples. This can happen with lower overall read depth (in this case this is whole genome sequencing at ~8X depth), but can be due to other factors like divergence between sample and reference. We also have indels, and SNPs with more than two alleles. Many population genomic analyses have been developed on the assumption that we are looking at two and only two alleles at a site. So first lets filter the VCF to a smaller set of usable sites.

#### NOTE:

* Yesterday we made VCFs for 50 individuals. For the cool analyses we want to run today, we are going to need more samples. We have a large VCF file containing many more individuals, you can copy it to from `/mnt/data/vcf/Chinook_GWAS.vcf.gz` to `~/vcf`

________________


We're going to use _bcftools_ a very fast program for processing and filtering VCF files. Here are what we want to filter for:

* At least 80/100 samples genotyped (which is 160 alleles since they're diploid).
* Variant call quality (QUAL) > 30 - a 99.9% chance there is a variant at the site given genotype calls across samples
* Only one alternate allele.
* No indels.
* At least 2 copies of the alternate allele


```bash
cd ~/vcf
bcftools view \
    -c 2 \
    -i 'INFO/AN >= 160 && QUAL > 30' \
    -m 2 \
    -M 2 \
    -v snps \
    -O z \
    Chinook_GWAS.vcf.gz > Chinook_GWAS_filtered.vcf.gz

#Lets index the vcf file for future use
tabix Chinook_GWAS_filtered.vcf.gz
```

#### Coding challenge

* How many sites remain in the filtered VCF? How many were removed? How many on each chromosome? Don't forget about `grep -v` ad `wc -l`!

________________

A common first pass analysis is to use structure to look at clustering in your data. Admixture is a widely used program but orders of magnitude faster than STRUCTURE (see lecture). We're going use that, but before that we have to convert our VCF to the specialized format. We can do that with Plink - Plink is a large set of tools for manipulating genetic data, running GWAS and calculating various stats. It's geared towards human data, so sometimes you have to force it to work with non-human data. For example, it assumes you have human chromosomes (eg 23 of them) and will complain if it doesn't see them.


```bash
cd ~/
mkdir analysis

#we had a bug in our pipeline that appended some extra characters to the beginning of sample names
#you can check sample names by greping for "#CHROM", which is the first string of the sample header line
zgrep "#CHROM" vcf/Chinook_GWAS_filtered.vcf.gz
#we could use a specialty software like bedtools reheader to fix this, but lets just use basic bash commands

zcat vcf/Chinook_GWAS_filtered.vcf.gz | sed 's/-e Chinook/Chinook/g' | bgzip > vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz
#the key command here is the sed 's/find/replace/g' , zcat is uncompression to standard out and bgzip is recompressing

plink=/mnt/software/plink
$plink --make-bed \
	--vcf vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz \
	--out vcf/Chinook_GWAS_filtered_fixedsamps \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr
```
This produces files with the suffix .nosex, .log, .fam, .bim, .bed. We can use these in Admixture.

NOTE: When inferring patterns of population structure (i.e. admixture/pca) its good practice to filter your VCF for linkage (i.e. remove highly linked sites). If you can't filter for linkage, subsetting the site also helps (i.e. selecting every 10th site). Not only does this make your inference draw from independent data, but it also keeps run times down!

To prune for LD, we'll ask plink to slide across the genome (10 snps at a time),and in windows of 100snps, calculate LD between each snp, removing those with an LD (r2) > .5

```bash
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps \
	--indep-pairwise 100 10 0.5 \
	--out vcf/Chinook_GWAS_filtered_fixedsamps \
	--make-founders \
	--allow-extra-chr
#allow extra chromosome is another way we force plink to work with non-human data (i.e. allow chromosomes that don't have the typical human chr names)

#this produces two files, .in (to include) and .out (to exclude)
#now lets actually extracts the snps that remain after LD pruning

$plink \
--bfile vcf/Chinook_GWAS_filtered_fixedsamps \
--extract vcf/Chinook_GWAS_filtered_fixedsamps.prune.in \
--make-bed \
--out vcf/Chinook_GWAS_filtered_fixedsamps_LDpruned \
--allow-extra-chr

#great! we have our files in the necessary format, now lets run admixture
admixture=/mnt/software/admixture
$admixture vcf/Chinook_GWAS_filtered_fixedsamps_LDpruned.bed 2
```
Uh oh that doesn't work, it produces this error message.
```bash
****                   ADMIXTURE Version 1.3.0                  ****
****                    Copyright 2008-2015                     ****
****           David Alexander, Suyash Shringarpure,            ****
****                John  Novembre, Ken Lange                   ****
****                                                            ****
****                 Please cite our paper!                     ****
****   Information at www.genetics.ucla.edu/software/admixture  ****

Random seed: 43
Point estimation method: Block relaxation algorithm
Convergence acceleration algorithm: QuasiNewton, 3 secant conditions
Point estimation will terminate when objective function delta < 0.0001
Estimation of standard errors disabled; will compute point estimates only.
Invalid chromosome code!  Use integers.
```
Our chromosomes are named chr_1 and chr_2, not integers like admixture is expecting. Like many programs, this is coded for human data where chromosomes are known and numbered clearly. We need to rename the chromosome column of the vcf so that they're integers. In this case, that means removing "chr_" from any line that starts with that, although it would depend on how your chromosomes are named. We'll also have to redo the LD pruning.

```bash
#numeric chr names
zcat vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz |\
	sed s/^chr_//g |\
	gzip > vcf/Chinook_GWAS_filtered_fixedsamps_numericChr.vcf.gz

#sed s/^chr_//g <- removes "chr_" from the beginning of lines

#make new bed from vcf
$plink --make-bed \
	--vcf vcf/Chinook_GWAS_filtered_fixedsamps_numericChr.vcf.gz \
	--out vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr

#redo LD analysis
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--out vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--allow-extra-chr \
	--make-founders \
	--indep-pairwise 100 10 0.5

#make new bed with snps in low LD
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--extract vcf/Chinook_GWAS_filtered_fixedsamps_numericChr.prune.in \
	--make-bed \
	--out vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned \
	--allow-extra-chr

#run admixture
$admixture --cv vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.bed 2 |\
tee analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.out; \

#NOTE: "tee" takes the output of a command and saves it to a file, while
# also printing letting it print to the screen. So we can watch the progress while also
# saving that output.
```
This works! With 100 samples and ~31000 SNPs across two chromosomes it finishes in less than a minute. We only ran it for one value of K (2) but we should also test different K values and select the best K value. Common practice is to run from 1 to number of populations (10). Lets do that but skip some in between. This will take a couple minutes.

```bash
for K in 1 3 4 10; \
do
$admixture --cv vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.bed $K |\
tee analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.${K}.out; \
done

#now move the admixture output files to the analysis folder
mv Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*.P Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*.Q analysis/

```
The best K value for Admixture is typically the K value with the lowest cross-validation (CV) error. The CV error are in the .out files we saved. One easy way to look at all those scores is to print all the .out files and then keep only the lines that include "CV" using grep.


```bash
cat analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*out | grep CV
# CV error (K=1): 0.49626
# CV error (K=10): 0.79357
# CV error (K=2): 0.51421
# CV error (K=3): 0.54330
# CV error (K=4): 0.57121
```


This shows that the lowest CV error is with K=1, but actually K=2 is a close second. To see how this lines up lets look at the .Q file, which shows group assignment for each sample. The Q file doesn't include sample names so we can put those together using "paste. One way to asses if there is meaningful population structure is to check whether each population grouping has individuals of major ancestry.

```bash
paste <(cut -d" " -f1 vcf/Chinook_GWAS_filtered_fixedsamps.fam) analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.Q

#what was that <( ) notation?
#this is another cool bash trick called *subprocessing*.
#Instead of making an intermediate file of sample names to join with our K values...
#...we can use subprocessing to take the output of the cut command and pass that to paste.

#Chinook.p1.i0	0.999990 0.000010
#Chinook.p1.i1	0.999990 0.000010
#Chinook.p1.i2	0.999990 0.000010
#Chinook.p1.i3	0.999990 0.000010
#Chinook.p1.i4	0.999990 0.000010
#Chinook.p1.i5	0.999990 0.000010
#Chinook.p1.i6	0.999990 0.000010
#Chinook.p1.i7	0.999990 0.000010
#Chinook.p1.i8	0.999990 0.000010
#Chinook.p1.i9	0.999990 0.000010
#Chinook.p10.i0	0.000010 0.999990
#Chinook.p10.i1	0.000010 0.999990
#Chinook.p10.i2	0.000010 0.999990
```

While K=1 has the lowest cross validation error, clearly we can see that individuals from population 1 and population 10 (furthest apart in sampling space) beling to different groupings.

We're going to plot these results, but before we leave the command line, lets also one more analysis. Lets run a *PCA*. This is a very nice model free approach to visualizing population structure, and is a good complement to model based structure analyses, like admixture.

We can do this with just one line of code in plink.

```bash
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr \
	--pca \
	--allow-extra-chr \
	--out analysis/Chinook_GWAS_filtered_fixedsamps_numericChr
```

This was on our full dataset, but we learned above it might be a good idea to prune for linkage. Lets redo the analysis, but now with the pruned set...

```bash
#we already have an LD pruned bed, so subbing in that input file
$plink \
	--bfile vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned \
	--pca \
	--allow-extra-chr \
	--out analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned
````

___________________

### 2. Patterns of Population Structure 


We're working with Rstudio on our desktops, so download the "vcf" and "analysis" directories to your laptop. The rest of this tutorial should be run in your Rstudio IDE.

NOTE: This tutorial is based on Rstudio 1.2.1335 and R 3.6.1, the latest version of both. Almost all steps should work identically on older versions, but there may be issues installing some packages. In this case, I recommend updating your version of R unless you have a specific reason not to.


The first step to any organized R project is to create a new Rstudio project. A project keeps all your different scripts and results together in a single directory. It also separates saved variables, so when you are switching between different projects you aren't accidentally using the same variables between them.

In Rstudio, select File->New Project

Then select "New Directory".

![](rstudio_project_1.jpeg)

Click on "New Project".

![](rstudio_project_2.jpeg)

Enter directory name "biol525d" and put it somewhere you can get to. In my case, I put it on the Desktop directory. Finally, click "Create Project".

![](rstudio_project_3.jpeg)

After it has been created, move your "analysis" and "vcf" directory into your biol525d project directory so you have easy access to those files.


You're now in your Rstudio Project and the next step is to install the tidyverse package, which includes a suite of tools for manipulating data and plotting. The examples today will be focused on tidyverse solutions to problems. One key feature of the tidyverse world is the use of "%>%" to pipe data between commands. This functions similar to "\|" in the commandline and helps put together strings of commands that are work together and are easy to understand.


``` r
#install.packages("tidyverse")
library(tidyverse)
```

We calculated cluster assignment for each sample at different values of K. Lets try loading up a single output file. We'll also need to keep track of the sample names.

``` r
samplelist <- read.table("analysis/Chinook_GWAS_fiiltered_fixedsamps_LDpruned.fam", header=F)
class(samplelist) #its useful to keep checking how your data is being treated - here its a DF

#all of these columns except the first aren't useful to us - all we need is the order of sample names. 
samplelist <- samplelist[,1] #this is called *INDEXING*
class(samplelist) #after subset, its just a vector of "character" sample names, so its no longer a DF

#read in admixture results
read_delim("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.Q",
                  col_names = paste0("Q",seq(1:2)),
                  delim=" ")
```

For this last command, we're using the _read\_delim_ function which requires column names. Since the data file doesn't have any column names, we have to specify them and we're using a combination of paste and seq to produce "Q1", "Q2". We could have hard coded it c("Q1","Q2") but this way it works for an arbitrary number of columns just by changing the second value in seq().

Now we could work on each value of K individually, but its easier to load all of them at once. One problem is that they each have different numbers of columns. The solution is converting from a wide to long format. In a long format, each row has a single data point instead of multiple. The tidyverse tools are set up to prefer long data ((and there are other reasons)[https://sejdemyr.github.io/r-tutorials/basics/wide-and-long/#a-case-for-long-data]) so lets do that.

Its possible to load multiple files in a single call, but for transparency lets use a loop. We first make an empty dataframe that we're going to fill, then we loop through our output files, convert them to long format and add them to the master set.

```r
all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in c(1:4,10)){
  data <- read_delim(paste0("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.",k,".Q"),
                  col_names = paste0("Q",seq(1:k)),
                  delim=" ")
  data$sample <- samplelist
  data$k <- k

  #Each time a new .Q file is read in, the following step converts it from wide to long format.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data) #append this iteration to the full dataset
}

all_data
# A tibble: 2,500 x 4
#   sample            k Q     value
#   <chr>         <dbl> <chr> <dbl>
# 1 Chinook.p1.i0     1 Q1        1
# 2 Chinook.p1.i1     1 Q1        1
# 3 Chinook.p1.i2     1 Q1        1
# 4 Chinook.p1.i3     1 Q1        1
# 5 Chinook.p1.i4     1 Q1        1
# 6 Chinook.p1.i5     1 Q1        1
# 7 Chinook.p1.i6     1 Q1        1
# 8 Chinook.p1.i7     1 Q1        1
# 9 Chinook.p1.i8     1 Q1        1
#10 Chinook.p1.i9     1 Q1        1
# … with 1,990 more rows

```

Now to plotting. We first try plotting just K=2 by filtering using the filter() command.
```r

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack")
```

![Image](/pages/topic_9_pics/structure_1.png)

Hurray, it works! Although the base plot is pretty ugly. Lets fix a few things:
* *mutate(sample=gsub("Chinook.p","",sample))* <= remove the redunant "Chinook.p" from sample labels
* *xlab("Sample") + ylab("Ancestry")* <= Change the axis labels to be more accurate
* *theme_bw()* <= Remove the grey background.
* *theme(axis.text.x = element_text(angle = 60, hjust = 1))* <= Rotate the x-axis labels so they don't overlap
* *scale_fill_brewer(palette="Set1",labels=c("1","2"),name="K")* <= Change the fill color and legend labels

``` r
all_data %>%
  filter(k == 2) %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))

```

![Image](/pages/topic_9_pics/structure_2.png)

This is much better, but its still annoying that population 10, comes between population 1 and 2 (by default R won't sort vectors numerically if it has a character in it). We can get a little hacky with this and mutate the "sample" column one more time, removing the i character and telling R to explicitly treat sample as a factor (but first sorted numerically).
```R

all_data %>%
  filter(k == 2) %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  mutate(sample=as.factor(as.numeric(gsub("i","",sample)))) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1","2"))

```


![Image](/pages/topic_9_pics/structure_3.png)

We can also plot all the different K values together using facet_wrap().

```R

all_data %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  mutate(sample=as.factor(as.numeric(gsub("i","",sample)))) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette="Spectral") +
  facet_wrap(~k,ncol=1)
  
```

![](structure_4.jpeg)

That looks awesome. Make sure to take some time to look at how modeling an increasing number of ancestral populations leads to different inferences about population structure. 

Remember though, K=1 (no structure) was our most likely scenario! This actually isn't what we might have inferred from looking at the structure plot alone...at K=2 to 4, we can find individuals that fit in majority to all of the different modelled pops, in a way that might make sense geographically. Hmm. 

Lets take a look at the PCA and compare.


```r
pca_output<-read.table("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr.eigenvec")
head(pca_output)
```

Plotting challenge 1
--------------------

Using some of the tricks we learned above, remove the redundant sample column, and then name the remaining columns (call the first "sample", the second "PC1", the third "PC2", all the way to the last column "PC20"! Save the output to "pca_named"

SOLUTION (check w a neighbour first!):

	pca_named<-pca_output[,-1]
	names(pca_named)<-c("sample",paste0("PC",seq(1:20)))
  {: .spoiler}
  
Now that you have your named column, plot your PCA!

```r
pca_named %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  mutate(sample=as.factor(as.numeric(gsub("i","",sample)))) %>%
  ggplot(aes(PC2,PC1, color=sample)) +
  geom_point() +
  theme_bw()

```
![](pca_1.jpeg)

Woah. We can kinda begin to make up some clustering, but theres something funky going on with the legend. Its giving every individual a distinct color because we've said to color by sample (aes(color=sample)). We would like to color by population though, so lets split up those labels by the "." to get two columns, one for population and one for individual.

```r
pca_named %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  separate(col=sample,into=c("Pop","Ind")) %>%
  mutate(Pop=as.factor(as.numeric(Pop))) %>%
  ggplot(aes(PC2,PC1, color=Pop)) +
  geom_point() +
  theme_bw()

#what did we just do?  
?tidyr:::separate

```
![](pca_2.jpeg)


Right away we can see that while there is _some_ clusering of individausl by population, that there is nearly as much clustering of populations with each other. We might not have expected this from looking only at the structure plot. Lets visualize some other PCs of genotype variation. 

Lets also make use of the actually eigenvalues.

```r
eigenvals<-read.table("analysis/Chinook_GWAS_fiiltered_fixedsamps.eigenval",col.names = "var")
#the eigenvals are the variance explained by each PC. We want to convert these to fraction of the total variance explained.

eigenvals$frac<-round( (eigenvals$var / sum(eigenvals$var)*100), digits=2) #percent of total variance

pc_1v2 <- pca_named %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  separate(col=sample,into=c("Pop","Ind")) %>%
  mutate(Pop=as.factor(as.numeric(Pop))) %>%
  ggplot(aes(PC2,PC1, color=Pop)) +
  geom_point() +
  ylab(paste0("PC1 (", eigenvals$frac[1], "%)")) +
  xlab(paste0("PC2 (", eigenvals$frac[2], "%)")) +
  theme_bw()

pc_1v3 <- pca_named %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  separate(col=sample,into=c("Pop","Ind")) %>%
  mutate(Pop=as.factor(as.numeric(Pop))) %>%
  ggplot(aes(PC3,PC1, color=Pop)) +
  geom_point() +
  ylab(paste0("PC1 (", eigenvals$frac[1], "%)")) +
  xlab(paste0("PC3 (", eigenvals$frac[3], "%)")) +
  theme_bw()
  

pc_2v3 <- pca_named %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  separate(col=sample,into=c("Pop","Ind")) %>%
  mutate(Pop=as.factor(as.numeric(Pop))) %>%
  ggplot(aes(PC2,PC3, color=Pop)) +
  geom_point() +
  ylab(paste0("PC3 (", eigenvals$frac[3], "%)")) +
  xlab(paste0("PC2 (", eigenvals$frac[2], "%)")) +
  theme_bw()


#plot them together
library(ggpubr)
ggarrange(pc_1v2, pc_1v3, pc_2v3, ncol=2,nrow=2, common.legend=T)

```
![](pca_3.jpeg)

We also have the output from a PCA after LD thinning. Lets compare the PC1 v PC2 plot


```r
pca_pruned<-read.table("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.eigenvec")
pca_LDpruned_named<-pca_pruned[,-1]
names(pca_LDpruned_named)<-c("sample",paste0("PC",seq(1:20)))

eigenvals_pruned<-read.table("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.eigenval",col.names="var")
eigenvals_pruned$frac<-round((eigenvals_pruned$var / sum(eigenvals_pruned$var)*100), digits=2)


pc_1v2_pruned <- pca_LDpruned_named %>%
  mutate(sample=gsub("Chinook.p","",sample)) %>%
  separate(col=sample,into=c("Pop","Ind")) %>%
  mutate(Pop=as.factor(as.numeric(Pop))) %>%
  ggplot(aes(PC2,PC1, color=Pop)) +
  geom_point() +
  ylab(paste0("PC1 (", eigenvals_pruned$frac[1], "%)")) +
  xlab(paste0("PC2 (", eigenvals_pruned$frac[2], "%)")) +
  theme_bw() +
  ggtitle("LDpruned")

#plot the pruned and unpruned together
ggarrange(pc_1v2, pc_1v2_pruned,nrow=2, common.legend=T)

```
![](pca_4.jpeg)

What happens when we infer patterns of structure without LD pruning?

**What conclusions would you draw about patterns of population structure in our data after today?**

______________________

Plotting challenge 2

-   For a single PC1 vs PC2 plot, add individual sample labels to each point on the plot.

*Hint*:
  * Try the ggrepel package
  {: .spoiler}



### 3. Genome-Wide Association Study


Last topic we used bcftools to filter snps, and plink & admixture to investigate population structure. We will plot these results to get a better visual sense of this data, but not quite yet! Before we move to R, we are going to run some analyses to infer potential patterns of selection and adaptation in our dataset. 

Remember, we want to test hypothesis about how our populations of Chinook salmon might have adapted to different temperature environments. We could do this a bunch of different ways, but lets try two. The approach we'll take is to see if we can identify loci that show particularly extreme differentiation between environments of varying temperatures. We'll use the Fst statistic to find these loci that diverge in allele frequency between extreme temperatures. We've already identified a polymorphic phenotype across this temperature gradient, so it might be interesting to see how well environmentally associated allelic variation lines up with alleles that associate with our phenotype. We can find these phenotypically-associated alleles with a GWA. 


## Fst ##

Since our PCA and structure plots were still a bit mysterious, we also might curious about what Fst among populations shows - for example, when calculating the mean Fst across all pairwise comparisons. 

We are also interested whether this any signatures of selection driven by temperature adaptation, so we might as well do these calculations at a fine-scale - in windows across the genome - to see if we can find any potentially adaptive regions. While there are tons of programs that calclulate Fst (plink, Pixy, SNPrelate), here we're going to use VCFtools.


```bash
mkdir analysis/fst_comparisons

#For vcftools --wier-fst-pop, we need a file containing sample names, for each pop

for pop in `cut -d"." -f2 vcf/Chinook_GWAS_filtered_fixedsamps.fam | uniq`
do
cut -d" " -f1 vcf/Chinook_GWAS_filtered_fixedsamps.fam | grep -w "$pop" > ~/analysis/fst_comparisons/${pop}.samples
done

```

Our input data is all set up for VCFtools but we have to consider how we can get all pairwise comparisons efficiently.
It helps to write out some _psuedocode_.

Comparisons of interest:  
For population 1, 2:10  
For population 2, 3:10  
For population 3, 4:10  
...  
For population 9, 10  

Breaking this down helps outline two major steps 1) that **for each population of interest** (left side), we want to compare **each population that comes after it**. This is a nested loop structure.

```bash
#how can we set our comparison pop to be dependent on the intial focal pop?
#dont forget how to do basic arithmetic in bash!
pop2=`echo $((1+1))`
echo $pop2

#but we want to compare our focal pop not just to the next one, but all pops that follow. seq helps with this.

pop_comp=`seq $((1+1)) 10`
echo $pop_comp

#Now we can put this all together in a nested loop to perform all population comparisons
#testing:
for pop1 in {1..9}
do
        for pop2 in `seq $((pop1+1)) 10`
        do
        echo $pop1 $pop2 
	done
done

```

Now we can use this same logic to run all pairwise population Fst comparisons

```bash
for i in {1..9}
do
	for k in `seq $((i+1)) 10`
	do

	vcftools \
	--gzvcf vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz \
	--weir-fst-pop ~/analysis/fst_comparisons/p$i.samples \
	--weir-fst-pop ~/analysis/fst_comparisons/p$k.samples \
	--out ~/analysis/fst_comparisons/pop${i}_pop${k}_10kb \
	--fst-window-size 10000 \
	--fst-window-step 10000

	done
done

#list all the log files from the vcftools fst calc runs
ls analysis/fst_comparisons/*.log

#lets check out one
less analysis/fst_comparisons/pop1_pop10_10kb.log
#what we're interested in is the weighted fst estimate
#this accounts for differences in coverage across sites

grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log
```

### Code Break Q! ###

Make a results file of the pairwise weighted fst estimates based on the grep command above. Use pipes and basic UNIX commands like _tr_, _cut_,and _sed_ to split the output into a space seperated file with three columns: 1) pop A, 2) pop B, and 3) Fst. Save it as analysis/fst_comparisons/weighted_fst_pairwise.txt

<details>
<summary markdown="span">**Answer**
</summary>
```bash
grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log | tr ":" "\t"  | sed 's|analysis/fst_comparisons/||g' | sed 's|_10kb.log||g' | cut -d$'\t' -f1,3 | tr "_" "\t" > analysis/fst_comparisons/weighted_fst_pairwise.txt
```
</details>

If you're not sure how that worked out, break down each part of the pipe, step by step and see how it changes as a result of each new command


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
