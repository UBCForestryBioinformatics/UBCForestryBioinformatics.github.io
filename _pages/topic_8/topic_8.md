---
title: "Topic 8 - SNP Filtering and Analysis"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---

### Accompanying material
[Lecture Slides](/pages/topic_8/topic_8.pdf)

In this tutorial we're going to use SNPs called with GATK to analyse patterns of population structure in the Chinook genome and conduct a GWAS. The specific data that you apply to analyses is obviously important (garbage in, garbage out). We're going to play around with some SNP filters to test their effects on downstream analyses with our data. 

# 1. Filtering SNPs for analysis

As we talked about in the lecture, the filters you choose to apply to your data can have a big impact on downstream analysis. Last topic we called variants across the two chromosomes. If you look at the VCF, you'll notice there are a lot of sites only genotyped in a small subset of the samples. This can happen with lower overall read depth (in this case this is whole genome sequencing at ~8X depth), but can be due to other factors like divergence between sample and reference. We also have indels, and SNPs with more than two alleles. Many population genomic analyses have been developed on the assumption that we are looking at two and only two alleles at a site. So first lets filter the VCF to a smaller set of usable sites.

#### NOTE:

* Yesterday we made VCFs for several individuals. For the cool analyses we want to run today, we are going to need more samples. We have a large VCF file containing many more individuals, you can copy it to from `/mnt/data/vcf/Chinook_GWAS.vcf.gz` to `~/vcf`

```bash
mkdir ~/vcf
cp /mnt/data/vcf/Chinook_GWAS.vcf.gz ~/vcf/

## Now build an index for your VCF file to help downstream software navigate the VCF efficiently
tabix ~/vcf/Chinook_GWAS.vcf.gz
```

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

### Coding challenge

* How many sites remain in the filtered VCF? How many were removed? How many on each chromosome? Don't forget about `grep -v` ad `wc -l`!

* Can you think of any additional filters that may be good to apply to your data? Try editing the BCFtools command to add additional filters that you think would be appropriate.

* What are the downsides of being too stringent with filtering? What are the downsides of being too inclusive? 

### Convert VCF to Plink format

Usually, the first analysis researchers will do with population genomic data such as we have in our VCF is to characterise any structure in the data. In particular, whether there is population structure. 

We're going to do that too using a couple of popular pieces of software, but before that we have to convert our VCF to the specialized formats those packages require. We can do that with Plink - Plink is a large set of tools for manipulating genetic data, running GWAS and calculating various stats. It's geared towards human data, so sometimes you have to force it to work with non-human data. For example, it assumes you have human chromosomes (eg 23 of them) and will complain if it doesn't see them.


```bash
cd ~/
mkdir analysis

# In our pipeline we accidentally appended extra characters to the beginning of sample names
# you can check sample names by greping for "#CHROM", which is the first string of the sample header line
zgrep "#CHROM" vcf/Chinook_GWAS_filtered.vcf.gz
# we could use a specialty software like bedtools reheader to fix this, but lets just use basic bash commands

zcat vcf/Chinook_GWAS_filtered.vcf.gz | sed 's/-e Chinook/Chinook/g' | bgzip > vcf/Chinook_GWAS_filtered_fixedsamps.vcf.gz

#the key command here is the sed 's/find/replace/g' , zcat is uncompression to standard out and bgzip is recompressing

plink=/mnt/software/plink


# Downstream software we want to run assumes numeric chr names
# We'll just convert those here (e.g. chr_1 > 1)

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

```


This produces files with the suffix .nosex, .log, .fam, .bim, .bed. 

### LD Pruning

When examining patterns of population structure (i.e. admixture/pca) its good practice to filter SNPs on the basis of linkage disequilibrium (LD). LD is a statistical association between alleles where particular configurations of SNPs are more common than expected by chance. LD arises due to restricted recombination among sites and is a fascinating topic in its own right, but for the sake of our analyses it is simply a statistical nuisance. Sites exhibiting high LD are statistically non-independant, so ideally we would remove them from analyses that assume independence. If you can't filter for LD, subsetting sites also helps (e.g. selecting every 10th site). This has the added benefit of keeping run times down!

To prune for LD, we'll ask ```plink``` to slide across the genome (10 snps at a time),and in windows of 100snps, calculate LD between each snp, removing those with LD (measured using the r2 statstic) > 0.5

```bash

#LD pruning analysis
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

#allow extra chromosome is another way we force plink to work with non-human data (i.e. allow chromosomes that don't have the typical human chr names)

#the above produces two files, .in (to include) and .out (to exclude)
#now lets actually extracts the snps that remain after LD pruning

#great! we have our files in the necessary format, now lets run some population structure analyses
```

**Questions:**
* How many SNPs were filtered out by the LD pruning?
* Can you think of any analyses where LD pruning would be disadvantageous?


# 2. Analysis of Population Structure - ADMIXTURE

With our newly pruned SNP data, we are going to run our data through ```ADMIXTURE```, a wildly popular piece of software to characterise patterns of population structure in genomic data. The aim of ```ADMIXTURE``` is to assign a proportion of an individual's genome to a set of distinct ancestral populations.

If you have not heard of ```ADMIXTURE```, do not worry. It is very likely that you have seen the output of the graphs ```ADMIXTURE``` produces some time in the past

When running ADMIXTURE we need to specify the number (K) of ancestral populations that our samples are derived from (don't worry if this sounds unintuitive - it is!). Determining the appropriate "K" in ADMIXTURE is possible with some statistical analyses. We are not going to run those here, we'll just run ```K=2```.

So, let's go ahead and throw our data into the software...

```bash
admixture=/mnt/software/admixture
#run admixture
$admixture --cv vcf/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.bed 2 |\
tee analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.out; \

#NOTE: "tee" takes the output of a command and saves it to a file, while
# also printing letting it print to the screen. So we can watch the progress while also
# saving that output.
```

With 100 samples and ~31000 SNPs across two chromosomes it should finish in less than a minute. We only ran it for one value of K (2) but we could also test different K values and select the best K value. Common practice is to run from 1 to number of populations (10). We'll skip this for the tutorial, but feel free to tinker with it at your leisure.

Let's just tidy up a little:

```bash
#now move the admixture output files to the analysis folder
mv Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*.P Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned*.Q analysis/

```


Let's look at how ```ADMIXTURE``` assigned our individuals to their ancestral groups. To do this, let's look at the .Q file, which shows group assignment for each sample. The Q file doesn't include sample names so we can put those together using "paste. One way to assess if there is meaningful population structure is to check whether each population grouping has individuals of major ancestry.

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

We can clearly see that individuals from population 1 and population 10 (furthest apart in sampling space) belong to different groupings.


We're going to plot these results, but we'll have to leave the command line to do that. If you want to do that analysis, make sure to download the output files from ADMIXTURE to your local machine.


# 3. Analysis of Population Structure - PCA

Another, extremely common way to quantify population structure among samples is with Prinical Components Analysis (PCA). Lets run a *PCA* on our data. A benefit of this approach is that it does not rely on a population genetic model like ADMIXTURE.

We can run a PCA on our data with just a single line of code in plink.

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
```

# 4. Analysis of Population Structure - Fst

Yet another standard method for characterising population structure among samples is through Fst (just don't ask me what it stands for!). 

Fst provides a quantitative measure of genetic differentiation among populations. We can use it to determine which populations are most genetically similar and dissimilar. 

Fortunately, there is a wide variety of packages that can compute Fst from VCF files. 

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

## Now grab out the weighted pairwise Fst for each comparison we made above
grep "Weir and Cockerham weighted Fst estimate:" analysis/fst_comparisons/*.log | tr ":" "\t"  | sed 's|analysis/fst_comparisons/||g' | sed 's|_10kb.log||g' | cut -d$'\t' -f1,3 | tr "_" "\t" > analysis/fst_comparisons/weighted_fst_pairwise.txt

#If you're not sure how that worked out, break down each part of the pipe, step by step and see how it changes as a result of each new command

```


___________________

# 5. Visualising Patterns of Population Structure in R 


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
### ADMIXTURE Result

We calculated cluster assignment for each sample at K =2 . Lets try loading up a single output file. We'll also need to keep track of the sample names.

``` r
samplelist <- read.table("analysis/Chinook_GWAS_fiiltered_fixedsamps_LDpruned.fam", header=F)
class(samplelist) #its useful to keep checking how your data is being treated - here its a DF

#all of these columns except the first aren't useful to us - all we need is the order of sample names. 
samplelist <- samplelist[,1] #this is called *INDEXING*
class(samplelist) #after subset, its just a vector of "character" sample names, so its no longer a DF

#read in admixture results
data <- read_delim("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr_LDpruned.2.Q",
  col_names = paste0("Q",seq(1:2)),
  delim=" ")
```

For this last command, we're using the _read\_delim_ function which requires column names. Since the data file doesn't have any column names, we have to specify them and we're using a combination of paste and seq to produce "Q1", "Q2". We could have hard coded it c("Q1","Q2") but this way it works for an arbitrary number of columns just by changing the second value in seq().

Currently our data are in a wide data format, but the tools we want to use in R prefer a long format. In a long format, each row has a single data point instead of multiple. The tidyverse tools are set up to prefer long data ((and there are other reasons)[https://sejdemyr.github.io/r-tutorials/basics/wide-and-long/#a-case-for-long-data]) so lets do that.

```r
 #converts it from wide to long format.
 data %>% gather(Q, value, -sample,-k) -> data
```

Now to plotting. We first try plotting our K=2 result.
```r
data %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) +
  geom_bar(stat="identity",position="stack")
```

![Image](/pages/topic_8_pics/structure_1.png)

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

![Image](/pages/topic_8_pics/structure_2.png)

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


![Image](/pages/topic_8_pics/structure_3.png)

### A cballenge for another time:
* Can you run ADMIXTURE on those data and make plots comparing the effects of different choices of K?


### Plotting PCA results

Lets take a look at our PCA results and compare with the ADMIXTURE graphs.

```r
pca_output<-read.table("analysis/Chinook_GWAS_filtered_fixedsamps_numericChr.eigenvec")
head(pca_output)
```

### A cballenge for another time:
--------------------

Using some of the tricks we learned above, we'll remove the redundant sample column, and then name the remaining columns (call the first "sample", the second "PC1", the third "PC2", all the way to the last column "PC20"! Save the output to "pca_named"

```R
pca_named<-pca_output[,-1]
	names(pca_named)<-c("sample",paste0("PC",seq(1:20)))
  {: .spoiler}
```

Now that you have your named columns, plot your PCA!

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


Right away we can see that while there is _some_ clusering of individuals by population, that there is nearly as much clustering of populations with each other. We might not have expected this from looking only at the structure plot. Lets visualize some other PCs of genotype variation. 

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

### Compare with LD pruned analysis...

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


