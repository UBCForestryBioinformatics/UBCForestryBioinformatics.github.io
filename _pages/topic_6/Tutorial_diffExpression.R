rm(list=ls())
#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("DESeq2")
library("DESeq2")

directory <- "~/UBC/ZoologyStuff/BIOL525D/Topic_X_RNA/read_counts/"
setwd(directory)


# Set the prefix for each output file name
outputPrefix <- "RNA_tutorial"

# List the names for each of the count files
sampleFiles<- c("warm_sample_01.read_counts.c.txt",
                "warm_sample_03.read_counts.c.txt",
                "warm_sample_02.read_counts.c.txt",
                "warm_sample_04.read_counts.c.txt",
                "warm_sample_05.read_counts.c.txt",
                "warm_sample_06.read_counts.c.txt",
                "cold_sample_07.read_counts.c.txt",
                "cold_sample_08.read_counts.c.txt",
                "cold_sample_09.read_counts.c.txt",
                "cold_sample_10.read_counts.c.txt",
                "cold_sample_11.read_counts.c.txt",
                "cold_sample_12.read_counts.c.txt")

# List the sample names
sampleNames<- c("sample_01",
                "sample_02",
                "sample_03",
                "sample_04",
                "sample_05",
                "sample_06",
                "sample_07",
                "sample_08",
                "sample_09",
                "sample_10",
                "sample_11",
                "sample_12")

# A simple little function to generate
sampleCondition <- c(rep(c("warm","cold"), 6))

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments = c("warm","cold")



# Built the DESeq dataset object using the counts data
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

dds <- estimateSizeFactors(ddsHTSeq)

# Take a look at the raw counts
View(counts(ddsHTSeq))
# Take a look at the normalised counts
View(counts(dds, normalized=TRUE))

str(ddsHTSeq)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)


#guts
dds <- DESeq(ddsHTSeq)
res <- results(dds)
plot(res$log2FoldChange, -log10(res$padj))

# copied from: https://benchtobioinformatics.wordpress.com/category/dexseq/
# order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]



# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))



# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))



ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))


# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()


# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t'))

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()



# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))

#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(png, paste0(outputPrefix, "-clustering.png"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))


# set condition
condition <- treatments
scores <- data.frame(pc$x)
scores$condition <- rep(c("warm","cold"),each =6)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))

