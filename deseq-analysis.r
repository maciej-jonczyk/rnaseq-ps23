# Time-course analysis DESeq2

https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/
https://www.r-bloggers.com/2021/04/note-for-deseq2-time-course-analysis/
https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/

library("DESeq2")

# Reading counts and variables, according to
https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input

# Files prepared in bash shell
# counts
cts <- as.matrix(read.table("counts4deseq-sample", sep="\t", header=T, row.names="Geneid"))
# variables
coldata=read.table("samples-rob4deseq", sep="\t", header=T, row.names=1)
# check - critical
all(rownames(coldata) == colnames(cts))

# Coding variables as factors and rearranging (default 1st alphabetical level is set as reference)
coldata$wiek <- factor(coldata$wiek)
coldata$linia <- factor(coldata$linia, , levels = c("wt","vp5"))

# Adding grouping variable - will be used for averaging for plotting
coldata$lw=paste(coldata$linia, coldata$wiek, sep="_")
coldata$lw <- factor(coldata$lw)

# Reading range data, from files prepared in bash shell
# Description in "counts4r.sh" file, for this analysis entire "dlchrom" was used 
# (without excluding molecules not present in count file as it don't contain chromosome info)
dozakresow=read.table("geny4ranges",sep=" ",header=T)
head(dozakresow)
library("GenomicRanges")
library("IRanges")
zakresy=IRanges(start=dozakresow$start, end=dozakresow$end, names=dozakresow$names)
head(zakresy)
chr=read.table("dlchrom", sep=" ", header=F, row.names=1)
head(chr)
x=chr$V2
class(x)
names(x)=rownames(chr)
head(x)
gzakresy=GRanges(seqnames=dozakresow$seqnames, ranges=zakresy, strand=dozakresow$strand, seqlengths=x)
seqinfo(gzakresy)

# DESeq Data Set FOR THREE REPS
library("DESeq2")
dds=DESeqDataSetFromMatrix(countData=cts, colData=coldata, rowRanges=gzakresy, design= ~ linia + wiek + linia:wiek)
all(rownames(coldata) == colnames(cts))

# examining data
seqinfo(gzakresy)
summary(gzakresy)
head(gzakresy)
head(rowRanges(dds), 3)
summary(dds)
assayNames(dds)
head(counts(dds), 2)
head(assay(dds), 3)

# SUBSETTING TO 2 REPLICATIONS
# I don't know how to do it using row/col names
coldata2p=coldata[-c(3,6,9,12,15,18,21,24), ]
rownames(coldata2p)
head(coldata, 3)
cts2powt=cts[,-c(3,6,9,12,15,18,21,24)]
colnames(cts2powt)
head(cts2powt, 3)

# TWO REPLICATIONS, control of input data
all(rownames(coldata2p) == colnames(cts2powt))

# TWO REPLICATIONS, Construction of DESeqDataSet using counts, variables and ranges (chromosome/contig and transcript lengths)
dds2=DESeqDataSetFromMatrix(countData=cts2powt, colData=coldata2p, rowRanges=gzakresy, design= ~ linia + wiek + linia:wiek)

# examining data
assayNames(dds2)
head(assay(dds2), 3)
head(colData(dds2), 3)
nrow(colData(dds2))
design(dds2)

# light filtering, removal of genes without counts or with only one count
nrow(dds2)
keep <- rowSums(counts(dds2)) > 1
dds2.f = dds2[keep,]
nrow(dds2.f)

## More general code for filtering
keep <- rowSums( counts(dds) >= X ) >= Y
dds <- dds[keep,]
## This requires genes to have Y or more samples with counts of X or more.
# It therefore filters out genes that have less than Y samples with counts of X or more.

# According to author, filtering at this stage only makes analysis faster. During test and adjusted p-val calculation independent filtering is done.
# So genes with low counts are filtered-out.
# https://support.bioconductor.org/p/65256/#65260

# Likelihood Ratio Test (LRT) analysis in DESeq2, only different shape of profiles gives significant result
# According to
https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
dds2.f=DESeq(dds2.f, test="LRT", reduced = ~ linia + wiek)
summary(dds2.f)
summary(dds2)
# object with result statistics
res2=results(dds2.f)
# sizeFactors are computed automatically: vector assigns to each column of the count matrix a value, the size factor, such that count values 
# in the columns can be brought to a common scale by dividing by the corresponding size factor
summary(res2)
dim(res2)
head(res2)

# filtering results
res2Filt <- res2[which(res2$padj < 0.05), ]
dim(res2)
dim(res2Filt)

# order filtered results by increasing adjusted p-value
res2FiltOrd <- res2Filt[order(res2Filt$padj),]
head(res2FiltOrd, 2)
tail(res2FiltOrd, 2)
head(res2FiltOrd[order(res2FiltOrd$padj),], 4)
# this file is not used further, leaved only for plotCounts tests

# saving objects for 3 and 2 replications
save(chr, coldata2p, cts2powt, dds2, dds2.f, dozakresow, gzakresy, res2, res2Filt, res2FiltOrd, zakresy, file="prb2powt-ranges.RDa")
save(chr, coldata, cts, dds, dozakresow, gzakresy, zakresy, file="prb3powt-ranges.RDa")
rm(list=ls())

# loading saved objects for TWO replications only - test before analysis of two reps of real experiment
load("prb2powt-ranges.RDa")

#################################################################
# Required for Explorative Analysis if done before DE-analysis, in DE-analysis it is done by DESeq function
dds <- estimateSizeFactors(dds)
#################################################################

# four best results
head(res2[order(res2$padj),], 4)
# time-course plot for the best resultlibrary(ggplot2)
# THIS COMMAND (USE OF "which.min(res2FiltOrd$padj)") RETURNS WRONG GENE
prb <- plotCounts(dds2.f, which.min(res2FiltOrd$padj), intgroup = c("wiek","linia"), returnData = TRUE) # for gene with minimal p-value
prb$wiek <- as.numeric(as.character(prb$wiek))
png("min.png")
ggplot(prb, aes(x = wiek, y = count, color = linia, group = linia)) + geom_point() + stat_summary(fun=mean, geom="line") +   scale_y_log10()
dev.off()
head(res2FiltOrd, 2)

# THE BELOW COMMANDS GIVE THE SAME CORRECT RESULT
prb <- plotCounts(dds2.f, "Zm00001eb371740", intgroup = c("wiek","linia"), returnData = TRUE) # for gene with minimal p-value
prb$wiek <- as.numeric(as.character(prb$wiek))
png("recz.png")
ggplot(prb, aes(x = wiek, y = count, color = linia, group = linia)) + geom_point() + stat_summary(fun=mean, geom="line") +   scale_y_log10()
dev.off()
head(res2, 2)
prb <- plotCounts(dds2.f, which.min(res2$padj), intgroup = c("wiek","linia"), returnData = TRUE) # for gene with minimal p-value
prb$wiek <- as.numeric(as.character(prb$wiek))
png("min-res2.png")
ggplot(prb, aes(x = wiek, y = count, color = linia, group = linia)) + geom_point() + stat_summary(fun=mean, geom="line") +   scale_y_log10()
dev.off()
res2FiltOrd["Zm00001eb371740",]
res2["Zm00001eb371740",]
# IT SEEMS THAT FOR "plotCounts" ONE HAVE TO USE DIRECT RESULT FROM "results" FUNCTION, HERE "res2"

head(rowRanges(dds2), 3)
seqinfo(gzakresy)

# selecting counts fo DE genes only, according to
https://statisticsglobe.com/subset-data-frame-and-matrix-by-row-names-in-r
ids005=rownames(res2Filt)
dds2sig=dds2[rownames(dds2) %in% ids005, ]
# the result is DESeqDataSet
dim(res2)
dim(res2Filt)
class(res2Filt)

# Next analysis - fuzzy clustering in Mfuzz, "mfuzz-analysis.r"

#############################################################################################################################
## Just illustrative - test for one time-point
#res22 <- results(dds, name="liniavp5.wiek22", test="Wald")
## Najlepszy wynik
#res22[which.min(res22$padj),]

# Comparisons and interpretation, s 52
# https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf

# NOTE, below code executed on older data!
# heatmap
install.packages("pheatmap") # if not installed already
library(pheatmap)
betas <- coef(dds)
colnames(betas)
# 20 best
topGenes <- head(order(res$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE)

# plots for the best result
head(results(dds))
summary(res)
plotCounts(dds, gene="Zm00001eb403230", intgroup="vp5")
plotCounts(dds, gene="Zm00001eb403230", intgroup="linia")
plotCounts(dds, gene="Zm00001eb403230", intgroup="wiek")
plotCounts(dds, gene="Zm00001eb403230", intgroup="wiek:linia")

# volcano-plot, seems different than in the tutorial - logFC in results has nothing to do with profile differences
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# PCA
# wg https://lashlock.github.io/compbio/R_presentation.html
# variance stabilization (only for explorative analysis)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="linia")

# To check
# https://github.com/csgillespie/bmc-microarray

# Próba symulowanej analizy dla więcej niż dwóch grup
# https://stackoverflow.com/questions/6422273/how-to-randomize-or-permute-a-dataframe-rowwise-and-columnwise
# https://www.geeksforgeeks.org/how-to-create-a-matrix-with-random-values-in-r/
# https://www.c-sharpcorner.com/article/matrix-in-r-operation-on-matrix-adding-two-matrix/
