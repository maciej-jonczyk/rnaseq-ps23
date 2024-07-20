# Based on https://bioconductor.org/packages/release/workflows/html/rnaseqGene.html
# Uses files made in counts4r.sh
# read-in ranges
library("GenomicRanges")
library("IRanges")
# space is a field separator
# first, construct IRanges object
dozakresow=read.table("geny4ranges",sep=" ",header=T)
zakresy=IRanges(start=dozakresow$start, end=dozakresow$end, names=dozakresow$names)
# chromosome length, don't have to be sorted
chr=read.table("counts.dlchrom", sep=" ", header=F, row.names=1)
head(chr)
# vector of lengths
x=chr$V2
class(x)
# assigning chromosome names
names(x)=rownames(chr)
gzakresy=GRanges(seqnames=dozakresow$seqnames, ranges=zakresy, strand=dozakresow$strand, seqlengths=x)
seqinfo(gzakresy)

# read-in counts
cts <- as.matrix(read.table("counts4r", sep="\t", header=T, row.names="Geneid"))
# read-in sample data
coldata=read.table("samples4r", sep="\t", header=T, row.names=1)
# check order of data
all(rownames(coldata) == colnames(cts))
# make factors
coldata$tk <- factor(coldata$tk)
coldata$ln <- factor(coldata$ln)
coldata$time <- factor(coldata$time)
# load package
library("DESeq2")
# create dataset, design included only for completness, different will be used
dds=DESeqDataSetFromMatrix(countData=cts, colData=coldata, rowRanges=gzakresy, design= ~ tk+ln+time)

# Prefiltering
nrow(dds)
keep <- rowSums(counts(dds)) > 1
# removing rows of the DESeqDataSet that have no counts, or only a single count across all samples.

############################################
# optionally, stronger filter
# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 3
############################################

length(keep)
head(keep)
dds.f <- dds[keep,]
dim(dds.f)

# testing two normalization methods
vsd <- vst(dds.f, blind = FALSE)
head(assay(vsd), 3)
rld <- rlog(dds.f, blind = FALSE)

# If fonts are not found -> https://askubuntu.com/a/1205053/150869
df <- bind_rows(
  as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")  
lvls <- c("vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)
library("ggplot2")
library("dplyr")
library("DESeq2")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

# In further analyses I use vsd-normalized data
# Clustering
library("pheatmap")
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$tk, vsd$ln, vsd$time, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Clustering with Poisson distance
install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds.f)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( vsd$tk, vsd$ln, vsd$time, sep = " - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
         
# PCA in DESeq2
plotPCA(vsd, intgroup = c("tk", "ln", "time"))
plotPCA(vsd, intgroup = c("ln", "time"))

# Multiple plots - useful for comparisons of coloring by factors
# Using multiplot.r function (available in this repo, source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
# File for plots
png("sam-pc34.png", width=90, height=30, units="cm", res=100)
# Saving individual plotPCA plots
x<-plotPCA(vsd.sam, intgroup = c("ln","powt"), ntop=30000, pcs=3:4)
x2<-plotPCA(vsd.sam, intgroup = c("day"), ntop=30000, pcs=3:4)
x3<-plotPCA(vsd.sam, intgroup = c("time"), ntop=30000, pcs=3:4)
# Plotting function itself
multiplot(x, x2, x3, cols=3)
dev.off()

# PCA with ggplot
pcaData=plotPCA(vsd, intgroup=c("ln", "time"), returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=ln, shape=time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.print(pdf, 'pca-counts.pdf')

# GLM-PCA
# Attention - this is on non-transformed data!
install.packages("glmpca")
library("glmpca")
gpca <- glmpca(counts(dds.f), L=2)
gpca.dat <- gpca$factors
gpca.dat$ln <- dds$ln
gpca.dat$time <- dds$time
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = ln, shape = time)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.print(pdf, 'glmpca-counts.pdf')
# Different picture but still without patterns

# Clustering, only 50 genes with greatest wariance of counts
# wg https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#gene-clustering
# Installing required packages in bash
sudo apt install -y libpng-dev
sudo apt install libssl-dev
# package for data filtering
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("genefilter")
library("genefilter")
library("DESeq2")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("tk", "ln","time")])
pheatmap(mat, annotation_col = anno)
dev.print(pdf, 'heat-counts.pdf')

# saving raw dataset
save(dds, file="dds.RDa")
# saving coldata and two objects for constructing ranges
save(coldata, dozakresow, chr, file="coldata-ranges.RDa")
