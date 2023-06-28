# Based on https://bioconductor.org/packages/release/workflows/html/rnaseqGene.html
# read-in ranges
library("GenomicRanges")
library("IRanges")
# separtor to spacja
# najpierw obiekt IRanges
dozakresow=read.table("geny4ranges",sep=" ",header=T)
zakresy=IRanges(start=dozakresow$start, end=dozakresow$end, names=dozakresow$names)
# wczytanie długości chromosomów, nie muszą być posortowane
chr=read.table("counts.dlchrom", sep=" ", header=F, row.names=1)
head(chr)
x=chr$V2
class(x)
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
# spróbować z
# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 3
############################################

length(keep)
head(keep)
dds.f <- dds[keep,]
dim(dds.f)

# dwie metody normalizacji
vsd <- vst(dds.f, blind = FALSE)
head(assay(vsd), 3)
rld <- rlog(dds.f, blind = FALSE)

# Porównanie dwóch metod normalizacji (stosowane tylko do eksploracji)
# Jesli marudzi o braku czcionek -> https://askubuntu.com/a/1205053/150869
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

# Klastrowanie
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

# Klastrowanie prób z odległością Poisson
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
         
# PCA wbudowane w DESeq2
plotPCA(vsd, intgroup = c("tk", "ln", "time"))
plotPCA(vsd, intgroup = c("ln", "time"))

# PCA z użyciem ggplot
pcaData=plotPCA(vsd, intgroup=c("ln", "time"), returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=ln, shape=time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.print(pdf, 'pca-counts.pdf')

# GLM-PCA
# Uwaga - to robi się na danych niestransformowanych
install.packages("glmpca")
library("glmpca")
gpca <- glmpca(counts(dds.f), L=2)
gpca.dat <- gpca$factors
gpca.dat$ln <- dds$ln
gpca.dat$time <- dds$time
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = ln, shape = time)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.print(pdf, 'glmpca-counts.pdf')
# Inny obraz ale nie widać nowych wzorów

# Klastrowanie, tylko 50 genów o największej wariancji counts
# wg https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#gene-clustering
## W bashu doinstalowanie pakietu
sudo apt install -y libpng-dev
sudo apt install libssl-dev
# pakiet do filtrowania danych
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

# zapis surowego zbioru
save(dds, file="dds.RDa")
# zapis coldata i dwoch obiektow do ranges
save(coldata, dozakresow, chr, file="coldata-ranges.RDa")
