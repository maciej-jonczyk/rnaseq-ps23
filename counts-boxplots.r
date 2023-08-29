# boxplots for normalized data in DEseq
# Attention - coloring is misleading, leaf and sam data are separated WITHIN day!!!
# reading data
cts <- as.matrix(read.table("counts-z-etyk4r", sep="\t", header=T, row.names="Geneid"))
dim(cts)
class(cts)
str(cts)
# coldata
coldata=read.table("samples4r", sep="\t", header=T, row.names=1)
str(coldata)
# checking synchronisation of files
all(rownames(coldata) == colnames(cts))
# read in ranges
load("/media/mj/ANTIX-LIVE/2023_04_20.counts/coldata-ranges.RDa")
library("GenomicRanges")
library("IRanges")
# constructing ranges
zakresy=IRanges(start=dozakresow$start, end=dozakresow$end, names=dozakresow$names)
x=chr$V2
names(x)=rownames(chr)
gzakresy=GRanges(seqnames=dozakresow$seqnames, ranges=zakresy, strand=dozakresow$strand, seqlengths=x)
seqinfo(gzakresy)
summary(gzakresy)
head(gzakresy)
head(coldata)
coldata$tk <- factor(coldata$tk)
coldata$ln <- factor(coldata$ln)
coldata$time <- factor(coldata$time)
coldata$day <- factor(coldata$day)
library("DESeq2")
dds=DESeqDataSetFromMatrix(countData=cts, colData=coldata, rowRanges=gzakresy, design= ~ tk+ln+time)
nrow(dds)
# quite strong filter
keep <- rowSums(counts(dds) >= 10) >= 3
length(keep)
head(keep)
dds.f <- dds[keep,]
dim(dds.f)
vsd <- vst(dds.f, blind = FALSE)
head(assay(vsd), 3)
# save only data matrix
zlicz.vsd=assay(vsd)
dim(zlicz.vsd)
head (zlicz.vsd)
# matrix to dataframe
zlicz.vsd=as.data.frame(zlicz.vsd)
head(zlicz.vsd)
# rownames to new variable
zlicz.vsd$GeneId=rownames(zlicz.vsd)
head(zlicz.vsd)
# package for data reorganization
library(reshape2)
zlicz.vsd1kol=melt(zlicz.vsd,id.vars=c("GeneId"))
# check dimensions
dim(zlicz.vsd1kol)
dim(zlicz.vsd)
31430*96
head(zlicz.vsd1kol, n=1000)
# boxplot for vsd normalized data
boxplot(value~variable, data=zlicz.vsd1kol, col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'))
# manually maximize window with plot and saving to pdf
dev.print(pdf, 'box-counts-znorm.pdf')
# boxplot for raw data in log scale
boxplot(value~variable, data=zlicz1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'))
# saving as above
dev.print(pdf, 'box-counts-rawlog.pdf')
# Plot for FPKM normalized data
# check object for gene length
head(rowRanges(dds.f))
# normalization itself
dds.f.fpkm=fpkm(dds.f, robust = TRUE)
# check object type
class(dds.f.fpkm)
head(dds.f.fpkm, n=10)
# add column with ids as previously
dds.f.fpkm=as.data.frame(dds.f.fpkm)
dds.f.fpkm$GeneId=rownames(dds.f.fpkm)
head(dds.f.fpkm)
# reorganization of data
dds.f.fpkm1kol=melt(dds.f.fpkm,id.vars=c("GeneId"))
# check data
dim(dds.f.fpkm)
head(dds.f.fpkm,n=2)
31430*96
dim(dds.f.fpkm1kol)
# change zeroes to NA
dds.f.fpkm1kol[dds.f.fpkm1kol == 0] <- NA
head(dds.f.fpkm1kol,n=100)
# boxplot and file with it
boxplot(value~variable, data=dds.f.fpkm1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'))
dev.print(pdf, 'box-counts-logfpkm.pdf')

