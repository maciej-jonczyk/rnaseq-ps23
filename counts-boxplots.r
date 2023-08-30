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

# Sepoarate plots for tissues
# VSD
colnames(zlicz.vsd)
x=select(zlicz.vsd,starts_with("s"))
head(x,2)
x$GeneId=rownames(x)
head(x,2)
zlicz.sam.vsd=x
head(zlicz.sam.vsd,2)
x=select(zlicz.vsd,starts_with("l"))
x$GeneId=rownames(x)
head(x,2)
zlicz.lisc.vsd=x
head(zlicz.lisc.vsd,2)
zlicz.lisc.vsd1kol=melt(zlicz.lisc.vsd,id.vars=c("GeneId"))
zlicz.sam.vsd1kol=melt(zlicz.sam.vsd,id.vars=c("GeneId"))
bp=boxplot(value~variable, data=zlicz.lisc.vsd1kol, col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'), xaxt="n", xlab=NULL)
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.7, bp$names, srt = 45, xpd = TRUE)
dev.print(pdf, 'vsd-lisc.pdf')
bp=boxplot(value~variable, data=zlicz.sam.vsd1kol, col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'), xaxt="n", xlab=NULL)
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.7, bp$names, srt = 45, xpd = TRUE)
dev.print(pdf, 'vsd-sam.pdf')
# FPKM
head(dds.f.fpkm,2)
x=select(dds.f.fpkm,starts_with("l"))
head(x,2)
x$GeneId=rownames(x)
head(x,2)
zlicz.lisc.fpkm=x
x=select(dds.f.fpkm,starts_with("s"))
head(x,2)
x$GeneId=rownames(x)
head(x,2)
zlicz.sam.fpkm=x
zlicz.sam.fpkm1kol=melt(zlicz.sam.fpkm,id.vars=c("GeneId"))
zlicz.lisc.fpkm1kol=melt(zlicz.lisc.fpkm,id.vars=c("GeneId"))
zlicz.sam.fpkm1kol[zlicz.sam.fpkm1kol == 0] <- NA
zlicz.lisc.fpkm1kol[zlicz.lisc.fpkm1kol == 0] <- NA
bp=boxplot(value~variable, data=zlicz.lisc.vsd1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'), xaxt="n", xlab=NULL)
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] + 3.7, bp$names, srt = 45, xpd = TRUE)
dev.print(pdf, 'fpkm-lisc.pdf')
bp=boxplot(value~variable, data=zlicz.sam.vsd1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'), xaxt="n", xlab=NULL)
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] + 3.7, bp$names, srt = 45, xpd = TRUE)
dev.print(pdf, 'fpkm-sam.pdf')
