# Fuzzy clustering of time profiles using Mfuzz package
# starting point is "dds2sig" DESeqDataSet from DESeq2 analysis (file deseq-analysis.r)

sig.fpkm=fpkm(dds2sig, robust=TRUE)
head(sig.fpkm, 5)

# examining data
class(sig.fpkm)
dim(sig.fpkm)
colnames(sig.fpkm)
head(exprs[,1:5])
head(sig.fpkm[,1:5])
dim(coldata2p)
rownames(coldata2p)
summary(coldata2p)
class(coldata2p)

# Below - test simplest ExpressionSet for mfuzz test
# According to
https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf
# eset <- new("ExpressionSet",exprs=sig.fpkm)
# alternative exp.set <- ExpressionSet(assayData = as.matrix(vsd.matrix))
# summary(eset)
# this object is NOT used further

# class(eset)
# head(eset)
# head(exprs(eset))

# does obiects are S4 class?
# isS4(eset)
isS4(dds2)

# check if elements of data are in concert
all(rownames(coldata2p)==colnames(sig.fpkm))

# examining data
str(coldata2p)
names(coldata2p)
summary(coldata2p)
sapply(coldata2p, class)
coldata2p[c(1,5), c("linia", "wiek")]
coldata2p[c(1,5), ]

# variable description
metadata=data.frame(labelDescription=c("maize inbred line", "age in days", "concatenated variable"), row.names=c("linia", "wiek", "lw"))
metadata
phenoData=new("AnnotatedDataFrame", data=coldata2p, varMetadata=metadata)
phenoData
head(pData(phenoData))
phenoData[c("m22.1", "m29.2"), "linia"]
summary(phenoData)
summary(coldata2p)
phenoData
sampleNames(phenoData)
pData(phenoData)
varMetadata(phenoData)
head(pData(phenoData))

# ExpressionSet with metadata
eset2=ExpressionSet(assayData=sig.fpkm, phenoData=phenoData)

# Averaging replicates according to
https://support.bioconductor.org/p/9156659/#9156662
average= t(avereps(t(exprs(eset2)), ID = eset2$lw)) # output is a matrix

# New ExpressionSet
# phenoData is not used as it now is not congruent with data values
eset2.average=ExpressionSet(assayData=average)

# log2 transformation, make distribution more normal-like
eset2.log <- eset2.average
exprs(eset2.log) <- log(exprs(eset2.average) + 1, 2)
## There is function for it, so equivalent is
# NO, following formula requires DESeqDataSet, not ExpressionSet: ntd <- normTransform(dds)

# standardisation
library(Mfuzz)
# THIS IS OBJECT FOR CLUSTERING
eset2.std=standardise(eset2.log)
head(exprs(eset2.std))
dim(exprs(eset2.std))

# check if there are NAs
any(is.na(exprs(eset2.std)))

# Variance Stabilising Transfomation vst in DESeq2
# using VST normalisation as advised in https://support.bioconductor.org/p/9156659/#9156692
# blind=FALSE as descibed in https://www.biostars.org/p/428369/#428372
vsd <- vst(dds2, blind=FALSE)
class(vsd)
# vst done on full dataset, so next I select only DE
class(ids005)
vsd.sig=vsd[rownames(vsd) %in% ids005, ]
dim(vsd.sig)
# check
all(rownames(vsd.sig)==rownames(eset2.std))
# make ExpressionSet
library(limma)
vsd.eset <- ExpressionSet(assayData=assay(vsd.sig), phenoData=phenoData)
# Averaging
average= t(avereps(t(exprs(vsd.eset)), ID = vsd.eset$lw)) # output is a matrix
vsd.eset.ave=ExpressionSet(assayData=average)

all(colnames(vsd.eset)==colnames(eset2.std))
# standardise
# THIS IS OBJECT FOR CLUSTERING
vsd.eset.std=standardise(vsd.eset.ave)
head(exprs(vsd.eset.std))

# Clustering

# estimate "m"
m1 <- mestimate(eset2.std) # *
cl1=mfuzz(eset2.std, c=16, m=m1)
mfuzz.plot(eset2.std, cl=cl1, mfrow=c(4,4))

## Deciding on cluster parameters ##
## clustering coefficient "m" estimation as shown above *
## the function partcoef can be used to test, whether random data is clustered for a particular setting of "m"
## Determininig optimal number of clusters ##
# For now decision made according to examination of "overlap.plot" (with "m" from "mestimate") and joining adjacent clusters
O=overlap(cl12)
Ptmp=overlap.plot(cl12, over=O, thres=0.05)

## TAKEN FROM partcoeff DESCRIPTION IN MANUAL
## codes from old project ##
## ctrl2cl_std is a ExpressionSet object after standardisation
eset2.stdR=randomise(eset2.std) # randomisation
cl9R=mfuzz(eset2.stdR, c=9, m=m) # clustering using "m" estimated for real data by "mestimate"
mfuzz.plot2(eset2.stdR, cl=cl9R, mfrow=c(3,3))
# uniform clustering, but does this "m" is minimal one making such clustering?
cl9R=mfuzz(eset2.stdR, c=9, m=1.25)
mfuzz.plot2(eset2.stdR, cl=cl9R, mfrow=c(3,3))
# shows cluster structure
# Test if random data are clustered
tmpR=partcoef(eset2.stdR, crange=seq(8,13,1), mrange=seq(1.25,2,0.05))
# decision based on part. coef and part. coef for uniform partition comparison
F <- tmpR[[1]];F.n <- tmpR[[2]];F.min <- tmpR[[3]]
F > 1.01 * F.min
     m:1.25 m:1.3 m:1.35 m:1.4 m:1.45 m:1.5 m:1.55 m:1.6 m:1.65 m:1.7 m:1.75
c:8    TRUE  TRUE   TRUE  TRUE  FALSE FALSE  FALSE FALSE  FALSE FALSE  FALSE
c:9    TRUE  TRUE   TRUE  TRUE  FALSE FALSE  FALSE FALSE  FALSE FALSE  FALSE
c:10   TRUE  TRUE   TRUE  TRUE  FALSE FALSE  FALSE FALSE  FALSE FALSE  FALSE
c:11   TRUE  TRUE   TRUE  TRUE  FALSE FALSE  FALSE FALSE  FALSE FALSE  FALSE
c:12   TRUE  TRUE   TRUE  TRUE  FALSE FALSE  FALSE FALSE  FALSE FALSE  FALSE
c:13   TRUE  TRUE   TRUE  TRUE  FALSE FALSE  FALSE FALSE  FALSE FALSE  FALSE
     m:1.8 m:1.85 m:1.9 m:1.95   m:2
c:8  FALSE  FALSE FALSE  FALSE FALSE
c:9  FALSE  FALSE FALSE  FALSE FALSE
c:10 FALSE  FALSE FALSE  FALSE FALSE
c:11 FALSE  FALSE FALSE  FALSE FALSE
c:12 FALSE  FALSE FALSE  FALSE FALSE
c:13 FALSE  FALSE FALSE  FALSE FALSE
# at selected granularity level 1.45 is the lowest "m" not clustering random data
# Selection can be fine-tuned using different values of "m"
tmpR=partcoef(eset2.stdR, crange=seq(8,13,1), mrange=seq(1.33,1.53,0.05))
F <- tmpR[[1]];F.n <- tmpR[[2]];F.min <- tmpR[[3]]
F > 1.01 * F.min
     m:1.33 m:1.38 m:1.43 m:1.48 m:1.53
c:8    TRUE   TRUE  FALSE  FALSE  FALSE
c:9    TRUE   TRUE  FALSE  FALSE  FALSE
c:10   TRUE   TRUE  FALSE  FALSE  FALSE
c:11   TRUE   TRUE   TRUE  FALSE  FALSE
c:12   TRUE   TRUE   TRUE  FALSE  FALSE
c:13   TRUE   TRUE   TRUE  FALSE  FALSE
# Note - there could be slight differences between "partcoef" execution with the same settings
cl9R=mfuzz(eset2.stdR, c=9, m=1.43)
mfuzz.plot2(eset2.stdR, cl=cl9R, mfrow=c(3,3)) # give uniform clustering
# CONCLUSION - in this case "mestimate" gave similar value as "partcoef"

# Choice of cluster number
# cselection
tmp <- cselection(eset2.std,m=m,crange=seq(8,16,1),repeats=10,visu=T)
# don't give empty clusters even for very high number of clusters(40)
# Dmin
tmp <- Dmin(eset2.std,m=m,crange=seq(6,20,1),repeats=10,visu=TRUE)
# there are gaps between 7-8 and 9-10 so either 7 or 9 clusters could be ok
# 9 clusters has been selected earlier using "overlap.plot"


# Extracting core genes for all and one sample
# genes with membership>=0.5
core10vst=acore(vsd.eset.std, cl10vst, min.acore=0.5)
# core10vst is a list with separate table for each cluster
class(core10vst)
summary(core10vst)
head(core10vst[[1]])
head(core10vst[[7]]$NAME)
# extracting IDs for a given cluster
ids.cl7=rownames(core10vst[[7]])
# counts with library-size normalization
dds2.fn=counts(dds2.f, normalized=TRUE)
class(dds2$sig)
class(dds2sig)
# extracting genes from dds
dds2.fn7=dds2.fn[rownames(dds2.fn) %in% ids.cl7, ]
class(dds2.fn7)
# dds2.fn7 is a matrix NOT dds
dim(dds2.fn7)

# saving session
save.image("prb-mfuzz190224.RDa")
savehistory("dds2exprset150224.r")

# For color bar for plots "marray" package is needed, available from bioconductor
# color bar is produced by modified function, saved in colorbar.R
