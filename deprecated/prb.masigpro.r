# Installation of required "mclust" package
# sudo apt install --no-install-recommends r-base-dev

install.packages("mclust")

# maSigPro and edgeR installation, according to bioconductor

# https://ucdavis-bioinformatics-training.github.io/2019-March-Bioinformatics-Prerequisites/thursday/linear_models.html

#####################################################################################################
# The first step is the count normalization. According to published literature TMM and cpm are used #
#####################################################################################################

# references
# https://support.bioconductor.org/p/69433/
# https://academic.oup.com/bioinformatics/article/30/18/2598/2475510?login=true
# https://www.researchgate.net/publication/341829510_Analysis_of_RNAseq_Time_Course_Data_with_maSigPro

cts <- as.matrix(read.table("../expr-anal/counts4deseq-sample", sep="\t", header=T, row.names="Geneid")) # load counts
library(edgeR)
cts.e=DGEList(counts=cts) # make an edgeR object

keep=rowSums(cpm(cts.e)>1)>=3 # filtering, minimum 1 count in at least 3 libraries (according to edgeR manual)
table(keep)
cts.e.f=cts.e[keep, keep.lib.sizes=F] # cleaned dataset

# TMM first, BTW. it is requred in edgeR analysis
cts.tmm <- calcNormFactors(cts.e.f)
# cpm, needed for maSigPro and visualisation
cts.cpm <- cpm(cts.tmm, normalized.lib.size=TRUE)

# MDS plots with colouring
plotMDS(cts.cpm, col=rep(1:2, each=12)) # here we have 2 lines, each with 12 samples

# calculation of theta, according to https://support.bioconductor.org/p/105249/
# design is taken into account, so experiment data is needed
# make model.matrix
coldata=read.table("../expr-anal/samples-rob4deseq", sep="\t", header=T, row.names=1) # the same what used in DESeq2 in other analysis
# design according to p 47 https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
design=model.matrix(~ linia + wiek + linia:wiek, coldata) # dispersion difference for model with and without interaction is negligible and close to default (10)
estimateDisp(cts.e.f, design)
# manual division, theta = 1 / common.dispersion
mytheta=1/0.1098658
theta
# 9.102014 for model with interaction
# Without interaction, only for comparison
design=model.matrix(~ linia + wiek, coldata)
estimateDisp(cts.e.f, design)
1/0.1103785
# 9.059735

#################################
# maSIgPro time-course analysis #
#################################

library(maSigPro)
# edesign for maSigPro, based on coldata (experiment data)
coldata
Time=coldata[,2]
Replicate=rep(1:8, each=3)
Control=rep(0:1, each=12)
Mutant=rep(1:0, each=12)
edesign=cbind(Time, Replicate, Control, Mutant)
edesign

# check edesign rows and counts' columns
all(rownames(edesign) == colnames(cts.cpm))
mydesign <- make.design.matrix(edesign, degree = 3) # because here are 4 time points
# minimum number of observations (cpms) in group
# 2 groups? Acording to p 3 maSigPro userguide and NBdata (example RNA-seq analysis in userguide)
mymin.obs=(3+1)*2+1 # (degree+1)xGroups+1. According to p.vector documentation
myp.vec=p.vector(cts.cpm, mydesign, Q=0.05, MT.adjust="BH", counts=T, theta=mytheta, min.obs=mymin.obs) # for RNA-seq counts=T must be given
mytfit=T.fit(myp.vec, step.method="backward", nvar.correction=T) # nvar.correction - correction for number of variables in the model. Here FDR is not appropriate because p-values comes in and out 

# check influential genes
influ=mytfit$influ.info
influ2=t(influ) # for exporting outside R and inspection
write.table(influ2, file="influ.txt", sep="\t", quote=F, row.names=T) # gees with NAa, for removal
cts.cpm.noinflu=subset(cts.cpm, !rownames(cts.cpm) %in% colnames(mytfit$influ.info)) # cleande data

#########################################################
# Optional theta check after influential genes removal  #
# firstly, count matrix must be retrieved
x=cts.e.f$counts
class(x)
x2=subset(x, !rownames(x) %in% colnames(mytfit$influ.info))
# make edgeR object with influential genes removed
cts.e.f.noinflu=DGEList(counts=x2)
# theta
estimateDisp(cts.e.f.noinflu, design)
1/0.1085462
# 9.212667
# difference at decimal point                           #
#########################################################

# Analysis after influential genes' removal
myp.vec=p.vector(cts.cpm.noinflu, mydesign, Q=0.05, MT.adjust="BH", counts=T, theta=mytheta, min.obs=mymin.obs)
mytfit=T.fit(myp.vec, step.method="backward", nvar.correction=T)
sigs=get.siggenes(mytfit, rsq=0.6, vars="groups")
# number of significant genes
sigs$sig.genes$Control$g
sigs$sig.genes$MutantvsControl$g
# less number of significant than after T.fit
# number of significant genes depends on rsq (r squared) used in get.siggenes
# high number of significant for Control vs flat profile

# clustering and plots
# load package for Mclust method
library(mclust)
# Clustering with saving to object
# see.genes overwtites cluster plots so I use modified function
source("see.gen2plot.r") # read function
seegen.res=seegen2pl(sigs$sig.genes$MutantvsControl, show.fit = T, dis=mydesign$dis, cluster.method="Mclust", k.mclust=T, cluster.data = 1)

##############################
# Plots for genes and groups #
##############################

# Plots for genes from selected cluster
# subset with ids for choosen cluster
x1=subset(seegen.res$cut, seegen.res$cut=="1")
# retrieval from dataset
genes1=cts.cpm.noinflu[names(x1),]
# plots
PlotProfiles(genes1, repvect=edesign[,2], cond=rownames(edesign))
PlotGroups(genes1, edesign = edesign, show.fit = T, dis = mydesign$dis, groups.vector = mydesign$groups.vector)

# Plots for another cluster
x2=subset(seegen.res$cut, seegen.res$cut=="2")
genes2=cts.cpm.noinflu[names(x2),]
x11() # new window for plot
PlotProfiles(genes2, repvect=edesign[,2], cond=rownames(edesign))
PlotGroups(genes2, edesign = edesign, show.fit = T, dis = mydesign$dis, groups.vector = mydesign$groups.vector)

# Validity of plots checked with see.genes result

# Plot for one gene
rownames(genes1[1:10,]) # quick look at significant genes
# PlotGroups plots median for genes so it is sensible for one gene or one cluster
Zm00001eb011690=cts.cpm.noinflu[rownames(cts.cpm.noinflu)=="Zm00001eb011690",]
PlotGroups(Zm00001eb011690, edesign = edesign, show.fit = T, dis = mydesign$dis, groups.vector = mydesign$groups.vector)
PlotProfiles(Zm00001eb011690, repvect=edesign[,2], cond=rownames(edesign))

# Potentially search for other clustering software (mfuzz?)

# https://btep.ccr.cancer.gov/docs/data-visualization-with-r/Lesson6_V2/

# Selecting MutantvsControl with rsq=0.8
x=get.siggenes(mytfit, rsq=0.8, vars="groups")
x2=subset(x$sig.genes$MutantvsControl$sig.profiles, !rownames(x$sig.genes$MutantvsControl$sig.profiles) %in% rownames(x$sig.genes$Control$sig.profiles))
dim(x2)
# [1] 22 24

# Saving all significant with rsq=0.8
write.table(x$sig.genes$MutantvsControl$sig.profiles, file="mvsc08", quote=F)
write.table(x$sig.genes$Control$sig.profiles, file="ctrl08", quote=F)


