# Analiza time-course w DESeq2

# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/08b_time_course_analyses.html
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
# https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/

# Wczytanie danych (wczesniej zrobione w bash)
# Wg https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
# counts - tylko counts, 1. kolumna to Geneid, dalej po kolei próby z krótką czytelną nazwą
cts <- as.matrix(read.table("counts4deseq-sample", sep="\t", header=T, row.names="Geneid"))
# coldata - tylko potrzebne kolumny czyli próba i czynniki
coldata=read.table("samples-rob4deseq", sep="\t", header=T, row.names=1)

# Konieczne sprawdzenie, czy kolejność prób w 1. wierszu counts jest ta sama co w pierwszej kolumnie coldata
all(rownames(coldata) == colnames(cts))

# Zakodowanie czynnikiow jako factor i uporządkowanie (domyślnie 1. alfabetycznie ustawiony poziom jest brany jako odniesienie)
coldata$wiek <- factor(coldata$wiek)
coldata$linia <- factor(coldata$linia, , levels = c("wt","vp5"))

library("DESeq2")
# konstrukcja zbioru do analizy
dds=DESeqDataSetFromMatrix(countData=cts, colData=coldata, design= ~ linia + wiek + linia:wiek)

# Oglądanie poszczególnych składników
assayNames(dds)
head(assay(dds), 3)
head(colData(dds), 3)
colSums(assay(dds))
# suma counts
nrow(colData(dds))
design(dds)
# po DESeq również nazwy wyników
resultsNames(dds)

# filtrowanie (nie używałem na razie w analizie)
# Wg https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pre-filtering-the-dataset
# usuniete wiersze bez counts albo z tylko jedym zliczeniem
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds.f = dds[keep,]
nrow(dds.f)

# wg autora, filtrowanie na tym etapie tylko oprzyśpiesza analizę. I tak jest robione independent filtering przy wyliczaniuadj p-val
# https://support.bioconductor.org/p/65256/#65260
# Będę używał takiego filtrowania


# Analiza
# Wg https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
# Likelyhood Ratio Test - istotne tylko te o różnym kształcie profilu w czasie
dds=DESeq(dds, test="LRT", reduced = ~ linia + wiek)
res=results(dds)
# Automatycznie obliczane sizeFactors: vector assigns to each column of the count matrix a value, the size factor, such that count values 
# in the columns can be brought to a common scale by dividing by the corresponding size factor

# Wybór istotnych. Tylko wg adjp, bo w tym przypadku logFC nie ma sensu
# To eksportuje tylko tabelę z wynikami a NIE wartości (A w ogóle, to jakie wartości są potrzebne? counts, fpkm?)
resFilt <- res[which(res$padj < 0.05), ]
# Uporządkowanie względem rosnących adj p-val
resFiltOrd <- resFilt[order(resFilt$padj),]
# Eksport do csv
write.csv(resFiltOrd, file = "results.csv")

#################################################################
# Required for Explorative Analysis if done before DE-analysis, in DE-analysis it is done by DESeq function
dds <- estimateSizeFactors(dds)
#################################################################

# 4 najlepsze wyniki
head(res[order(res$padj),], 4)
# Wykres time-course dla najlepszego wyniku wg adj p-val
library(ggplot2)
prb <- plotCounts(dds, which.min(res$padj), intgroup = c("wiek","linia"), returnData = TRUE)
prb$wiek <- as.numeric(as.character(prb$wiek))
ggplot(prb, aes(x = wiek, y = count, color = linia, group = linia)) + geom_point() + stat_summary(fun=mean, geom="line") +   scale_y_log10()

# test dla pojedynczego punktu czasowego
res22 <- results(dds, name="liniavp5.wiek22", test="Wald")
# Najlepszy wynik
res22[which.min(res22$padj),]

# Porównania i ich interpretacja, s 52
# https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf

# heatmapa
install.packages("pheatmap") # if not installed already
library(pheatmap)
betas <- coef(dds)
colnames(betas)
# 20 najlepszych
topGenes <- head(order(res$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE)

# wykresy dla najlepszego
head(results(dds))
summary(res)
plotCounts(dds, gene="Zm00001eb403230", intgroup="vp5")
plotCounts(dds, gene="Zm00001eb403230", intgroup="linia")
plotCounts(dds, gene="Zm00001eb403230", intgroup="wiek")
plotCounts(dds, gene="Zm00001eb403230", intgroup="wiek:linia")

# volcano-plot, nie wygląda jak w tutorialu - podana w wynikach logFC nie odnosi się do
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# PCA
# wg https://lashlock.github.io/compbio/R_presentation.html
# stabilizacja wariancji (takie rzeczy tylko do eksploracji = wykresów)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="linia")

# Do sprawdzenia
# https://github.com/csgillespie/bmc-microarray

# Próba symulowanej analizy dla więcej niż dwóch grup
# https://stackoverflow.com/questions/6422273/how-to-randomize-or-permute-a-dataframe-rowwise-and-columnwise
# https://www.geeksforgeeks.org/how-to-create-a-matrix-with-random-values-in-r/
# https://www.c-sharpcorner.com/article/matrix-in-r-operation-on-matrix-adding-two-matrix/
