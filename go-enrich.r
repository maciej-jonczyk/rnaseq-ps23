# Workflow for GO category enrichment
# goseq package is used
https://bioconductor.org/packages/release/bioc/manuals/goseq/man/goseq.pdf
https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf

## Preparation of files ##
# The starting point is file with counts created with featureCounts in command-line.
# The required action is retrieval of 1st (gene IDs, ) and 6th (length) columns.
# Below example.
# cut -f1,6 -d"        " ../2023_08_24.counts/counts1p-popr.txt | tail -n +3 > x
# cut -f1 -d"  " x > assayed.genes
# cut -f2 -d"  " x > length.vec
# Once created files are valid for a given genome version.
# GO categories, created from gaf (gene association file) here from Maizegamer project. Should be valid for every gaf file.
# tail -n +2 3.1_B73v5.MaizeGDB.CLEANED.gaf | cut -f2,5 -d"        " > ../prb-analfun/go.cats

## Analysis itself ##
# Assayed genes, ie. all known maize genes. Reading in file with IDs one in a line, prepared in command-line.
assayed.genes=scan("assayed.genes", what="character")
# check
class(assayed.genes)
str(assayed.genes)
dim(assayed.genes)
length(assayed.genes)
[1] 44303
# including non-coding genes

# Length vector
length.vec=scan("length.vec")
# check
str(length.vec)
dim(length.vec)
length(length.vec)

# gene to GO categories
go.cats=read.table("go.cats", sep="\t", header=T)
# check
head(go.cats, 2)
str(go.cats)

# Three above files are valid for future analyses as long as genome version (assayed.genes, length.vec) or GO annotation (go.cats) don't change.

# Vector of DE genes, here I used clustring results from maSigPro.
# This package will NOT be used, but procedure will be the same.
# Files "genesN" contains genes from individual clusters detected by maSigPro. Appropriate IDs were retrieved from cts.cpm.moinflu object.
# The values here are cpm but they are not used anyways, only gene list is needed

# Vector of DE gene IDs
de.genes1=rownames(genes1)
# check
class(de.genes1)
length(de.genes1)
[1] 746
# Named vector with 1/0 coding, 1 = DE, 0 = no DE, code from documentation
gene.vector1=as.integer(assayed.genes%in%de.genes1)
names(gene.vector1)=assayed.genes
# check
gene.vector1[1:10]
assayed.genes[1:10]
de.genes1[1:10]

# Fitting the probability weighting function (pwf), saving plot to file
png("nullp1.png")
pwf=nullp(gene.vector1, bias.data=length.vec)
dev.off()

# If the plot fits badly to data -> problems with data
# according to
https://support.bioconductor.org/p/9156310/#9156340

# To recreate the plot if needed
plotPWF(pwf)

# Enrichment, both over and undrrepresentation
res1=goseq(pwf, gene2cat=go.cats, method="Wallenius", use_genes_without_cat=FALSE)
# Wallenius is default method, genes without categories are excluded
# check
head(res1)
tail(res1)

# There is an option to restrict test to only selected branches (BP, MF, CC)
# p 13 frm https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf

# CAUTION, the result is not corrected for Multiple testing
# Correction according to deseq manual
enriched1=res1$category[p.adjust(res1$over_represented_pvalue, method="BH")<.05]
# here I got zero rows as raw p-values were very high

# Proceed with uncorrected results - ONLY for workflow checking
res1.nocorr=subset(res1, over_represented_pvalue<=0.05, select=category)
res1.nocorr.vec=res1.nocorr$category

# Retrieval of categories description
library(GO.db)

for(go in res1.nocorr.vec[1:10]){
print(GOTERM[[go]])
cat("--------------------------------------\n")
}

# Exporting GO categories to external file
sink('my_list.txt')
print(my_list)
sink()
