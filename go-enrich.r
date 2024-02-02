# Workflow for GO category enrichment
# goseq package is used

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
# Files "genesN" contains 
