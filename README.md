# rnaseq-ps23
RNA-seq pipeline, from raw reads to biological interpretaion of time-course experiment.
Whenewer possible bash command-line tools are preffered over R.

*alternatywy.r* - just some links to web pages with potentially usefu information.

## Analyses on simulated or test data
The first trials of read-processing are given in *rnaseq_pipeline.sh*

*prb.masigpro.r* - test of time-course analysis in maSigPro R package.

*see.gen2plot.r* - modification of the see.genes function, used to prevent over-writing of plots in R (used from command-line)
In the end this package will not be used as it don't allow to compare multiple groups and the contact with developers is poor.

*prb-expr-deseq.r* - test of differential gene expression analysis in DESeq2 package in R. This package will be most probably used for real analysis.

*go-enrich.r* - GO over and underrepresentation test in goseq

- [ ] other explorative analyses

## Analyses on real data
*rnaseq-ost.sh* - the clean read-processing workflow, based on *rnaseq_pipeline.sh*

*rnaseq72.sh* - the clean read-processing workflow for whole first replication of experiment

*qc-pca-all1powt* - quality control plots from BAM files (PCA and correlation)

*counts4r.sh* - construction of files used for DESeqDataSet in DESeq2 package in R

*bam4labels.sh* - commands used to change header of the counts file to match sample labels

*boxplot-r.r* - simple boxplot of raw counts

*counts-boxplots.r* - more advanced boxplots of read counts, involving transformation and FPKM calculation

*explorative-plots.r* - exploring patterns in count data from first sequencing results (PCA and clustering)

- [ ] correct data after new sequencing of "U3" sample - just (02.02.24) sent to sequencing along with 2nd replication of experiment
- [ ] full analysis of 2nd replication, plus joint PCA/clustering analysis
- [ ] update sample names in count file to include replication number
- [ ] check, if it is possible to make DE analysis on two replications (DEseq2)
