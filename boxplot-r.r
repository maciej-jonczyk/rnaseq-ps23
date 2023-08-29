# boxplot for RAW counts by sample
# package for reshaping data
install.packages("reshape2")
library(reshape2)
zlicz=read.table("counts-z-etyk4r", sep="\t", header=T)
# checking data
str(zlicz)
dim(zlicz)
# should have 44303 rows (genes) and 97 columns (id + samples)
# reshape to have one variable and samples as groups
zlicz1kol=melt(zlicz,id.vars=c("Geneid"))
# checking data
dim(zlicz1kol)
# should have 44303*96=4253088 rows and 3 columns
head(zlicz1kol)
# now grouping column is "variable" and data column is "value"
# most basic boxplot
boxplot(value~variable, data=zlicz1kol)
# raw counts don't show anything interesting
# removal of zero counts
zlicz1kol[zlicz1kol == 0] <- NA
# checking data
head(zlicz1kol, n=1000)
# log-scale
boxplot(value~variable, data=zlicz1kol, log="y")
# this one gives nice result but the data should be coorected for library size - FPKM?
boxplot(value~variable, data=zlicz1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'))
