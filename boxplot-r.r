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
# Attention, colors are misleading - one-tissue scheme was used for data with both tissues interleaved within day

# Separate boxplots for tissues
# using dplyr package
library(dplyr)
# column with gene names
Geneid=zlicz$Geneid
# selecting columns for leaves based on start string
x=select(zlicz,starts_with("l"))
# binding with Geneid
zlicz.lisc=cbind(Geneid,x)
dim(zlicz.lisc)
# The same for sam
x=select(zlicz,starts_with("s"))
zlicz.sam=cbind(Geneid,x)
dim(zlicz.sam)

# reshaping data
zlicz.lisc1kol=melt(zlicz.lisc,id.vars=c("Geneid"))
# checking
dim(zlicz.lisc1kol)
dim(zlicz.lisc)
44303*48
# the same for sam
zlicz.sam1kol=melt(zlicz.sam,id.vars=c("Geneid"))
dim(zlicz.sam1kol)
dim(zlicz.sam)

# zeroes to NA
zlicz.lisc1kol[zlicz.lisc1kol == 0] <- NA
# boxplot with tilted labels, according to
https://stackoverflow.com/a/61676627/1040763
# plot without label
bp=boxplot(value~variable, data=zlicz.lisc1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'), xaxt="n", xlab=NULL)
# note that both xlabel and xtitle was removed
# manually maximize plot window
# construct label
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] + 0.5, bp$names, srt = 45, xpd = TRUE)
# value added to par depends on vertical scale
dev.print(pdf, 'surowe-lisc.pdf')

# the same for sam
zlicz.sam1kol[zlicz.sam1kol == 0] <- NA
head(zlicz.sam1kol,10)
bp=boxplot(value~variable, data=zlicz.sam1kol, log="y", col=c('#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#cfd8dc', '#607d8b', '#263238', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#81d4fa', '#03a9f4', '#01579b', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#fff59d', '#fdd835', '#f57f17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17', '#ccff90', '#76ff03', '#64dd17'), xaxt="n", xlab=NULL)
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] + 0.5, bp$names, srt = 45, xpd = TRUE)
dev.print(pdf, 'surowe-sam.pdf')
