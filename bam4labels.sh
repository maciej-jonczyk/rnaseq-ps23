# this code semi-automatically changes header in count file to match samle labels
head -n2 counts1p-popr.txt | tail -n1 > xnagl
cut -f1-6 -d"	" xnagl > x
# xetyk made manually from .ods file with labels (used for multiBamSummary)
tr '\n' '\t' < xetyk > x2
paste x x2 | sed 's/\t$//'> x3
# to check with original column names
cat x3 xnagl > x5
libreoffice --calc x5
tail -n +3 counts1p-popr.txt > x6
cat x3 x6 > counts-z-etyk
# to check if header paired with data properly
head -n2 counts-z-etyk > x7
libreoffice --calc x7
# directory cleaning
rm x*
# for boxplot in R. only id and counts
cut --complement -f2-6 -d"	" counts-z-etyk > counts-z-etyk4r
