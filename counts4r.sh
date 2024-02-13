# This file explains construction of files used to make DESeqDataSet in DESeq2in R
# before each session
export LC_ALL=C
# file after featureCounts - SEPARATE FOR EXPERIMENT
tail -n +2 counts-230420.txt > x
head -n1 x | tr '\t' '\n' | cat -n
# TAB is a field separator
cut --complement -f2-6 -d"     " x > x2


# file with sample identification - SEPARATE FOR EXPERIMENT
less proby-labels
# cut two parts of id separated by dot, to make sorted sampling time variable
# TAB is a separator
cut -f5 -d"    " proby-labels | cut -f1-2 -d"." > x
# MANUAL addition of numbers and export to TAB-separated csv
libreoffice --calc x
# merging numbers and id
tr '\t' '.' < x.csv > x3
# adding column name with id and transposition
cat <(echo Geneid) x3 | tr '\n' '\t' > x4
# MANUAL deleting TAB from end of line and moving to new line
# removal of original header from counts
tail -n +2 x2 > x5
# adding new header
cat x4 x5 > counts4r

# file with experiment description - SEPARATE FOR EXPERIMENT 
# adding modyfied column with samples
paste x3 proby-labels > x
# sampling time
cut -f3 -d"." x3 > x2
# sticking together
paste x x2 > x3
# header congruent with columns
echo id litera tk ln dl-time dl-nazwa time | tr ' ' '\t' > x2
cat x2 x3 > samples4r

############ I don't do that #############
# ranges, NOT from count result, it may have multiple ranges per gene
tail -n +2 counts-230420.txt | cut -f1-6 -d"  " > table4ranges
##########################################

# Code and data below are THE SAME FOR EXPERIMENTS FOR GENOME IN A GIVEN VERSION
# Order in counts congruent with this - they used the same GTF

# ranges from GTF
# ranges without ID, it may be done without length calculation - either way it will be not congruent with width from IRanges
awk -v FS="\t" -v OFS=" " '$3=="gene"{print $1,$4,$5,$7,$5-$4}'  ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf > x
# ID
awk -v FS="\t" -v OFS=" " '$3=="gene"'  ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf | cut -f2 -d'"' > x2
paste -d" " x2 x > x3
# header
echo names seqnames start end strand length > x4
cat x4 x3 > geny4ranges

# chromosome lengths - for all
zgrep -F "sequence-region"  ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3.gz | cut -f2,4 -d"     " > x
tr -s " " < x | cut -f2,4 -d" " > x2
sed '/^[0-9]/s/^/chr/' x2 > dlchrom
# selecting only those present in counts
cut -f2 -d"   " counts-230420.txt | tail -n +3 | cut -f1 -d";" | sort -u > counts.chrom
grep -Fwf counts.chrom dlchrom > counts.dlchrom

# File import to DeseqDataSet described in xplorative-plots.r and counts-boxplots.r
