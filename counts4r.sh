# przed każdą sesją
export LC_ALL=C
# plik po featureCounts - OSOBNY DLA DOSWIADCZENIA
tail -n +2 counts-230420.txt > x
head -n1 x | tr '\t' '\n' | cat -n
# separatorem jest TAB
cut --complement -f2-6 -d"     " x > x2


# plik z identyfikacją prób - OSOBNY DLA DOSWIADCZENIA
less proby-labels
# wycięcie dwóch części id rozdzielanych kropką, żeby zrobić terminy po cyframi kolei
# pierwotny separator to TAB
cut -f5 -d"    " proby-labels | cut -f1-2 -d"." > x
# ręczne dodanie cyfr i eksport do csv z TAB jako separatorem
libreoffice --calc x
# scalenie liczb do id
tr '\t' '.' < x.csv > x3
# dodanie nazwy kolumny z id i transpozycja
cat <(echo Geneid) x3 | tr '\n' '\t' > x4
# ręczne usunięcie TABa z końca i przejście do nowego wiersza
# usunięcie oryginalnego nagłówka z counts
tail -n +2 x2 > x5
# dodanie nowego nagłówka
cat x4 x5 > counts4r

# plik z opisem eksperymentu - OSOBNY DLA DOSWIADCZENIA 
# dodanie zmodyfoikowanej kolumny z próbami
paste x3 proby-labels > x
# terminy
cut -f3 -d"." x3 > x2
# połączenie
paste x x2 > x3
# nagłówek zgodny z kolumnami
echo id litera tk ln dl-time dl-nazwa time | tr ' ' '\t' > x2
cat x2 x3 > samples4r

############ tego nie robić #############
# zakresy, NIE z wyniku zliczeń, może zawierać wiele zakresów dla jednego genu
tail -n +2 counts-230420.txt | cut -f1-6 -d"  " > table4ranges
##########################################

# Ponizsze WSPOLNE DLA EKSPERYMENTOW DLA GENOMU W DANEJ WERSJI
# Kolejnośc w zliczeniach zgodna z tym bo one korzystały z tego GTF

# zakresy z GTF
# zakresy bez ID, można bez obliczania długości - i tak nie zgodzi się z width z IRanges
awk -v FS="\t" -v OFS=" " '$3=="gene"{print $1,$4,$5,$7,$5-$4}'  ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf > x
# ID
awk -v FS="\t" -v OFS=" " '$3=="gene"'  ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf | cut -f2 -d'"' > x2
paste -d" " x2 x > x3
# nagłówek
echo names seqnames start end strand length > x4
cat x4 x3 > geny4ranges

# długości chromosomów - wszystkich
zgrep -F "sequence-region"  ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3.gz | cut -f2,4 -d"     " > x
tr -s " " < x | cut -f2,4 -d" " > x2
sed '/^[0-9]/s/^/chr/' x2 > dlchrom
# wybór tylko tych w zliczeniach
cut -f2 -d"   " counts-230420.txt | tail -n +3 | cut -f1 -d";" | sort -u > counts.chrom
grep -Fwf counts.chrom dlchrom > counts.dlchrom
