# plik po featureCounts
tail -n +2 counts-230420.txt > x
head -n1 x | tr '\t' '\n' | cat -n
cut --complement -f2-6 -d"     " x > x2


# plik z identyfikacją prób
less proby-labels
# wycięcie dwóch części id rozdzielanych kropką, żeby zrobić terminy po cyframi kolei
cut -f5 -d"    " proby-labels | cut -f1-2 -d"." > x
# ręczne dodanie cyfr i eksport do csv
libreoffice --calc x
# scalenie liczb do id
tr '\t' '.' < x.csv > x3
# dodanie nazwy kolumny z id i transpozycja
cat <(echo Geneid) x3 | tr '\n' '\t' > x4
# ręczne usunięcie TABa z końca i przejście do nowego wiersza
# usunięcie oryginalnego nagłówka
tail -n +2 x2 > x5
# dodanie nowego nagłówka
cat x4 x5 > counts4r

# plik z opisem eksperymentu
# dodanie zmodyfoikowanej kolumny z próbami
paste x3 proby-labels > x
# terminy
cut -f3 -d"." x3 > x2
# połączenie
paste x x2 > x3
# nagłówek zgodny z kolumnami
echo id litera tk ln dl-time dl-nazwa time | tr ' ' '\t' > x2
cat x2 x3 > samples4r
