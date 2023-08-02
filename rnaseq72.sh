# Analiza pozostałych terminów z pierwszego powtórzenia. Sekwencjonowane były razem.
# Wnioski z analizy fastqc surowych danych:
# 1. G3 R1 i R2 - mają o połowę mniej mniej unikalnych odczytów niż reszta z liści
# 2. A4 R1, dziwny pik w GC przy 65%
# 3. A4 R1, bardzo dużo over-represented sequences
# 4. A4 R1, bardzo dziwny wynik "Per base sequence content" (dobrze widoczny w wyniku fastqc). W tym przypadku również R2 dał dziwny wynik
# Odpowiedź z CeNT
# Jeśli chodzi o próbkę A4, w niej zachowały się dimery adapterów, ale ponieważ nie było ich dużo, żeby nie wprowadzać dodatkowego różnicowania
# próbek nie odczyszczaliśmy ich dodatkowo. Wg obliczeń powinny stanowić ok. 3%, ale widzę, że namnożyły się podczas klastrowania wydajniej. 
# Dosekwencjonujemy dla tej próbki jeszcze ok. 5 MR odczytów, powinno to wyrównać ilość danych, a po trimmingu adapterów statystyki powinny 
# wrócić do normy.
# Próbka G3 faktycznie odstaje. W mojej ocenie trudno stwierdzić czy jest to fenotypowy czy techniczny efekt, dlatego jak najbardziej powtórzymy 
# tę próbkę. Sprawdzimy jeszcze raz RNA przed wykonaniem biblioteki, żeby ocenić czy nie stało się coś z materiałem wejściowym.

# Trymowanie, komputer "feniks"
# 12 lipca 23
# trymowanie wszystkiego, również prób do powtorki: A4 i G3
# poziom 2023_07_10.trim
for i in A1 A3 A4 B1 B3 B4 C1 C3 C4 D1 D3 D4 E1 E3 E4 F1 F3 F4 G1 G3 G4 H1 H3 H4 I1 I3 I4 J1 J3 J4 K1 K3 K4 L1 L3 L4 M1 M3 M4 N1 N3 N4 O1 O3 O4 P1 P3 P4 R1 R3 R4 S1 S3 S4 T1 T3 T4 U1 U3 U4 W1 W3 W4 X1 X3 X4 Y1 Y3 Y4 Z1 Z3 Z4 ; do java -jar /home/mj/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 ../sekw-surowe/rna/ps23/2023_07_10/raw_data/${i}_R1_001.fastq.gz ../sekw-surowe/rna/ps23/2023_07_10/raw_data/${i}_R2_001.fastq.gz -baseout ${i}.fastq.gz ILLUMINACLIP:/home/mj/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:true MINLEN:40; done

for i in A1 A3 A4 B1 B3 B4 C1 C3 C4 D1 D3 D4 E1 E3 E4 F1 F3 F4; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_07_10.fastqc; done &
for i in G1 G3 G4 H1 H3 H4 I1 I3 I4 J1 J3 J4 K1 K3 K4 L1 L3 L4; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_07_10.fastqc; done &
for i in M1 M3 M4 N1 N3 N4 O1 O3 O4 P1 P3 P4 R1 R3 R4 S1 S3 S4; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_07_10.fastqc; done &
for i in T1 T3 T4 U1 U3 U4 W1 W3 W4 X1 X3 X4 Y1 Y3 Y4 Z1 Z3 Z4; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_07_10.fastqc; done &

shutdown -h +10

# to było źle zrobione, fastqc działało w tle a shutdown razem z tym i wyłączyło 10 min po rozpoczęciu
# Więc fastqc musiałem zrobić jeszcze raz

# Instalacja multiqc, w Trisquel musi być przez Conda
conda create --name mqc multiqc
conda activate mqc
# wywołanie z opcją zapewniającą dynamiczne wykresy przy dużej liczbie prób
multiqc --interactive .
conda deactivate
