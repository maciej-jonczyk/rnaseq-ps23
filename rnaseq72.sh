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
# Trzeba by opóźnić ostatnie fastqc i wywołać je nie w tle, wtedy shutdown musiałoby poczekać.

# Instalacja multiqc, w Trisquel musi być przez Conda
conda create --name mqc multiqc
conda activate mqc
# wywołanie z opcją zapewniającą dynamiczne wykresy przy dużej liczbie prób
multiqc --interactive .
conda deactivate

# Po przejrzeniu wyniku multiqc
# G3 R1 i R2 - nadal o połowę mniej unikalnych odczytów niż reszta z liści
# A4 R1 - zniknął pik w wykresie %GC
# A4 R1 - over-represented sequences porównywalne z resztą prób
# A4 - "Per base sequence content" porównywalne do innych

# Mapowanie tych 72 prób, wszystkich, po uzupełnieniu danych przez CeNT dorobię co trzeba.
# Włączam na dell, i tak zajmie to prawie 4 doby
# skrypt map-72.sh
# Mapowanie, poziom 2023_07_10.map
# dane sekw. są na dysku 3Tb podłączonym przez unitek
# dane genomu (indeks i NAMv5 są na dysku w dell)
# Pętla obejmuje razem mapowanie i bamqc - w ten sposób łatwiej zobaczyć czy każda część procedury zakończy się pomyślnie i ile będzie trwała
ulimit -n 10000
for i in A B C D E F G H I J K L M N O P R S T U W X Y Z; do \
for j in 1 3 4; do \
~/bin/STAR_2.7.10b/Linux_x86_64/STAR --runThreadN 8 --readFilesIn ../2023_07_10.trim/${i}${j}_1P.fastq.gz ../2023_07_10.trim/${i}${j}_2P.fastq.gz --genomeDir ~/star-index --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${i}${j} --readFilesCommand zcat
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam --threads 8 -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done

shutdown -h +10

# Na dell przewidywana długość trwania to prawie 4 doby. Na jedno mapowanie: ok. 1 h 10 min + ok 3 min na bamqc
# Pomysły na przyspieszenie, nzal. od kompa.
# 1. W STAR użycie opcji genomeLoad z wartością LoadAndKeep     ... load genome into shared and keep it in memory after run
# Ładowanie trwa tylko 43 s, więc razem dla mapowania 72 PE dałoby 52 min mniej
# 2. W bamqc można wskazać liczbę procesorów ale i tak używa jednego - wywoływać po mapowaniu w tylu kopiach ile jest procesorów.
# Dla 72 plików bam zamiast 216 min powinno dać 27 min

# Program się zawiesił - zgubione połączenie z dyskiem
# Dalsze mapowanie od H1, używam opcji ładowania genomu i trymania go w pamięci, wydzielone bamqc w osobną pętlę na wypadek jakby interferowalo z tym.
# Mapowanie, poziom 2023_07_10.map
ulimit -n 10000
for i in H I J K L M N O P R S T U W X Y Z; do \
for j in 1 3 4; do \
~/bin/STAR_2.7.10b/Linux_x86_64/STAR --limitBAMsortRAM 31619156 --genomeLoad LoadAndKeep --runThreadN 8 --readFilesIn ../2023_07_10.trim/${i}${j}_1P.fastq.gz ../2023_07_10.trim/${i}${j}_2P.fastq.gz --genomeDir ~/star-index --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${i}${j} --readFilesCommand zcat; done; done
for i in H I J K L M N O P R S T U W X Y Z; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam --threads 8 -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done

shutdown -h +10

# Ostatecznie bez przechowywania genomu w pamięci - jest jej zbyt mało
# Mapowanie, poziom 2023_07_10.map
ulimit -n 10000
for i in H I J K L M N O P R S T U W X Y Z; do \
for j in 1 3 4; do \
~/bin/STAR_2.7.10b/Linux_x86_64/STAR --runThreadN 8 --readFilesIn ../2023_07_10.trim/${i}${j}_1P.fastq.gz ../2023_07_10.trim/${i}${j}_2P.fastq.gz --genomeDir ~/star-index --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${i}${j} --readFilesCommand zcat; done; done

# bez --limitBAMsortRAM 30000000000 --genomeLoad LoadAndKeep
# The 3Tb disk has bad-sectors so I have to repeat mapping for G4 and P3 samples, as above

# New QC analysis, thea last command have run in foreground - in this way shutdown must wait for QC completion
for i in H I; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in J K; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in L M; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in N O; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in P R; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in S T; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in U W; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done &
for i in X Y Z; do \
for j in 1 3 4; do \
~/bin/BamQC/bin/bamqc ${i}${j}Aligned.sortedByCoord.out.bam -f ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_07_10.bamqc/; done; done
shutdown -h +10

# Visual inspection of bamqc results - they look okay
# Further commands run on other disk - Seagate...
# Quality filtering, sorting and indexing,2023_07_10.map level
for i in A B C D E F G H I J K L M N O P R S T U W X Y Z; do \
for j in 1 3 4; do \
samtools view -@8 -bq10 ${i}${j}Aligned.sortedByCoord.out.bam -o q10-srt-idx/${i}${j}q10.bam && samtools sort -@8 -o q10-srt-idx/${i}${j}q10srt.bam q10-srt-idx/${i}${j}q10.bam && cd q10-srt-idx && samtools index -@8 ${i}${j}q10srt.bam && rm ${i}${j}q10.bam && echo zrobione ${i}${j} && cd .. ; done; done
# In the same script multiBamSummary for each new relication was run
cd q10-srt-idx
multiBamSummary bins --bamfiles A1q10srt.bam B1q10srt.bam C1q10srt.bam D1q10srt.bam E1q10srt.bam F1q10srt.bam G1q10srt.bam H1q10srt.bam I1q10srt.bam J1q10srt.bam K1q10srt.bam L1q10srt.bam M1q10srt.bam N1q10srt.bam O1q10srt.bam P1q10srt.bam R1q10srt.bam S1q10srt.bam T1q10srt.bam U1q10srt.bam W1q10srt.bam X1q10srt.bam Y1q10srt.bam Z1q10srt.bam --labels l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 -p 8 -o ../../2023_07_10.pca/macierz1.npz
multiBamSummary bins --bamfiles A3q10srt.bam B3q10srt.bam C3q10srt.bam D3q10srt.bam E3q10srt.bam F3q10srt.bam G3q10srt.bam H3q10srt.bam I3q10srt.bam J3q10srt.bam K3q10srt.bam L3q10srt.bam M3q10srt.bam N3q10srt.bam O3q10srt.bam P3q10srt.bam R3q10srt.bam S3q10srt.bam T3q10srt.bam U3q10srt.bam W3q10srt.bam X3q10srt.bam Y3q10srt.bam Z3q10srt.bam --labels l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 -p 8 -o ../../2023_07_10.pca/macierz3.npz
multiBamSummary bins --bamfiles A4q10srt.bam B4q10srt.bam C4q10srt.bam D4q10srt.bam E4q10srt.bam F4q10srt.bam G4q10srt.bam H4q10srt.bam I4q10srt.bam J4q10srt.bam K4q10srt.bam L4q10srt.bam M4q10srt.bam N4q10srt.bam O4q10srt.bam P4q10srt.bam R4q10srt.bam S4q10srt.bam T4q10srt.bam U4q10srt.bam W4q10srt.bam X4q10srt.bam Y4q10srt.bam Z4q10srt.bam --labels l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 -p 8 -o ../../2023_07_10.pca/macierz4.npz
shutdown -h +10

# PCA plots, this time on "feniks" so in conda, 2023_07_10.pca level
conda activate deeptools
plotPCA -in macierz1.npz -l l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 --colors '#f06292' '#f06292' '#f06292' '#4dd0e1' '#4dd0e1' '#4dd0e1' '#aed581' '#aed581' '#aed581' '#ffe082' '#ffe082' '#ffe082' '#c2185b' '#c2185b' '#c2185b' '#0097a7' '#0097a7' '#0097a7' '#689f38' '#689f38' '#689f38' '#ffa000' '#ffa000' '#ffa000'  --markers 'v' 'o' 'x' --transpose -o pca1.svg
plotPCA -in macierz3.npz -l l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 --colors '#f06292' '#f06292' '#f06292' '#4dd0e1' '#4dd0e1' '#4dd0e1' '#aed581' '#aed581' '#aed581' '#ffe082' '#ffe082' '#ffe082' '#c2185b' '#c2185b' '#c2185b' '#0097a7' '#0097a7' '#0097a7' '#689f38' '#689f38' '#689f38' '#ffa000' '#ffa000' '#ffa000'  --markers 'v' 'o' 'x' --transpose -o pca3.svg
plotPCA -in macierz4.npz -l l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 --colors '#f06292' '#f06292' '#f06292' '#4dd0e1' '#4dd0e1' '#4dd0e1' '#aed581' '#aed581' '#aed581' '#ffe082' '#ffe082' '#ffe082' '#c2185b' '#c2185b' '#c2185b' '#0097a7' '#0097a7' '#0097a7' '#689f38' '#689f38' '#689f38' '#ffa000' '#ffa000' '#ffa000'  --markers 'v' 'o' 'x' --transpose -o pca4.svg
conda deactivate

# As in the 2nd 24 hours tissue separates dataset most.
# In the 3rd 24 hours the leaf samples from night cluster together, but the % variability is really low
# In the 4th 24 hours the sample l.s3.10 clearly outlies. It is L	leaf	s311	10-17 so it has L4 id.

# Clustering
conda activate deeptools
plotCorrelation --corData macierz1.npz -c pearson -p heatmap -o pear1.svg
plotCorrelation --corData macierz3.npz -c pearson -p heatmap -o pear3.svg
plotCorrelation --corData macierz4.npz -c pearson -p heatmap -o pear4.svg
conda deactivate
# Don't show anything interesting

# In the meantime at Dell multiBamSummary for all samples is running. After that fragment counting for all samples.
# For sample N highly sequenced original file is used as it don't disturb the analysis.
# 2023_07_10.map/q10-srt-idx/ level
multiBamSummary bins --bamfiles A1q10srt.bam B1q10srt.bam C1q10srt.bam D1q10srt.bam E1q10srt.bam F1q10srt.bam G1q10srt.bam H1q10srt.bam I1q10srt.bam J1q10srt.bam K1q10srt.bam L1q10srt.bam M1q10srt.bam N1q10srt.bam O1q10srt.bam P1q10srt.bam R1q10srt.bam S1q10srt.bam T1q10srt.bam U1q10srt.bam W1q10srt.bam X1q10srt.bam Y1q10srt.bam Z1q10srt.bam ../../2023_04_20.map/q10-srt-idx/Aq10srt.bam ../../2023_04_20.map/q10-srt-idx/Bq10srt.bam ../../2023_04_20.map/q10-srt-idx/Cq10srt.bam ../../2023_04_20.map/q10-srt-idx/Dq10srt.bam ../../2023_04_20.map/q10-srt-idx/Eq10srt.bam ../../2023_04_20.map/q10-srt-idx/Fq10srt.bam ../../2023_04_20.map/q10-srt-idx/Gq10srt.bam ../../2023_04_20.map/q10-srt-idx/Hq10srt.bam ../../2023_04_20.map/q10-srt-idx/Iq10srt.bam ../../2023_04_20.map/q10-srt-idx/Jq10srt.bam ../../2023_04_20.map/q10-srt-idx/Kq10srt.bam ../../2023_04_20.map/q10-srt-idx/Lq10srt.bam ../../2023_04_20.map/q10-srt-idx/Mq10srt.bam ../../2023_04_20.map/q10-srt-idx/Nq10srt.bam ../../2023_04_20.map/q10-srt-idx/Oq10srt.bam ../../2023_04_20.map/q10-srt-idx/Pq10srt.bam ../../2023_04_20.map/q10-srt-idx/Rq10srt.bam ../../2023_04_20.map/q10-srt-idx/Sq10srt.bam ../../2023_04_20.map/q10-srt-idx/Tq10srt.bam ../../2023_04_20.map/q10-srt-idx/Uq10srt.bam ../../2023_04_20.map/q10-srt-idx/Wq10srt.bam ../../2023_04_20.map/q10-srt-idx/Xq10srt.bam ../../2023_04_20.map/q10-srt-idx/Yq10srt.bam ../../2023_04_20.map/q10-srt-idx/Zq10srt.bam A3q10srt.bam B3q10srt.bam C3q10srt.bam D3q10srt.bam E3q10srt.bam F3q10srt.bam G3q10srt.bam H3q10srt.bam I3q10srt.bam J3q10srt.bam K3q10srt.bam L3q10srt.bam M3q10srt.bam N3q10srt.bam O3q10srt.bam P3q10srt.bam R3q10srt.bam S3q10srt.bam T3q10srt.bam U3q10srt.bam W3q10srt.bam X3q10srt.bam Y3q10srt.bam Z3q10srt.bam A4q10srt.bam B4q10srt.bam C4q10srt.bam D4q10srt.bam E4q10srt.bam F4q10srt.bam G4q10srt.bam H4q10srt.bam I4q10srt.bam J4q10srt.bam K4q10srt.bam L4q10srt.bam M4q10srt.bam N4q10srt.bam O4q10srt.bam P4q10srt.bam R4q10srt.bam S4q10srt.bam T4q10srt.bam U4q10srt.bam W4q10srt.bam X4q10srt.bam Y4q10srt.bam Z4q10srt.bam --labels l.a5.18d1 l.a5.02d1 l.a5.10d1 l.s0.18d1 l.s0.02d1 l.s0.10d1 l.s8.18d1 l.s8.02d1 l.s8.10d1 l.s3.18d1 l.s3.02d1 l.s3.10d1 s.a5.18d1 s.a5.02d1 s.a5.10d1 s.s0.18d1 s.s0.02d1 s.s0.10d1 s.s8.18d1 s.s8.02d1 s.s8.10d1 s.s3.18d1 s.s3.02d1 s.s3.10d1 l.a5.18d2 l.a5.02d2 l.a5.10d2 l.s0.18d2 l.s0.02d2 l.s0.10d2 l.s8.18d2 l.s8.02d2 l.s8.10d2 l.s3.18d2 l.s3.02d2 l.s3.10d2 s.a5.18d2 s.a5.02d2 s.a5.10d2 s.s0.18d2 s.s0.02d2 s.s0.10d2 s.s8.18d2 s.s8.02d2 s.s8.10d2 s.s3.18d2 s.s3.02d2 s.s3.10d2 l.a5.18d3 l.a5.02d3 l.a5.10d3 l.s0.18d3 l.s0.02d3 l.s0.10d3 l.s8.18d3 l.s8.02d3 l.s8.10d3 l.s3.18d3 l.s3.02d3 l.s3.10d3 s.a5.18d3 s.a5.02d3 s.a5.10d3 s.s0.18d3 s.s0.02d3 s.s0.10d3 s.s8.18d3 s.s8.02d3 s.s8.10d3 s.s3.18d3 s.s3.02d3 s.s3.10d3 l.a5.18d4 l.a5.02d4 l.a5.10d4 l.s0.18d4 l.s0.02d4 l.s0.10d4 l.s8.18d4 l.s8.02d4 l.s8.10d4 l.s3.18d4 l.s3.02d4 l.s3.10d4 s.a5.18d4 s.a5.02d4 s.a5.10d4 s.s0.18d4 s.s0.02d4 s.s0.10d4 s.s8.18d4 s.s8.02d4 s.s8.10d4 s.s3.18d4 s.s3.02d4 s.s3.10d4 -p 8 -o ../../2023_07_10.pca/macierz1powt.npz
~/bin/subread-2.0.6-Linux-x86_64/bin/featureCounts -C  -B  -Q 10  -F GTF  -s 2  -T 8  --ignoreDup  -p  --countReadPairs  -t exon  -g gene_id  -a ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../counts1p.txt  A1q10srt.bam B1q10srt.bam C1q10srt.bam D1q10srt.bam E1q10srt.bam F1q10srt.bam G1q10srt.bam H1q10srt.bam I1q10srt.bam J1q10srt.bam K1q10srt.bam L1q10srt.bam M1q10srt.bam N1q10srt.bam O1q10srt.bam P1q10srt.bam R1q10srt.bam S1q10srt.bam T1q10srt.bam U1q10srt.bam W1q10srt.bam X1q10srt.bam Y1q10srt.bam Z1q10srt.bam ../../2023_04_20.map/q10-srt-idx/Aq10srt.bam ../../2023_04_20.map/q10-srt-idx/Bq10srt.bam ../../2023_04_20.map/q10-srt-idx/Cq10srt.bam ../../2023_04_20.map/q10-srt-idx/Dq10srt.bam ../../2023_04_20.map/q10-srt-idx/Eq10srt.bam ../../2023_04_20.map/q10-srt-idx/Fq10srt.bam ../../2023_04_20.map/q10-srt-idx/Gq10srt.bam ../../2023_04_20.map/q10-srt-idx/Hq10srt.bam ../../2023_04_20.map/q10-srt-idx/Iq10srt.bam ../../2023_04_20.map/q10-srt-idx/Jq10srt.bam ../../2023_04_20.map/q10-srt-idx/Kq10srt.bam ../../2023_04_20.map/q10-srt-idx/Lq10srt.bam ../../2023_04_20.map/q10-srt-idx/Mq10srt.bam ../../2023_04_20.map/q10-srt-idx/Nq10srt.bam ../../2023_04_20.map/q10-srt-idx/Oq10srt.bam ../../2023_04_20.map/q10-srt-idx/Pq10srt.bam ../../2023_04_20.map/q10-srt-idx/Rq10srt.bam ../../2023_04_20.map/q10-srt-idx/Sq10srt.bam ../../2023_04_20.map/q10-srt-idx/Tq10srt.bam ../../2023_04_20.map/q10-srt-idx/Uq10srt.bam ../../2023_04_20.map/q10-srt-idx/Wq10srt.bam ../../2023_04_20.map/q10-srt-idx/Xq10srt.bam ../../2023_04_20.map/q10-srt-idx/Yq10srt.bam ../../2023_04_20.map/q10-srt-idx/Zq10srt.bam A3q10srt.bam B3q10srt.bam C3q10srt.bam D3q10srt.bam E3q10srt.bam F3q10srt.bam G3q10srt.bam H3q10srt.bam I3q10srt.bam J3q10srt.bam K3q10srt.bam L3q10srt.bam M3q10srt.bam N3q10srt.bam O3q10srt.bam P3q10srt.bam R3q10srt.bam S3q10srt.bam T3q10srt.bam U3q10srt.bam W3q10srt.bam X3q10srt.bam Y3q10srt.bam Z3q10srt.bam A4q10srt.bam B4q10srt.bam C4q10srt.bam D4q10srt.bam E4q10srt.bam F4q10srt.bam G4q10srt.bam H4q10srt.bam I4q10srt.bam J4q10srt.bam K4q10srt.bam L4q10srt.bam M4q10srt.bam N4q10srt.bam O4q10srt.bam P4q10srt.bam R4q10srt.bam S4q10srt.bam T4q10srt.bam U4q10srt.bam W4q10srt.bam X4q10srt.bam Y4q10srt.bam Z4q10srt.bam  --verbose
shutdown -h +10

# Moving count results to new directory
mkdir 2023_07_10.counts
mv counts1p.txt* 2023_07_10.counts/
# statistics
multiqc .
# conclusions, even though few samples looked like outliers at the sount level differences are not visible.

# Tissue strongly divide dataset what hides inter-tissue variability.
# To cope with this separate analyses are made by tissue.
# poziom 2023_07_10.map/q10-srt-idx/
multiBamSummary bins --bamfiles A1q10srt.bam B1q10srt.bam C1q10srt.bam D1q10srt.bam E1q10srt.bam F1q10srt.bam G1q10srt.bam H1q10srt.bam I1q10srt.bam J1q10srt.bam K1q10srt.bam L1q10srt.bam ../../2023_04_20.map/q10-srt-idx/Aq10srt.bam ../../2023_04_20.map/q10-srt-idx/Bq10srt.bam ../../2023_04_20.map/q10-srt-idx/Cq10srt.bam ../../2023_04_20.map/q10-srt-idx/Dq10srt.bam ../../2023_04_20.map/q10-srt-idx/Eq10srt.bam ../../2023_04_20.map/q10-srt-idx/Fq10srt.bam ../../2023_04_20.map/q10-srt-idx/Gq10srt.bam ../../2023_04_20.map/q10-srt-idx/Hq10srt.bam ../../2023_04_20.map/q10-srt-idx/Iq10srt.bam ../../2023_04_20.map/q10-srt-idx/Jq10srt.bam ../../2023_04_20.map/q10-srt-idx/Kq10srt.bam ../../2023_04_20.map/q10-srt-idx/Lq10srt.bam A3q10srt.bam B3q10srt.bam C3q10srt.bam D3q10srt.bam E3q10srt.bam F3q10srt.bam G3q10srt.bam H3q10srt.bam I3q10srt.bam J3q10srt.bam K3q10srt.bam L3q10srt.bam A4q10srt.bam B4q10srt.bam C4q10srt.bam D4q10srt.bam E4q10srt.bam F4q10srt.bam G4q10srt.bam H4q10srt.bam I4q10srt.bam J4q10srt.bam K4q10srt.bam L4q10srt.bam --labels a5.18d1 a5.02d1 a5.10d1 s0.18d1 s0.02d1 s0.10d1 s8.18d1 s8.02d1 s8.10d1 s3.18d1 s3.02d1 s3.10d1 a5.18d2 a5.02d2 a5.10d2 s0.18d2 s0.02d2 s0.10d2 s8.18d2 s8.02d2 s8.10d2 s3.18d2 s3.02d2 s3.10d2 a5.18d3 a5.02d3 a5.10d3 s0.18d3 s0.02d3 s0.10d3 s8.18d3 s8.02d3 s8.10d3 s3.18d3 s3.02d3 s3.10d3 a5.18d4 a5.02d4 a5.10d4 s0.18d4 s0.02d4 s0.10d4 s8.18d4 s8.02d4 s8.10d4 s3.18d4 s3.02d4 s3.10d4 -p 8 -o ../../macierz1powt-lisc.npz
multiBamSummary bins --bamfiles M1q10srt.bam N1q10srt.bam O1q10srt.bam P1q10srt.bam R1q10srt.bam S1q10srt.bam T1q10srt.bam U1q10srt.bam W1q10srt.bam X1q10srt.bam Y1q10srt.bam Z1q10srt.bam ../../2023_04_20.map/q10-srt-idx/Mq10srt.bam ../../2023_04_20.map/q10-srt-idx/Nq10srt.bam ../../2023_04_20.map/q10-srt-idx/Oq10srt.bam ../../2023_04_20.map/q10-srt-idx/Pq10srt.bam ../../2023_04_20.map/q10-srt-idx/Rq10srt.bam ../../2023_04_20.map/q10-srt-idx/Sq10srt.bam ../../2023_04_20.map/q10-srt-idx/Tq10srt.bam ../../2023_04_20.map/q10-srt-idx/Uq10srt.bam ../../2023_04_20.map/q10-srt-idx/Wq10srt.bam ../../2023_04_20.map/q10-srt-idx/Xq10srt.bam ../../2023_04_20.map/q10-srt-idx/Yq10srt.bam ../../2023_04_20.map/q10-srt-idx/Zq10srt.bam M3q10srt.bam N3q10srt.bam O3q10srt.bam P3q10srt.bam R3q10srt.bam S3q10srt.bam T3q10srt.bam U3q10srt.bam W3q10srt.bam X3q10srt.bam Y3q10srt.bam Z3q10srt.bam M4q10srt.bam N4q10srt.bam O4q10srt.bam P4q10srt.bam R4q10srt.bam S4q10srt.bam T4q10srt.bam U4q10srt.bam W4q10srt.bam X4q10srt.bam Y4q10srt.bam Z4q10srt.bam --labels a5.18d1 a5.02d1 a5.10d1 s0.18d1 s0.02d1 s0.10d1 s8.18d1 s8.02d1 s8.10d1 s3.18d1 s3.02d1 s3.10d1 a5.18d2 a5.02d2 a5.10d2 s0.18d2 s0.02d2 s0.10d2 s8.18d2 s8.02d2 s8.10d2 s3.18d2 s3.02d2 s3.10d2 a5.18d3 a5.02d3 a5.10d3 s0.18d3 s0.02d3 s0.10d3 s8.18d3 s8.02d3 s8.10d3 s3.18d3 s3.02d3 s3.10d3 a5.18d4 a5.02d4 a5.10d4 s0.18d4 s0.02d4 s0.10d4 s8.18d4 s8.02d4 s8.10d4 s3.18d4 s3.02d4 s3.10d4 -p 8 -o ../../macierz1powt-sam.npz
shutdown -h +10
# It run ca 5 h on dell

# PCA, tweaking with colors
plotPCA -in macierz1powt-lisc.npz -l a5.18d1      a5.02d1 a5.10d1 s0.18d1 s0.02d1 s0.10d1 s8.18d1 s8.02d1 s8.10d1 s3.18d1 s3.02d1 s3.10d1      a5.18d2 a5.02d2 a5.10d2 s0.18d2 s0.02d2 s0.10d2 s8.18d2 s8.02d2 s8.10d2 s3.18d2 s3.02d2 s3.10d2 a5.18d3 a5.02d3 a5.10d3 s0.18d3 s0.02d3 s0.10d3 s8.18d3 s8.02d3 s8.10d3 s3.18d3 s3.02d3      s3.10d3 a5.18d4 a5.02d4 a5.10d4 s0.18d4 s0.02d4 s0.10d4 s8.18d4 s8.02d4 s8.10d4 s3.18d4 s3.02d4 s3.10d4 --colors  '#cfd8dc'     '#607d8b'       '#263238'       '#cfd8dc'       '#607d8b'    '#263238'       '#cfd8dc'       '#607d8b'       '#263238'       '#cfd8dc'       '#607d8b'       '#263238'       '#81d4fa'       '#03a9f4'       '#01579b'       '#81d4fa'       '#03a9f4'    '#01579b'       '#81d4fa'       '#03a9f4'       '#01579b'       '#81d4fa'       '#03a9f4'       '#01579b'       '#fff59d'       '#fdd835'       '#f57f17'       '#fff59d'       '#fdd835'    '#f57f17'       '#fff59d'       '#fdd835'       '#f57f17'       '#fff59d'       '#fdd835'       '#f57f17'       '#ccff90'       '#76ff03'       '#64dd17'       '#ccff90'       '#76ff03'    '#64dd17'       '#ccff90'       '#76ff03'       '#64dd17'       '#ccff90'       '#76ff03'       '#64dd17' --markers 'v' 'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'  'd'     'd'     'v'     'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'     'd'     'd'     'v'     'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'  'd'     'd'     'v'     'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'     'd'     'd' --transpose --plotHeight 20 --plotWidth 20 -o pca1p-lisc.svg
plotPCA -in macierz1powt-sam.npz -l a5.18d1      a5.02d1 a5.10d1 s0.18d1 s0.02d1 s0.10d1 s8.18d1 s8.02d1 s8.10d1 s3.18d1 s3.02d1 s3.10d1      a5.18d2 a5.02d2 a5.10d2 s0.18d2 s0.02d2 s0.10d2 s8.18d2 s8.02d2 s8.10d2 s3.18d2 s3.02d2 s3.10d2 a5.18d3 a5.02d3 a5.10d3 s0.18d3 s0.02d3 s0.10d3 s8.18d3 s8.02d3 s8.10d3 s3.18d3 s3.02d3      s3.10d3 a5.18d4 a5.02d4 a5.10d4 s0.18d4 s0.02d4 s0.10d4 s8.18d4 s8.02d4 s8.10d4 s3.18d4 s3.02d4 s3.10d4 --colors  '#cfd8dc'     '#607d8b'       '#263238'       '#cfd8dc'       '#607d8b'    '#263238'       '#cfd8dc'       '#607d8b'       '#263238'       '#cfd8dc'       '#607d8b'       '#263238'       '#81d4fa'       '#03a9f4'       '#01579b'       '#81d4fa'       '#03a9f4'    '#01579b'       '#81d4fa'       '#03a9f4'       '#01579b'       '#81d4fa'       '#03a9f4'       '#01579b'       '#fff59d'       '#fdd835'       '#f57f17'       '#fff59d'       '#fdd835'    '#f57f17'       '#fff59d'       '#fdd835'       '#f57f17'       '#fff59d'       '#fdd835'       '#f57f17'       '#ccff90'       '#76ff03'       '#64dd17'       '#ccff90'       '#76ff03'    '#64dd17'       '#ccff90'       '#76ff03'       '#64dd17'       '#ccff90'       '#76ff03'       '#64dd17' --markers 'v' 'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'  'd'     'd'     'v'     'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'     'd'     'd'     'v'     'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'  'd'     'd'     'v'     'v'     'v'     's'     's'     's'     'o'     'o'     'o'     'd'     'd'     'd' --transpose --plotHeight 20 --plotWidth 20 -o pca1p-sam.svg

# in case of leaves two 'cold' periods are separated from two 'regrowth' ones.
# for sam this is not visible
# Also it seems that 's3' line is somewhat separated by 2nd PC

# 24.08.23 - CeNT dosłał A4 uzupełnione o 20% odczytów (sklejone stare fastq z 10.07 i nowe 20%) i nowe sekwencjonowanie G3 (G3_S9)
# Analiza fastqc od nich OK
# Trymowanie i Fastqc
for i in A4 G3_S9 ; do java -jar /home/mj/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 ../../sekw-surowe/rna/ps23/2023_08_24/raw_data/${i}_R1_001.fastq.gz ../../sekw-surowe/rna/ps23/2023_08_24/raw_data/${i}_R2_001.fastq.gz -baseout ${i}.fastq.gz ILLUMINACLIP:/home/mj/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:true MINLEN:40; done
~/bin/FastQC/fastqc A4_1P.fastq.gz --outdir=../2023_08_24.fastqc &
~/bin/FastQC/fastqc A4_2P.fastq.gz --outdir=../2023_08_24.fastqc &
~/bin/FastQC/fastqc G3_S9_1P.fastq.gz --outdir=../2023_08_24.fastqc &
~/bin/FastQC/fastqc G3_S9_2P.fastq.gz --outdir=../2023_08_24.fastqc &

# Mapowanie
for i in A4 G3_S9 ; do ~/bin/STAR_2.7.10b/Linux_x86_64/STAR --runThreadN 24 --readFilesIn ../2023_08_24.trim/${i}_1P.fastq.gz ../2023_08_24.trim/${i}_2P.fastq.gz --genomeDir ../../star-index --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${i} --readFilesCommand zcat; done

# bamqc
~/bin/BamQC-master/bin/bamqc G3_S9Aligned.sortedByCoord.out.bam --threads 24 -f ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_08_24.bamqc/ &
~/bin/BamQC-master/bin/bamqc A4Aligned.sortedByCoord.out.bam --threads 24 -f ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_08_24.bamqc/ &
# Wynik OK

# Filtrowanie i sortowanie
for i in A4 G3_S9 ; do samtools view -@24 -bq10 ${i}Aligned.sortedByCoord.out.bam -o q10-srt-idx/${i}q10.bam && samtools sort -@24 -o q10-srt-idx/${i}q10srt.bam q10-srt-idx/${i}q10.bam && cd q10-srt-idx && samtools index -@24 ${i}q10srt.bam && rm ${i}q10.bam && echo zrobione ${i} && cd .. ; done

# Zrobić multiqc na surowych, zestaw jak niżej

# Multiqc na wszystkich plikach fastqc (zamiana starego A4 na nowy, dodanie G3_S9 obok starego G3) wyglada OK.
# Analiza z poziomu nadrzędnego do katalogów
multiqc --interactive 2023_04_20.fastqc 2023_07_10.fastqc/ 2023_08_24.fastqc

# Zrobić multibamsummary osobno dla tkanek, z podmienionym A4 nowym zamiast starego, w pierwszej analizie dla liści dodać G3_S9 bez uzuwania G3.
# Usunąć G3 po sprawdzeniu wyniku z powyższego.

# Dalej Counts dla całego zbioru, będzie trzeba wymienić pliki po kolei bo nie są alfabetycznie
