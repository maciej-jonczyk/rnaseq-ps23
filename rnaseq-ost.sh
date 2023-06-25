# CeNT daje wyniki z fastqc
# Uwaga, ścieżki i liczba procesorów/RAM mogą się różnić - robiłem na różnych komputerach!
# Trymowanie i fastqc robione na dell!
# Trymowanie - zwykle zostaje ,,universal adapter''. Też poli-A i poli-T ale tego nie wyrzucam.
# z poziomu 2023_04_20.trim
for i in A B C D E F G H I J K L M N O P R S T U W X Y Z; do /home/nev/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 ../2023_04_20/raw_data/${i}_R1_001.fastq.gz ../2023_04_20/raw_data/${i}_R2_001.fastq.gz -baseout ${i}.fastq.gz ILLUMINACLIP:/home/nev/bin/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:true MINLEN:40; done
# fastq dla wytrymowanych. Nie ma hyperthreadingu więc podział na tyle procedur ile jest procesorów.
for i in A B C; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in D E F; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in G H I; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in J K L; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in M N O; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in P R S; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in T U W; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &
for i in X Y Z; do ~/bin/FastQC/fastqc ${i}_[12]P.fastq.gz --outdir=../2023_04_20.fastqc; done &

# Ręczne sprawdzenie wyników fastqc po trymowaniu

# Mapowanie, poziom 2023_04_20.map
ulimit -n 10000
for i in A B C D E F G H I J K L M N O P R S T U W X Y Z; do STAR --runThreadN 24 --readFilesIn ../2023_04_20.trim/${i}_1P.fastq.gz ../2023_04_20.trim/${i}_2P.fastq.gz --genomeDir ../star-index --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${i} --readFilesCommand zcat; done

# bamqc
for i in B C D E F G H I J K L M N O P R S T U W X Y Z; do ~/bin/BamQC/bin/bamqc ${i}Aligned.sortedByCoord.out.bam --threads 24 -f ../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../2023_04_20.bamqc/; done

# Ręczne sprawdzenie wyników bamqc po trymowaniu
# Ostrzeżenie w insert size spodziewane - sprawa uwzględniania intronów w analizie

# Odfiltrowanie MAQC >= 10. FeatureCounts tego nie potrzebuje ale do PCA lepiej same dobre dane wziąć.
# Dalej sortowanie (niekoniecznie potrzebne, najwyżej do IGV), usunięcie pośrednich bam-ów
# Dalej indeksowanie - Musi być zrobione przed multiBamSummary
for i in A B C D E F G H I J K L M N O P R S T U W X Y Z; do samtools view -@24 -bq10 ${i}Aligned.sortedByCoord.out.bam -o q10-srt-idx/${i}q10.bam && samtools sort -@24 -o q10-srt-idx/${i}q10srt.bam q10-srt-idx/${i}q10.bam && cd q10-srt-idx && samtools index -@24 ${i}q10srt.bam && rm ${i}q10.bam && echo zrobione ${i} && cd .. ; done

# PCA i klastrowanie, poziom 2023_04_20.map/q10-srt-idx

multiBamSummary bins --bamfiles Aq10srt.bam Bq10srt.bam Cq10srt.bam Dq10srt.bam Eq10srt.bam Fq10srt.bam Gq10srt.bam Hq10srt.bam Iq10srt.bam Jq10srt.bam Kq10srt.bam Lq10srt.bam Mq10srt.bam Nq10srt.bam Oq10srt.bam Pq10srt.bam Rq10srt.bam Sq10srt.bam Tq10srt.bam Uq10srt.bam Wq10srt.bam Xq10srt.bam Yq10srt.bam Zq10srt.bam --labels l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 -p 24 -o ../../2023_04_20.pca/macierz.npz

# PCA, poziom 2023_04_20.pca
plotPCA -in macierz.npz -l l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 --colors '#ff9e80' '#ff9e80' '#ff9e80' '#ff6e40' '#ff6e40' '#ff6e40' '#ff3d00' '#ff3d00' '#ff3d00' '#d50000' '#d50000' '#d50000' '#84ffff' '#84ffff' '#84ffff' '#40c4ff' '#40c4ff' '#40c4ff' '#00b0ff' '#00b0ff' '#00b0ff' '#2962ff' '#2962ff' '#2962ff'  --markers 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' --transpose -o prb-pca.svg

# z różnymi sybolami dla terminów
plotPCA -in macierz.npz -l l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 --colors '#ff9e80' '#ff9e80' '#ff9e80' '#ff6e40' '#ff6e40' '#ff6e40' '#ff3d00' '#ff3d00' '#ff3d00' '#d50000' '#d50000' '#d50000' '#84ffff' '#84ffff' '#84ffff' '#40c4ff' '#40c4ff' '#40c4ff' '#00b0ff' '#00b0ff' '#00b0ff' '#2962ff' '#2962ff' '#2962ff'  --markers 'v' 'o' 'x' --transpose -o prb-pca2.svg

# Pomyśleć o zmianie kolorów dla linii - większy kontrast ale podobne bo rozdział na tkanki

# Klastrowanie
plotCorrelation --corData macierz.npz -c pearson -p heatmap -o pear.svg
# można dodać --removeOutliers ale zwykle to nie pomaga

# zliczenie fragmentów
~/bin/subread-2.0.6-Linux-x86_64/bin/featureCounts -C  -B  -Q 10  -F GTF  -s 2  -T 4  --ignoreDup  -p  --countReadPairs  -t exon  -g gene_id  -a ~/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o counts-230420.txt  /media/mj/a59e16f5-75c1-41cb-98e5-c1db006c8e7c/2023_04_20.map/q10-srt-idx/*.bam  --verbose
# Statystyki zliczeń, poziom /media/mj/Seagate Backup Plus Drive/anal-sekw-ps23/2023_04_20.counts
# Rozdzielenie reverse stranded (2, to jest OK) i stranded (1)
mkdir xcounts2 xcounts1
mv *str1* xcounts1
mv *txt* xcounts2
# Iw każdym z katalogów
multiqc .

###################################################################################################################################################
###################################################################################################################################################
# Downsampling próby N
# Wydobycie liczby counts z reads forward (bo reverse mają tyle samo). Ręcznie z raportów html fastqc.
# Zapis do 2023_04_20.readnum i dodanie w calc ID literowego
# Wyliczenie śr liczby odczytów bez próby N

# ***********IE, za dużo RAMu wymaga, nawet na superkomp ******************
# Downsampling w seqtk
# instalacja dependencies
sudo apt install make
sudo apt install gcc
sudo apt install zlib1g
sudo apt install zlib1g-dev
# instalacja seqtk, w bin
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
#****************************************************************************8

# downsampling ze skryptami z BBmap
# pobranie binaries z sourceforge i rozpakowanie
# Uwaga, inne ścieżki bo robiłem na Dell
# reformat.sh marudzi co do plików PE więc najpierw re-parowanie
~/bin/bbmap/repair.sh in=../2023_04_20.trim/N_1P.fastq.gz in2=../2023_04_20.trim/N_2P.fastq.gz out=N_1P.rep.fastq.gz out2=N_2P.rep.fastq.gz
Input:                          170098586 reads                 22036125152 bases.
Result:                         170098586 reads (100.00%)       22036125152 bases (100.00%)
Pairs:                          170098586 reads (100.00%)       22036125152 bases (100.00%)
Singletons:                     0 reads (0.00%)         0 bases (0.00%)

Time:                           581.135 seconds.
Reads Processed:        170m    292.70k reads/sec
Bases Processed:      22036m    37.92m bases/sec
# kontrolne fastqc
~/bin/FastQC/fastqc *.fastq.gz --outdir=.
# downsampling
~/bin/bbmap/reformat.sh in=N_1P.rep.fastq.gz in2=N_2P.rep.fastq.gz out=N_1P.dwn.fastq.gz out2=N_2P.dwn.fastq.gz samplereadstarget=45780325 sampleseed=1 overwrite=t
Set INTERLEAVED to false
Input is being processed as paired
Input:                          170098586 reads                 22036125152 bases
Output:                         91560650 reads (53.83%)         11861138140 bases (53.83%)
                                                                                                                                                                                             
Time:                           395.597 seconds.                                                                                                                                             
Reads Processed:        170m    429.98k reads/sec                                                                                                                                            
Bases Processed:      22036m    55.70m bases/sec 
# i znów fastqc ale tylko dla downsampled
~/bin/FastQC/fastqc *dwn.fastq.gz --outdir=.
# mapowanie
~/bin/STAR_2.7.10b/Linux_x86_64/STAR --runThreadN 8 --readFilesIn N_1P.dwn.fastq.gz N_2P.dwn.fastq.gz --genomeDir /media/mj/SAMSUNG/z-toshiba/star-index/ --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Ndwn --readFilesCommand zcat
# filtrowanie jakości, sortowanie, indeksowanie
samtools view -@8 -bq10 NdwnAligned.sortedByCoord.out.bam -o Ndwnq10.bam && samtools sort -@8 -o Ndwnq10srt.bam Ndwnq10.bam && samtools index -@8 Ndwnq10srt.bam && rm Ndwnq10.bam
# analiza jkościowa bam
~/bin/BamQC/bin/bamqc Ndwnq10srt.bam --threads 8 -f /media/mj/SAMSUNG/z-toshiba/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o .
# porównanie jakości, po przeniesieniu plików fastqc dla N do pozostałych
multiqc . -o ../multiqcn/

# Macierz do PCA
# poziom /media/mj/Seagate Backup Plus Drive/anal-sekw-ps23/2023_04_20.map/q10-srt-idx
multiBamSummary bins --bamfiles Aq10srt.bam Bq10srt.bam Cq10srt.bam Dq10srt.bam Eq10srt.bam Fq10srt.bam Gq10srt.bam Hq10srt.bam Iq10srt.bam Jq10srt.bam Kq10srt.bam Lq10srt.bam Mq10srt.bam ../../Ndwn/Ndwnq10srt.bam Oq10srt.bam Pq10srt.bam Rq10srt.bam Sq10srt.bam Tq10srt.bam Uq10srt.bam Wq10srt.bam Xq10srt.bam Yq10srt.bam Zq10srt.bam --labels l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 -p 8 -o ../../2023_04_20.pcan/macierz.npz
# PCA, kolory z https://htmlcolors.com/
plotPCA -in macierz.npz -l l.a5.18 l.a5.02 l.a5.10 l.s0.18 l.s0.02 l.s0.10 l.s8.18 l.s8.02 l.s8.10 l.s3.18 l.s3.02 l.s3.10 s.a5.18 s.a5.02 s.a5.10 s.s0.18 s.s0.02 s.s0.10 s.s8.18 s.s8.02 s.s8.10 s.s3.18 s.s3.02 s.s3.10 --colors '#f06292' '#f06292' '#f06292' '#4dd0e1' '#4dd0e1' '#4dd0e1' '#aed581' '#aed581' '#aed581' '#ffe082' '#ffe082' '#ffe082' '#c2185b' '#c2185b' '#c2185b' '#0097a7' '#0097a7' '#0097a7' '#689f38' '#689f38' '#689f38' '#ffa000' '#ffa000' '#ffa000'  --markers 'v' 'o' 'x' --transpose -o prb-pca-n.svg
# korelacja Pearsona
plotCorrelation --corData macierz.npz -c pearson -p heatmap -o pear.svg
# Zdecydowanie rozdziela zbiór tkanka, pozostałe czynniki (linia, czas) wyjaśniają już tylko 3% zmienności
# zliczenie fragmentów, poziom /media/mj/Seagate Backup Plus Drive/anal-sekw-ps23/2023_04_20.map/q10-srt-idx
~/bin/subread-2.0.6-Linux-x86_64/bin/featureCounts -C  -B  -Q 10  -F GTF  -s 2  -T 8  --ignoreDup  -p  --countReadPairs  -t exon  -g gene_id  -a /media/mj/Seagate\ Backup\ Plus\ Drive/NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o ../../Ndwn/counts-230420n.txt  Aq10srt.bam Bq10srt.bam Cq10srt.bam Dq10srt.bam Eq10srt.bam Fq10srt.bam Gq10srt.bam Hq10srt.bam Iq10srt.bam Jq10srt.bam Kq10srt.bam Lq10srt.bam Mq10srt.bam ../../Ndwn/Ndwnq10srt.bam Oq10srt.bam Pq10srt.bam Rq10srt.bam Sq10srt.bam Tq10srt.bam Uq10srt.bam Wq10srt.bam Xq10srt.bam Yq10srt.bam Zq10srt.bam --verbose
# Statystyki zliczeń, poziom /media/mj/Seagate Backup Plus Drive/anal-sekw-ps23/Ndwn
mkdir xcounts4multiqc # muszę przenieść plik z summary, żeby tylko jego używał
cp -a counts-230420n.txt.summary xcounts4multiqc/
cd ./xcounts4multiqc/
multiqc .
