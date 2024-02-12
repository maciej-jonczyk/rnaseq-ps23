
# sprawdzić workflow
https://huoww07.github.io/Bioinformatics-for-RNA-Seq/

# dane z https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03751-1

# For convenience in some steps loops were used.

# Preprocessing of the raw reads
# 1. Quality control, FastQC program was downloaded from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# For analysis of all files in the current directory *fastq.gz expression can be used. Resulting files have added ".html" at the end.
# Rozbicie na dwie analizy, żeby trochę szybciej było
~/bin/FastQC/fastqc *1.fastq.gz --outdir=../surowe_qc &
~/bin/FastQC/fastqc *2.fastq.gz --outdir=../surowe_qc &

# Adapter and quality trimming. Trimmomatic version 0.39
# Download
git clone https://github.com/usadellab/Trimmomatic.git
# Installation
sudo apt-get install ant
cd Trimmomatic/
# Edit a build.xml file to avoid error. Accoding to https://github.com/usadellab/Trimmomatic/issues/24
# for ubuntu 22.04
In order to fix this, edit file build.xml and replace 1.5 with 1.7 on line 34
featherpad build.xml
save and exit the file
# Compilation
ant
# The program is available in 
/home/mj/bin/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar

# Program podumowujący
pip install multiqc
# do działania muszą być pliki zip po np. fastqc, same html nie są obsługiwane!
# Aktualizacja
pip install --upgrade multiqc

# Trymowanie w pliku rnaseq-ost

# O adapterach Illiminy
https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm
https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002905

# Wg analizy do artykulu z grantu. (Dla FAIRE-seq!) Dwie ostatnie opcje to quality-trimming
# dla paired

java -jar /home/sekwoja/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 18 -phred33 160518_SND393_A_L004_HDS-11_R1.fastq.gz 160518_SND393_A_L004_HDS-11_R2.fastq.gz -baseout 160518_SND393_A_L004_HDS-11minlen.fastq.gz ILLUMINACLIP:/home/sekwoja/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:1:true MAXINFO:40:0.6 MINLEN:40

# dla single end

java -jar /home/sekwoja/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 22 160927_SND393_A_L001_HDS-22_R1.fastq.gz 160927_SND393_A_L001_HDS-22minlen_R1.fastq.gz ILLUMINACLIP:/home/sekwoja/bin/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MAXINFO:40:0.6 MINLEN:40

# Info o zbiorze - sekwencjonowanie na HiSeqTM 2500 Sequencing System (Illumina, USA)
# Więcej wskazówek: https://github.com/usadellab/Trimmomatic
# Procedury z trymowania rna-seq do art. - w ucz_seqrna na ANTIX

# *************** DANE DO CWICZEN SA JUZ WYTRYMOWANE ********************

# Mapping

# Najpewniej STAR
Pobranie binarnego pliku z:
https://github.com/alexdobin/STAR
https://www.reneshbedre.com/blog/star-aligner.html
https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf

# gtf zrobiony z gff (z MaizeGDB) skryptem z pakietu AGAT

~/bin/AGAT/bin/agat_convert_sp_gff2gtf.pl --gff Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf

# Instalacja AGAT: robiłem ręcznie w g https://github.com/NBISweden/AGAT#old-school---manually
# Instrukcja do skryptów: https://agat.readthedocs.io/en/latest/tools/agat_sp_compare_two_annotations.html#

# index na dell
STAR --runThreadN 7 --runMode genomeGenerate --genomeDir star-index --genomeFastaFiles NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa --sjdbGTFfile NAMv5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf --sjdbOverhang 149

# Okazało się, że mój gtf zrobiony z gff różni się od tego z gramene (/pub/gramene/release-66/gtf/zea_mays) - biorę więc z gramene
# Wg mejla z helpdesku gramene gff z MaizeGDB i Gramene nie różnią się co do genów kodujących białka - tym bardziej trzeba używać gtf z gramene
#*******************************************************************************************************************************
# WAŻNE - gtf z gramene ma tylko numery chromosomów a fasta również "chr". W analzie naszych wyników używać poprawionego gtf.
# To wyszło dopiero przy oglądaniu wyników w IGV
# Oznacza to ponowne wykonanie indeksu STAR

sed '/^[0-9]/s/^/chr/' Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gtf > Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf
#*******************************************************************************************************************************

# Poprawiony gtf

# uruchamianie z /home/mj/bin/STAR_2.7.10b/Linux_x86_64/STAR
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir star-index --genomeFastaFiles NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa --sjdbGTFfile NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf --sjdbOverhang 149
# ok. 24 min

# Polecenie porównujące gtf-y
~/bin/AGAT/bin/agat_sp_compare_two_annotations.pl -gff1 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gtf -gff2 Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gtf -o cmp.gtf -v


# https://www.reneshbedre.com/blog/star-aligner.html
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#Sec3
# https://www.ebi.ac.uk/training/online/courses/functional-genomics-ii-common-technologies-and-data-analysis-methods/rna-sequencing/performing-a-rna-seq-experiment/data-analysis/read-mapping-or-alignment/
# https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/
# https://www.biostars.org/p/153684/
# https://mperalc.gitlab.io/bulk_RNA-seq_workshop_2021/index.html
# GATK, Picard, bamCoverage z deeptools - track do IGV
# do niektórych poleceń Picard potrzebne programy z kentUtils, pojedyncze skrypty: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

# dane z https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/ (oprócz gtf, który jest z gramene)

# Przed mapowaniem trzeba dać poniższe polecenie (pozwolenie na otwarcie większej niż omyślna (1024) liczby plików z https://github.com/alexdobin/STAR/issues/528). Inaczej daje błąd
ulimit -n 10000

# Mapowanie, od razu wszystko. Polecenie w jednym wierszu, z ukośnikami nie  działa
### Zrobić próbę bez sortowania - i tak jest filrowanie i sort później. Chyba, że któryś z programów do analizy jakościowej surowych potrzebuje sort
for i in {7454620..7454643}; do STAR --runThreadN 24 --readFilesIn ../surowe/ERR${i}_1.fastq.gz ../surowe/ERR${i}_2.fastq.gz --genomeDir ../../star-index --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ERR${i} --readFilesCommand zcat; done

# Jeśli w pętli to w prefixie trzeba dać ${i}, inaczej nadpisze wynik kolejnym (bo byłaby stała nazwa)
# ok. 8 minut na plik
# 24 pary odczytów przerobiło w niecałe 3 h 40 min


# Quality control of bam files. bamqc program was used, avaliable from https://github.com/s-andrews/BamQC
# Trzea zmodyfikować build.xml jak dla Trimmomatic
# Edit a build.xml file to avoid error. Accoding to https://github.com/usadellab/Trimmomatic/issues/24
# for ubuntu 22.04
In order to fix this, edit file build.xml and replace 1.5 with 1.7 on line 9 and 10
featherpad build.xml
save and exit the file
# Program available in ~/bin/BamQC-master/bin/bamqc
# Resulting files have added ".html" at the end.
~/bin/BamQC/bin/bamqc rnaseq_cw/map-star/*.bam --threads 24 -f NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o rnaseq_cw/qc-bam-raw/
# jest multithreading więc po jednym na raz robię
# ok. 8 min zajęło dla wszystkich 24

# Filtering, only reads with MAPQ score at least 10 were retained
# Sorting and indexing
# poziom katalogu "rnaseq_cw"
for i in {7454620..7454643}; do samtools view -@24 -bq10 map-star/ERR${i}Aligned.sortedByCoord.out.bam -o bam10mapq/ERR${i}q10.bam && samtools sort -@24 -o bam10mapq-srt/ERR${i}q10srt.bam bam10mapq/ERR${i}q10.bam && cd bam10mapq-srt && samtools index -@24 ERR${i}q10srt.bam && echo zrobione ${i} && cd .. ; done
# ??? Opcjonalnie sortować wg nazwy - tak zalecane do htseq-count (zliczanie odczytów dla transkryptów)

# po tym mozna znowu bamqc

# Dla jednego pliku 2-3 min

# W rna-seq nie usuwa się duplikatów

#************************** Powtarzalność **************************************
# PCA
# Instalacja deeptools
# wymaga pythona3, już jest na serwerze
# Instalacja pip
sudo apt-get install python3-pip
# instalacja deeptools
python3 -m pip install deeptools
# w razie jakby marudził, że za stare dpenedencies. Tu dla "packaging"
python3 -m pip install --upgrade packaging

# Instrukcja pakietu: https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html
# macierz do analiz (poszczególne elementy kodu - nazwy - wydobywane z tabelki, której źródłem było sdrf dla projektu - szczegóły: rnaseq_pca0902223.sh)
# Konieczny indeks -> samtools
multiBamSummary bins --bamfiles ERR7454620q10srt.bam ERR7454621q10srt.bam ERR7454622q10srt.bam ERR7454623q10srt.bam ERR7454624q10srt.bam ERR7454625q10srt.bam ERR7454626q10srt.bam ERR7454627q10srt.bam ERR7454628q10srt.bam ERR7454629q10srt.bam ERR7454630q10srt.bam ERR7454631q10srt.bam ERR7454632q10srt.bam ERR7454633q10srt.bam ERR7454634q10srt.bam ERR7454635q10srt.bam ERR7454636q10srt.bam ERR7454637q10srt.bam ERR7454638q10srt.bam ERR7454639q10srt.bam ERR7454640q10srt.bam ERR7454641q10srt.bam ERR7454642q10srt.bam ERR7454643q10srt.bam --labels vp15.1 vp15.2 vp15.3 vp22.1 vp22.2 vp22.3 vp29.1 vp29.2 vp29.3 vp36.1 vp36.2 vp36.3 wt15.1 wt15.2 wt15.3 wt22.1 wt22.2 wt22.3 wt29.1 wt29.2 wt29.3 wt36.1 wt36.2 wt36.3 -p 24 -o ../pca_bam/prb-seq.npz

# trwało to ok 1 godz

# PCA
plotPCA -in prb-seq.npz -l vp15.1 vp15.2 vp15.3 vp22.1 vp22.2 vp22.3 vp29.1 vp29.2 vp29.3 vp36.1 vp36.2 vp36.3 wt15.1 wt15.2 wt15.3 wt22.1 wt22.2 wt22.3 wt29.1 wt29.2 wt29.3 wt36.1 wt36.2 wt36.3 --colors '#ff9e80' '#ff9e80' '#ff9e80' '#ff6e40' '#ff6e40' '#ff6e40' '#ff3d00' '#ff3d00' '#ff3d00' '#d50000' '#d50000' '#d50000' '#84ffff' '#84ffff' '#84ffff' '#40c4ff' '#40c4ff' '#40c4ff' '#00b0ff' '#00b0ff' '#00b0ff' '#2962ff' '#2962ff' '#2962ff'  --markers 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'v' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' --transpose -o prb-pca.svg
# kolory np z https://htmlcolors.com/

# bardzo szybko robi

# heatmap
plotCorrelation --corData prb-seq.npz -c spearman -p heatmap -o spear.svg # mało czytelne
plotCorrelation --corData prb-seq.npz -c pearson -p heatmap -o pear.svg # lepsze, to robić
# można dodać --removeOutliers ale nie pomaga to

# bardzo szybko robi
#************************** Powtarzalność **************************************

# wielkość fragmentów
# poziom katalogu z bam po filtrowaniu i sortowaniu
bamPEFragmentSize -b ERR7454620q10srt.bam ERR7454621q10srt.bam ERR7454622q10srt.bam ERR7454623q10srt.bam ERR7454624q10srt.bam ERR7454625q10srt.bam ERR7454626q10srt.bam ERR7454627q10srt.bam ERR7454628q10srt.bam ERR7454629q10srt.bam ERR7454630q10srt.bam ERR7454631q10srt.bam ERR7454632q10srt.bam ERR7454633q10srt.bam ERR7454634q10srt.bam ERR7454635q10srt.bam ERR7454636q10srt.bam ERR7454637q10srt.bam ERR7454638q10srt.bam ERR7454639q10srt.bam ERR7454640q10srt.bam ERR7454641q10srt.bam ERR7454642q10srt.bam ERR7454643q10srt.bam -o frag-size.svg -p 24 --samplesLabel vp15.1 vp15.2 vp15.3 vp22.1 vp22.2 vp22.3 vp29.1 vp29.2 vp29.3 vp36.1 vp36.2 vp36.3 wt15.1 wt15.2 wt15.3 wt22.1 wt22.2 wt22.3 wt29.1 wt29.2 wt29.3 wt36.1 wt36.2 wt36.3
# b. podobna dla prób

# pokrycie, trzeba robić po kilka bam'ów inacej nieczytelne. Chyba NIE ma sensu dla rna-seq
plotCoverage -b ERR7454620q10srt.bam ERR7454621q10srt.bam ERR7454622q10srt.bam --label vp15.1 vp15.2 vp15.3 -o cover20-22.svg -p 24 -v

# Analiza jakości w picard, http://broadinstitute.github.io/picard/

# Przygotowanie plików
java -jar ~/bin/picard.jar CreateSequenceDictionary -R Zm-B73-REFERENCE-NAM-5.0.fa -O Zm-B73-REFERENCE-NAM-5.0.dict
# ok

java -jar ~/bin/picard.jar BedToIntervalList -I IN-BED -O OUT-intervalist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa -SD PLIK-dict
# tego nie robie, bo nie mam pozycji rRNA

# Analizy jakości, poziomy katalogów wyjściowych

java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I INFILE -O OUTFILE --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -CHART PLOTFILE -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa # --RIBOSOMAL_INTERVALS ZROBIONY-Z-BED

# TO MOŻE TRZEBA POPRAWIĆ - REFLAT ZROBIONY Z PLIKU GTF BEZ chr, muszę znaleźć polecenie robiące refFlat
# nie ma multithreadingu więc w kilku poleceniach robię

for i in {7454620..7454623} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454624..7454627} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454628..7454631} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454632..7454635} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454636..7454639} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454640..7454643} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &

# i dla surowych, można włączyć w drugiej zakładce - inny katalog wyjściowy

for i in {7454620..7454623} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454624..7454627} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454628..7454631} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454632..7454635} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454636..7454639} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in {7454640..7454643} ; do java -jar ~/bin/picard.jar CollectRnaSeqMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ERR${i}qual --REF_FLAT ../../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.refFlat -STRAND NONE  -CHART ERR${i}q10plot -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &

# ok. 4 min / plik
# ok. 30 min jak puszczone razem w dwóch zakładkach

# poniższe można podzielone na mniejsze części, bo jeden robi się ok 7 min
for i in 7454620 7454621 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454622 7454623 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454624 7454625 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454626 7454627 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454628 7454629 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454630 7454631 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454632 7454633 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454634 7454635 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454636 7454637 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454638 7454639 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454640 7454641 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454642 7454643 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../bam10mapq-srt/ERR${i}q10srt.bam -O ${i}alignqual -H ${i}hist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &

# i dla surowych, jw w drugiej zakładce
for i in 7454620 7454621 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454622 7454623 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454624 7454625 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454626 7454627 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454628 7454629 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454630 7454631 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454632 7454633 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454634 7454635 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454636 7454637 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454638 7454639 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454640 7454641 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &
for i in 7454642 7454643 ; do java -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics -I ../map-star/ERR${i}Aligned.sortedByCoord.out.bam -O ${i}rawalignqual -H ${i}rawhist -R ../../NAMv5/Zm-B73-REFERENCE-NAM-5.0.fa; done &

# włączone razem w dwóch zakładkach, dla odfiltrowanych ok 40 min, dla surowych ok 50 min

# Podsumowanie wyników - tabela zbiorcza
grep 'PF_BASES' ERR7454629qual > x # naglowek, wziac z któregokolwiek pliku
grep -E '^[0-9]{4,}' *qual | tr ":" "\t" > x2
cat x x2 > all_qual

# Zliczanie fragmentów: https://htseq.readthedocs.io/en/master/htseqcount.html
# Instalacja
python3 -m pip install HTSeq

# Zostaję przy featureCounts, https://pubmed.ncbi.nlm.nih.gov/24227677/, ponoć ~ 20x szybszy
# inaczej traktuje ambiguous: https://bioinformatics.cvr.ac.uk/featurecounts-or-htseq-count/

# Instalacja
# pobranie z http://subread.sourceforge.net) i rozpakowanie dystrybucji binarnej
# Instrukcje: https://subread.sourceforge.net/SubreadUsersGuide.pdf

# Zdaje się, że lepszy > 80% zliczonych (w próbnym projekcie E-MTAB-11224 mRNA było wzbogacane metodą oligo-dT. U nas też tak będzie)
# pozwala na multithreading, co na dell przy 8 procesorach dało 30 - 40 min / plik
for i in {21..43}; do ~/bin/subread-2.0.4-Linux-x86_64/bin/featureCounts -C  -B  -Q 10  -F GTF  -s 2  -T 8  --ignoreDup  -p  --countReadPairs  -t exon  -g gene_id  -a ../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o fc-counts-str2/ERR74546${i}counts.txt  bam10mapq-srt/ERR74546${i}q10srt.bam  --verbose; done

# Lepsze, wszystkie pliki na raz - od razu count matrix
~/bin/subread-2.0.4-Linux-x86_64/bin/featureCounts -C  -B  -Q 10  -F GTF  -s 2  -T 8  --ignoreDup  -p  --countReadPairs  -t exon  -g gene_id  -a ../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf -o fc-counts-str2/counts-all.txt  bam10mapq-srt/*.bam  --verbose
# ok 12 min / 24 pliki po 20 mln read-pairs

***************************** Wyjaśnienie opcji ********************************************************
# poziom rnaseq_cw
featureCounts -C # nie licz chimeric fragments (r1 i r2 mapujących się do różnych chr)
-B # tylko jak oba odczyty zmapowane
-Q 10 # min mapping quality, ustawiłem tyle ile jest domyślnie w htseq-count # W obu przypadkach zbędne bo już odfiltrowane bam'y
-F GTF # format anotacji, DEFAULT
-s 2 # strand-specific, reverse
-T 24 # liczba procesorów
--ignoreDup # ignorowanie duplikatów
-p # paired end
--countReadPairs # musi być, jak PE
-t exon # feature type DEFAULT
-g gene_id # feature type to group in meta-features DEFAULT
-a ../NAMv5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gtf
-o ERR7454620counts.txt # output file
bam10mapq-srt/ERR7454620q10srt.bam # input file
--verbose # tryb gadatliwy
*********************************************************************************************************

# Rezygnuję z htseq-counts, większość niezłapana.
# Też idzie długo ok 2 h 40 min
# gorzej niż featureCounts
__no_feature    1445948                   # ok 72 %                                                                                                                                                   
__ambiguous     35403 # ok 1,8 %

# Analiza ekspresji
# Instalacja nowego R
# https://cran.r-project.org/bin/linux/ubuntu/
# Instalacja wymaganych pakietów Ubuntu
# sudo apt-get install libcurl4-openssl-dev
# sudo apt install libxml2-dev

# >>>>>>>> New version here
# https://cran.r-project.org/src/base/R-4/
# it is proper way of installation in Trisquel, now it is R-4.3.2
# it seems taht r-base-dev is automatically installed, it is required for package installation

# Fixing access to R libraries
# https://stackoverflow.com/a/49366252/1040763

# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# i wg artykułu https://www.huber.embl.de/pub/pdf/nprot.2013.099.pdf
# Wczytanie danych do DESeq: DESeqDataSetFromHTSeq z https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
# Konstrukcja zbioru https://lashlock.github.io/compbio/R_presentation.html
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
