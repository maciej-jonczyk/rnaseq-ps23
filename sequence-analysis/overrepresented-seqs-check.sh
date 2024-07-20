# Extraction of overrepresented sequences (seqs), making fasta file and comparing seqs to databases (dbs)

# From FastQC result for each fastq file extract text file with stats
# here for four dirs = from four runs
for i in {A..Z}; do \
for j in {1..4}; do \
for k in 1 2; do \
unzip -j 2024_04_01.fastqc/${i}${j}_${k}P_fastqc.zip ${i}${j}_${k}P_fastqc/fastqc_data.txt -d /media/mj/ANTIX-LIVE/qc2024_04_01
mv /media/mj/ANTIX-LIVE/qc2024_04_01/fastqc_data.txt /media/mj/ANTIX-LIVE/qc2024_04_01/${i}${j}_${k}P_fastqc \
; done ; done ; done

for i in {A..Z}; do \
for k in 1 2; do \
unzip -j 2023_04_20.fastqc/${i}_${k}P_fastqc.zip ${i}_${k}P_fastqc/fastqc_data.txt -d /media/mj/ANTIX-LIVE/qc2023_04_20
mv /media/mj/ANTIX-LIVE/qc2023_04_20/fastqc_data.txt /media/mj/ANTIX-LIVE/qc2023_04_20/${i}_${k}P_fastqc \
; done ; done

for i in {A..Z}; do \
for j in {1..4}; do \
for k in 1 2; do \
unzip -j 2023_07_10.fastqc/${i}${j}_${k}P_fastqc.zip ${i}${j}_${k}P_fastqc/fastqc_data.txt -d /media/mj/ANTIX-LIVE/qc2023_07_10
mv /media/mj/ANTIX-LIVE/qc2023_07_10/fastqc_data.txt /media/mj/ANTIX-LIVE/qc2023_07_10/${i}${j}_${k}P_fastqc \
; done ; done ; done

for i in {A..Z}; do \
for j in {1..4}; do \
for k in 1 2; do \
unzip -j 2023_08_24.fastqc/${i}${j}_${k}P_fastqc.zip ${i}${j}_${k}P_fastqc/fastqc_data.txt -d /media/mj/ANTIX-LIVE/qc2023_08_24
mv /media/mj/ANTIX-LIVE/qc2023_08_24/fastqc_data.txt /media/mj/ANTIX-LIVE/qc2023_08_24/${i}${j}_${k}P_fastqc \
; done ; done ; done

cd qc
# extract ranges with overrepresented seqs
cat qc*/* | sed -n '/>>Over/,/>>/s/.*/&/p' > allover
# remove not needed lines
grep -v '[>#]' allover | sort -k1,1 -t'       ' -u > allover-uniq

# make fasta with row number as seq name
awk -v FS='\t' -v OFS='\t' '{print ">"NR,$1}' allover-uniq | tr '\t' '\n' 
> allover-uniq.fa

# Manual deleting poly-A and poly-T seqs. If file is big one can use grep with regex

# Exonerate, ungapped alignment
# 'blast' to NAMv5 cDNA db using Exonerate (optimized for short sequences)
~/bin/exonerate-2.2.0-x86_64/bin/exonerate --query allover-uniq.fa --target /home/mj/NAMv5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cdna.fa --showvulgar FALSE --showsugar FALSE --fsmmemory 20000 --ryo ">%S %V %em" --querytype dna --targettype dna > allover2cdna

# only hit info
grep '^>' allover2cdna > allover2cdna-hits

# only 'good' hits
awk '$11>=40 && $12>=40 && $13<=5' allover2cdna-hits > allover2cdna-hits-ok

# only gene names
cut -f5 -d' ' allover2cdna-hits | cut -f1 -d'_' | sort -u > allover-gene-ids

# retrieving gene function
grep -wf allover-gene-ids anno_fun_v45/anno-all-wide > allover-gene-fun

# Blast
# Installation according to
https://www.ncbi.nlm.nih.gov/books/NBK52640/

# Blast to NAMv5 cDNA db
# make db
~/bin/ncbi-blast-2.15.0+/bin/makeblastdb -in Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cdna.fa -dbtype nucl
# blast itself
~/bin/ncbi-blast-2.15.0+/bin/blastn -query allover-uniq.fa -db ~/NAMv5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cdna.fa -out=allover-blast-res -evalue=1e-10 -num_threads=4 -max_hsps 1 -outfmt="10 qseqid sseqid sgi evalue length pident gaps mismatch qstart sstart qend send"

# Gene function retrieval
awk -v FS=',' -v OFS='\t' '$5>=40' allover-blast-res > min40
cut -f2 -d',' min40 | cut -f1 -d'_' | sort -u > allover-blast-best40-id
grep -wf allover-blast-best40-id anno_fun_v45/anno-all-wide > allover-blast-gene-fun

# Blast to NCBI dbs
# list available dbs
~/bin/ncbi-blast-2.15.0+/bin/update_blastdb.pl --showall
# downloading needed dbs, if there is not enought space the program will stop
# if selected dbs exist they will be updated
~/bin/ncbi-blast-2.15.0+/bin/update_blastdb.pl --decompress refseq_rna 16S_ribosomal_RNA LSU_eukaryote_rRNA SSU_eukaryote_rRNA

# Downloading taxonomic db. needed for organism filtering by taxid
~/bin/ncbi-blast-2.15.0+/bin/update_blastdb.pl --decompress taxdb
# assignment of environmental variable is neccessary
# temporal
export BLASTDB='/media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA'
# permanent, NOT good idea as dbs on different disks could be used in which case variable value will be incorrect.
# adding line to ~/.bashrc file
BLASTDB='/media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA'; export BLASTDB

# Blast to entire dbs
~/bin/ncbi-blast-2.15.0+/bin/blastn -query allover-uniq.fa -db /media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA/16S_ribosomal_RNA -out=allover-blast-16S-res -evalue=1e-10 -num_threads=4 -max_hsps 10 -outfmt="10 qseqid sseqid sgi evalue length pident gaps mismatch qstart sstart qend send ssciname sallseqid sallacc stitle"

~/bin/ncbi-blast-2.15.0+/bin/blastn -query allover-uniq.fa -db /media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA/LSU_eukaryote_rRNA -out=allover-blast-LSU-res -evalue=1e-10 -num_threads=4 -max_hsps 10 -outfmt="10 qseqid sseqid sgi evalue length pident gaps mismatch qstart sstart qend send ssciname sallseqid sallacc stitle"

~/bin/ncbi-blast-2.15.0+/bin/blastn -query allover-uniq.fa -db /media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA/SSU_eukaryote_rRNA -out=allover-blast-SSU-res -evalue=1e-10 -num_threads=4 -max_hsps 10 -outfmt="10 qseqid sseqid sgi evalue length pident gaps mismatch qstart sstart qend send ssciname sallseqid sallacc stitle"

~/bin/ncbi-blast-2.15.0+/bin/blastn -query allover-uniq.fa -db /media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA/refseq_rna -out=allover-blast-refseq_rna-res -evalue=1e-10 -num_threads=4 -max_hsps 10 -outfmt="10 qseqid sseqid sgi evalue length pident gaps mismatch qstart sstart qend send ssciname sallseqid sallacc stitle"

# Blast to Zea mays seqs only
# here only ref_seq db is used as the other don't contain Zea seqs
~/bin/ncbi-blast-2.15.0+/bin/blastn -query allover-uniq.fa -taxids 4577 -db /media/mj/fe05a501-f1ce-408e-b69b-1e96e40ffbd4/blast-db-rRNA/refseq_rna -out=allover-blast-refseq_rna-res.qq -evalue=1e-10 -num_threads=4 -max_hsps 10 -outfmt="10 qseqid sseqid sgi evalue length pident gaps mismatch qstart sstart qend send ssciname sallseqid sallacc stitle"
# runs ca 50 min on machine with 4 threads and 16 Gb RAM

# Filtering and checking results - Blast to entire dbs

# concatenating results for rRNA dbs
cat allover-blast-16S-res allover-blast-LSU-res allover-blast-SSU-res > allres-rna
# extracting perfect hits (full length, 100% identity, 0 mismatches and 0 gaps) to rRNA dbs
export LC_ALL=C
awk -v FS=',' -v OFS=',' '$5==50 && $6==100 && $7==0 && $8==0' allres-rna > allres-rna-perf
# extracting query ids
cut -f1 -d',' allres-rna-perf | sort -u > ~/x
# selecting perfect hits to refseq_rna db
awk -v FS=',' -v OFS=',' '$5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res > ~/x2
# extracting only ribosomal hits
grep -i 'ribosomal' ~/x2 > ~/x3
# extracting ids for perfect hits and adding them to ids from blast to rRNA dbs
cut -f1 -d',' ~/x3 | sort -u >> ~/x
# retain only unique values
sort -u ~/x > ~/x4
# all ids of overrepresented seqs
grep '>' allover-uniq.fa | tr -d '>' > ~/x2
# ids of overrepresented seqs with non-ribosomal hits
grep -wvf ~/x4 ~/x2 > ids-no-ribo
# examining non-ribosomal hits - only 7 overrepresented seqs don't have perfect hit to rRNA or ribosomal protein
awk -v FS=',' -v OFS=',' '$1==23 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==25 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==23 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==46 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==55 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==74 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==97 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less
awk -v FS=',' -v OFS=',' '$1==114 && $5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res | less

# Filtering and checking results - Blast to refseq_rna db for maize

# extracting perfect hits
awk -v FS=',' -v OFS=',' '$5==50 && $6==100 && $7==0 && $8==0' allover-blast-refseq_rna-res.qq > perf.qq
# extracting only ribosomal hits
grep -i 'ribosomal' perf.qq > ~/x3
# extracting ids for perfect hits
cut -f1 -d',' ~/x3 | sort -u > ~/x4
# all ids of overrepresented seqs
grep '>' allover-uniq.fa | tr -d '>' > ~/x2
# ids of overrepresented seqs with non-ribosomal hits
grep -wvf ~/x4 ~/x2 > ids-no-ribo
# get ids
cat ids-no-ribo
# examining non-ribosomal hits - only 7 overrepresented seqs don't have perfect hit to rRNA or ribosomal protein
awk -v FS=',' -v OFS=',' '$1==11' perf.qq | less
# and similar for remaining ids

# Conclusion - the majority of hits are to rRNA
# I'm not excluding them, libraries were poly-A enriched so rRNA is not very abundant.
# Moreover removal of rRNA is not advised here (especially by Devon and ATpoint)
https://www.biostars.org/p/434543/
# Also posts by h.mon here
https://www.biostars.org/p/346252/#346276
https://www.biostars.org/p/346252/#346301
# most of remaining sequences align perfectly with 'NAD(P)H-quinone oxidoreductase' which according to literature is very big family of proteins

# NOTE, this analysis is not exhaustive as it bases on FastQC result which assess overrepresented seqs on subset of data
