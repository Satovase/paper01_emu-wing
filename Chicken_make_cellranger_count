# Chicken
# download fastqs from Sequence Read Archive
## chicken stage 24
for i in {67..74}
do
wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR145701${i}/SRR145701${i}
fastq-dump --gzip --split-files SRR145701${i}　>> report.txt
vdb-dump --id_range SRR145701${i} >> report.txt
mv SRR145701${i}_1.fastq.gz SRR145701${i}_S1_R1_001.fastq.gz
mv SRR145701${i}_2.fastq.gz SRR145701${i}_S1_R2_001.fastq.gz
done

mkdir SRR
mv SRR14570167 SRR14570168 SRR14570169 SRR14570170 SRR14570171 SRR14570172 SRR14570173 SRR14570174 SRR
cat SRR14570167_S1_R1_001.fastq.gz SRR14570168_S1_R1_001.fastq.gz > allfiles-cat1_S1_R1_001.fastq.gz
cat SRR14570169_S1_R1_001.fastq.gz SRR14570170_S1_R1_001.fastq.gz > allfiles-cat2_S1_R1_001.fastq.gz
cat SRR14570171_S1_R1_001.fastq.gz SRR14570172_S1_R1_001.fastq.gz > allfiles-cat3_S1_R1_001.fastq.gz
cat SRR14570173_S1_R1_001.fastq.gz SRR14570174_S1_R1_001.fastq.gz > allfiles-cat4_S1_R1_001.fastq.gz
cat allfiles-cat1_S1_R1_001.fastq.gz allfiles-cat2_S1_R1_001.fastq.gz > allfiles-cat5_S1_R1_001.fastq.gz
cat allfiles-cat3_S1_R1_001.fastq.gz allfiles-cat4_S1_R1_001.fastq.gz > allfiles-cat6_S1_R1_001.fastq.gz
cat allfiles-cat5_S1_R1_001.fastq.gz allfiles-cat6_S1_R1_001.fastq.gz > allfiles-cat7_S1_R1_001.fastq.gz
cat SRR14570167_S1_R2_001.fastq.gz SRR14570168_S1_R2_001.fastq.gz > allfiles-cat1_S1_R2_001.fastq.gz
cat SRR14570169_S1_R2_001.fastq.gz SRR14570170_S1_R2_001.fastq.gz > allfiles-cat2_S1_R2_001.fastq.gz
cat SRR14570171_S1_R2_001.fastq.gz SRR14570172_S1_R2_001.fastq.gz > allfiles-cat3_S1_R2_001.fastq.gz
cat SRR14570173_S1_R2_001.fastq.gz SRR14570174_S1_R2_001.fastq.gz > allfiles-cat4_S1_R2_001.fastq.gz
cat allfiles-cat1_S1_R2_001.fastq.gz allfiles-cat2_S1_R2_001.fastq.gz > allfiles-cat5_S1_R2_001.fastq.gz
cat allfiles-cat3_S1_R2_001.fastq.gz allfiles-cat4_S1_R2_001.fastq.gz > allfiles-cat6_S1_R2_001.fastq.gz
cat allfiles-cat5_S1_R2_001.fastq.gz allfiles-cat6_S1_R2_001.fastq.gz > allfiles-cat7_S1_R2_001.fastq.gz

for i in {1..6}
do
rm allfiles-cat${i}_S1_R1_001.fastq.gz
rm allfiles-cat${i}_S1_R2_001.fastq.gz
done

## chicken stage 27
for i in {75..78}
do
wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR145701${i}/SRR145701${i}
fastq-dump --gzip --split-files SRR145701${i} >> report.txt
vdb-dump --id_range SRR145701${i} >> report.txt
mv SRR145701${i}_1.fastq.gz SRR145701${i}_S1_R1_001.fastq.gz
mv SRR145701${i}_2.fastq.gz SRR145701${i}_S1_R2_001.fastq.gz
done

mkdir SRR
mv SRR14570175 SRR14570176 SRR14570177 SRR14570178 SRR
cat SRR14570175_S1_R1_001.fastq.gz SRR14570176_S1_R1_001.fastq.gz > allfiles-cat1_S1_R1_001.fastq.gz
cat SRR14570177_S1_R1_001.fastq.gz SRR14570178_S1_R1_001.fastq.gz > allfiles-cat2_S1_R1_001.fastq.gz
cat allfiles-cat1_S1_R1_001.fastq.gz allfiles-cat2_S1_R1_001.fastq.gz > allfiles-cat3_S1_R1_001.fastq.gz
cat SRR14570175_S1_R2_001.fastq.gz SRR14570176_S1_R2_001.fastq.gz > allfiles-cat1_S1_R2_001.fastq.gz
cat SRR14570177_S1_R2_001.fastq.gz SRR14570178_S1_R2_001.fastq.gz > allfiles-cat2_S1_R2_001.fastq.gz
cat allfiles-cat1_S1_R2_001.fastq.gz allfiles-cat2_S1_R2_001.fastq.gz > allfiles-cat3_S1_R2_001.fastq.gz

for i in {1..2}
do
rm allfiles-cat${i}_S1_R1_001.fastq.gz
rm allfiles-cat${i}_S1_R2_001.fastq.gz
done

# cellranger-7.1.0
wget https://ftp.ensembl.org/pub/release-105/fasta/gallus_gallus/dna
wget https://ftp.ensembl.org/pub/release-105/gtf/gallus_gallus 
cellranger mkgtf Gallus_gallus.GRCg6a.105.gtf.gz Gallus_gallus.GRCg6a.105.filtered.protein.gtf --attribute=gene_biotype:protein_cording
cellranger mkref --genome=Gallus.gallus_genome.1 --fasta=Gallus_gallus.GRCg6a.dna.toplevel.fa --genes=Gallus_gallus.GRCg6a.105.filtered.protein.gtf

## chicken stage 24
cellranger count --id=run_count_GgHH24sc_GSM5319711 --fastqs=GgHH24sc_GSM5319711 --sample=allfiles-cat7 --transcriptome=Gallus.gallus_genome.1

## chicken stage 27
cellranger count --id=run_count_GgHH27sc_GSM5319712 --fastqs=GgHH27sc_GSM5319712 --sample=allfiles-cat3 --transcriptome=Gallus.gallus_genome.1
