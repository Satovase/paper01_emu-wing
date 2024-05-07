# cellranger-7.1.0
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/342/905/GCF_003342905.1_droNov1
cellranger mkgtf GCF_003342905.1_droNov1_genomic.gtf droNov_refseq.filtered.gtf --attribute=gene_biotype:protein_cording
cellranger mkref --genome=droNov_refseq_genome --fasta=GCF_003342905.1_droNov1_genomic.fna --genes=droNov_refseq.filtered.gtf

# emu stage 25
cellranger count --id=run_count_EmuHH25sc_refseq \
   --fastqs=EmuHH25-ScRNA_fastqs \
   --sample=EmuHH25-ScRNA \
   --transcriptome=droNov_refseq_genome

# emu stage 20-21
cellranger count --id=run_count_EmuHH20-21sc_refseq \
   --fastqs=EmuHH20-21ScRNA_fastqs \
   --sample=EmuHH2021scRNA \
   --transcriptome=droNov_refseq_genome
