# cellranger-7.1.0
cellranger mkgtf genomic.gtf droNov_refseq.filtered.gtf --attribute=gene_biotype:protein_cording
cellranger mkref --genome=droNov_refseq_genome --fasta=GCF_003342905.1_droNov1_genomic.fna --genes=droNov_refseq.filtered.gtf
cellranger count --id=run_count_EmuHH20-21sc_refseq \
   --fastqs=EmuHH20-21ScRNA_fastqs \
   --sample=EmuHH2021scRNA \
   --transcriptome=droNov_refseq_genome >output.txt

