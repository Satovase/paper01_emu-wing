# paper01_emu-wing

This repository contains R codes used to perform analysis in the "Immobilization secondary to cell death of muscle precursors with a dual transcriptional signature contributes to the emu wing skeletal pattern" paper. doi: <https://doi.org/10.1038/s41467-024-52203-x>.

Raw single-cell transcriptome sequencing data (FASTQ) of emu is available in the DNA Data Bank of Japan (Accession numbers: DRA014432 and DRA017391; Bioproject: PRJDB16987 and PRJDB13845).

Raw single-cell transcriptome sequencing data (FASTQ) of stage 24 (SRA accession numbers: SRR14570167 to SRR14570174) and stage 27 (SRA accession numbers: SRR14570175 to SRR14570178) were downloaded from the Sequence Read Archive (SRA). These data were previously published in Feregrino C, Tschopp P. Assessing evolutionary and developmental transcriptome dynamics in homologous cell types. Dev Dyn. 2022 Sep;251(9):1472-1489. doi: <https://doi.org/10.1002/dvdy.384>.

### Preprocessing using cellranger
Emu_make_cellranger_count

Chicken_make_cellranger_count

### Creating figures using Seurat
Emu_stage20-21.R: Figure 3a, 3b, S9, S10, Supplymentary table 1

Emu_stage25.R: Figure 3d, 3e, 3g, 5a, S11, S12, Supplymentary table 2, Supplymentary table 3

Chicken_stage24.R:  Figure S13a, S13b, S13d, Supplymentary table 4

Chicken_stage27.R:  Figure S13e, S13f, S13h, Supplymentary table 5

# License
This project is covered under the Apache 2.0 License.
