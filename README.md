# Transcriptomic Analysis of Rheumatoid Arthritis

## Contents

- `Countmatrix/` - How many sequencing reads map to each gene for each sample.
- `Data_RA_processed_BAM/` - Compressed, binary version of a SAM file that stores the alignment of sequencing reads to a reference genome.
- `Data_RA_raw/` - Raw data from heathy controls and people with RA.
- 

## 🧬 Introduction

In 2020, approximately 17.6 million people worldwide were living with rheumatoid arthritis (RA) ([Black et al., 2023](https://doi.org/10.1016/S2665-9913(23)00211-4)), with an estimate of 500,000 new cases occurring each year ([Zhang et al., 2025](https://doi.org/10.1038/S41598-025-92150-1;SUBJMETA=4023,420,692;KWRD=PATHOGENESIS,RHEUMATOLOGY)). RA is a chronic autoimmune disease characterized by inflammation, driven by abnormal post-translational modification of citrullinated peptides, leading to synovial inflammation, bone and cartilage erosion and joint destruction ([Jahid et al., 2023](https://doi.org/10.31138/MJR.20230801.OO)). Despite advances in research, there is still no cure in sight ([Schett et al., 2021](https://doi.org/10.1038/S41584-020-00543-5;SUBJMETA=1670,1750,256,4023,420,498,692,700;KWRD=INFLAMMATION,PROGNOSIS,RHEUMATOID+ARTHRITIS)), however, early diagnosis can enable disease remission ([Baker et al., 2024](https://doi.org/10.1136/ARD-2024-226772)). Studies indicate that females are at higher risk, possibly due to oestrogen metabolites ([Pradeepkiran, 2019](https://doi.org/10.1016/J.JTAUTO.2019.100012); [Alpízar-Rodríguez et al., 2017](https://doi.org/10.1093/RHEUMATOLOGY/KEW318); [Jiang et al., 2024](https://doi.org/10.1136/RMDOPEN-2023-003338)). Other risk factors include genetics, smoking, obesity, lung disease and age, which affect gene expression through epigenetic changes ([Poudel et al., 2020](https://doi.org/10.1007/S11926-020-00933-4); [Ishikawa & Terao, 2020](https://doi.org/10.3390/CELLS9020475); [Yunt & Solomon, 2015](https://doi.org/10.1016/J.RDC.2014.12.004)). Compared to healthy individuals, RA patients show distinct alterations in molecular pathways due to the up- and downregulation of genes ([Hao et al., 2017](https://doi.org/10.7717/PEERJ.3078/SUPP-2)). These differences can be detected through RNA sequencing (RNA-seq), transcriptomic and Gene Ontology (GO) enrichment analysis, providing valuable insight for diagnosis and treatment. By researching RA at molecular level, this study contributes to earlier and more accurate diagnosis and to generate insights that may inform future diagnostic and therapeutic approaches. 

---

## 🔬🧪 Materials and Methods
#### Study Population

This study analysed transcriptomic data from synovial biopsy samples from 8 female individuals, including 4 RA patients (>12 months diagnosis) and 4 healthy controls, as described by [Platzer et al. (2019)](https://doi.org/10.1371/JOURNAL.PONE.0219698). The patients with RA tested positive for anti-cyclic citrullinated peptide antibodies (ACPA), while controls tested negative. Sample metadata is summarised in Table 1.

**Table 1.** Overview of sample metadata used in this study.
| Sample ID     | Age | Sex    | Condition                          |
|---------------|-----|--------|------------------------------------|
| SRR4785819    | 31  | Female | Normal                             |
| SRR4785820    | 15  | Female | Normal                             |
| SRR4785828    | 31  | Female | Normal                             |
| SRR4785831    | 42  | Female | Normal                             |
| **Average Age** | **29.8** | – | **Normal**                       |
| SRR4785979    | 54  | Female | Rheumatoid arthritis (established) |
| SRR4785980    | 66  | Female | Rheumatoid arthritis (established) |
| SRR4785986    | 60  | Female | Rheumatoid arthritis (established) |
| SRR4785988    | 59  | Female | Rheumatoid arthritis (established) |
| **Average Age** | **59.8** | – | **Rheumatoid arthritis (established)** |



#### RNA-seq Data Processing
Raw RNA sequencing data was processed and analysed using R (V4.4.1; [R Core Team 2021](https://www.R-project.org/)). RA and control samples were mapped against the GRCh38.p14 (GCA_000001405.29) human reference genome using Rsubread (V2.20.0; [Liao et al., 2019](https://doi.org/10.1093/nar/gkz114)). 

#### Gene Annotation and Filtering
An index was created, and paired-end reads were aligned. Sorted and indexed BAM files were created with Rsamtools (V2.22.0; [Morgan et al., 2023](https://bioconductor.org/packages/Rsamtools/)). Gene data (GFF3-file) were imported via Readr (V2.1.5; [Wickham et al., 2023](https://CRAN.R-project.org/package=readr)) and filtered using dplyr (V1.1.4; [Wickham et al., 2023](https://CRAN.R-project.org/package=dplyr)).

#### Differential Gene Expression Analysis
Differential expression between RA and control group was analysed using DESeq2 (V1.46.0; [Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8)), an adjusted p-value < 0.05 and |log2 fold change| >1 was considered significant.

#### Visualisation of Results
[Volcano plots](Volcano_plot/Volcano_plot_all_genes_RA.png) were generated using EnhancedVolcano (V1.24.0; [Blighe et al., 2023](https://bioconductor.org/packages/EnhancedVolcano/)). KEGGREST (V1.46.0; [Tenenbaum, 2023](https://bioconductor.org/packages/KEGGREST/)) and pathview (V1.46.0; [Luo & Brouwer, 2013](https://doi.org/10.1093/bioinformatics/btt285)) were used to make a specific KEGG-pathway. 

#### Pathway and Enrichment Analysis
GO enrichment analysis was carried out with clusterProfiler (V4.14.6l; [Wu et al., 2021](https://doi.org/10.1016/j.xinn.2021.100141)) with the org.Hs.eg.db (V3.20.0; [Carlson, 2023](https://bioconductor.org/packages/org.Hs.eg.db/)) database. Results were visualized with barplots and dotplots using ggplot2 (V3.5.2; [Wickham, 2016](https://ggplot2.tidyverse.org)).

---

## 📊 Results
#### Differential Gene Expression

Gene expression analysis revealed molecular changes in RA compared to controls. [KEGG-pathway](KEGG_pathway/hsa05323.pathview.png) analysis showed upregulation of immune-related genes such as **MHC class II molecules** and **CTLA-4**.  
Downregulated genes included **TGFβ** and **IL17**.

#### GO Enrichment Analysis

GO enrichment analysis identified multiple pathways to be significantly altered in RA. Among these, [upregulated](GO-analysis/Barplot_GO/GO_Upregulated_results_plot.png) **leukocyte-mediated immunity** (adjusted *p*-value = 2.06 x 10⁻²³, fold enrichment = 2.09) and **lymphocyte-mediated immunity** (adjusted *p*-value = 3.53 x 10⁻²⁵, fold enrichment = 2.30) were selected for further analysis due to their known relevance to RA pathology. Conversely, **axonogenesis** was significantly [downregulated](GO-analysis/Barplot_GO/GO_Downregulated_results_plot.png) (adjusted *p*-value = 3.13 x 10⁻³, fold enrichment = 1.42) and selected for additional examination. 

#### Key Genes in RA-Associated Pathways

Across the leukocyte mediated immunity pathway, genes associated with RA included **TNF** and **CD40/CD40LG**. Upregulation of **IL6** and **HLA-DRB1** was detected in the lymphocyte mediated immunity pathway. In contrast, the axonogenesis-related genes **GDNF** and **NGFR** were downregulated and selected due to their relevance in neurogenic inflammation and their role in RA pathogenesis.

---
## 💡 Conclusion 
 
This analysis shows significant alterations in immune-related and nervous system pathways in individuals with RA compared to healthy controls. The upregulation of **MHC class II molecules, CD28, CD80,** and **CTLA-4**, combined with the downregulation of **TGFβ** and **IL17**, reflects immune dysregulation in RA.

Key genes like **TNF, CD40/CD40LG, IL6,** and **HLA-DRB1**  are central to inflammation and RA pathogenesis. Additionally, the downregulation of **GDNF** and **NGFR** highlights a potential neurogenic component of RA, suggesting involvement of the nervous system in disease progression.
 
Due to the limited sample size and demographic variation, further research with a larger and more varied study population is needed to confirm these findings and ensure their reliability.
