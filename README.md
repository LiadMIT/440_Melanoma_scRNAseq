# 440_Melanoma_scRNAseq
Melanoma RNA seq for 20.440 project

In order to complete all preprocessing steps to begin producing main figures run steps 1-6:
1. Begin with 1_Data Loading of GEO dataset
2. Proceed with 2_Metadata_Mapping_Preprocessing.R
3. Proceed with 3_dimensionality_reduction_harmony_QC.R to cluster cells **supplementary Figure 1 UMAP** and output Silouette score

Harmony: Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0

4. To get the **supplementary Figure 1 heatmap**, run 4_Differential_gene_expression_sup_fig_1 immediately after 3_dimensionality_reduction_harmony_QC.R
5. To apply final cell typing and validation of cell typing for **supplementary figure 2**, run 6_Final_Cell_Type_Supp_Fig_2.R
6. Apply CellChat to uncover ligand-receptor interactions by running 7_CellChat_and_QC.R

CellChat: Jin, S., Plikus, M.V. & Nie, Q. CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nat Protoc 20, 180–219 (2025).       https://doi.org/10.1038/s41596-024-01045-4

To produce Main Figures:
1. Rerun 5_Final_Clustering_Fig_1.R to produce all sections of **figure 1**
2. Run 8_T_cell_state_Fig_2.R to produce all sections of **figure 2**
3. Run 9_Chemokines_Fig_3.R to produce all sections of **figure 3**
4. Run 12_Myeloid_Interactions_Fig_4.R to produce all sections of **figure 4**

To produce supplementary figures 3 and 4
1. Run 10_Supplementary_Fig_3.R to produce **supplementary figure 3**
2. Run 11_Supplementary_Figure_5.R to produce **supplementary figure 4**

To produce Supplementary tables run Tables.R

To investigate other CD45+ Cell interactions in an unsupervise CircosPlot, please run Unsupervised_CircosPlots.R



