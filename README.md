# 440_Melanoma_scRNAseq
Melanoma RNA seq for 20.440 project

In order to complete all preprocessing steps to begin producing main figures run:
1. Begin with 1_Data Loading of GEO dataset
2. Proceed with 2_Metadata_Mapping_Preprocessing
3. Proceed with 3_dimensionality_reduction_harmony_QC to cluster cells **supplementary Figure 1 UMAP** and output Silouette score
4. To get the **supplementary Figure 1 heatmap**, run 4_Differential_gene_expression_sup_fig_1 immediately after 3_dimensionality_reduction_harmony_QC
6. To apply final cell typing and validation of cell typing for **supplementary figure 2**, run 6_Final_Cell_Type_Supp_Fig_2
7. Apply CellChat to uncover ligand-receptor interactions by running 7_CellChat_and_QC

To produce Main Figures:
1. Rerun 5_Final_Clustering_Fig_1 to produce all sections of **figure 1**
2. Run 8_T_cell_state_Fig_2 to produce all sections of **figure 2**
3. Run 9_Chemokines_Fig_3 to produce all sections of **figure 3**
4. Run 12_Myeloid_Interactions_Fig_4 to produce all sections of **figure 3**

To produce supplementary figures 3 and 5
1. Run 10_Supplementary_Fig_3
2. Run 11_Supplementary_Figure_5



