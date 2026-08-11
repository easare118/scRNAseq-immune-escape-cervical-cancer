# Single-Cell RNA-Seq Analysis of Immune Escape Mechanisms in Cervical Cancer

Analysis of a publicly available cervical tumor scRNA-seq dataset (11,800+ cells) in R (Seurat), investigating PD-1/PD-L1 checkpoint-mediated immune escape and T cell exhaustion.

---

## 1. PCA — Elbow Plot

<img src="02_ElbowPlot.png" width="600">

*PCA standard deviation across the first 20 principal components. The gradual decline (rather than a sharp elbow) reflects the greater cellular heterogeneity of solid tumor tissue compared to blood. 15 PCs were retained for downstream clustering.*

---

## 2. Unsupervised Clustering

<img src="03_umap_clusters.png" width="600">

*UMAP projection showing 14 unsupervised clusters (Louvain algorithm, resolution 0.5) prior to cell-type annotation. Clusters group into distinct spatial "islands," consistent with separate major cell lineages.*

---

## 3. Marker Gene Expression (FeaturePlot)

<img src="04_FeaturePlot_umap.png" width="700">

*Expression of canonical lineage markers (EPCAM, PTPRC, CD3D, CD68, COL1A1, PECAM1) projected onto the UMAP. Used alongside the DotPlot below to assign cluster identities.*

---

## 4. Marker Gene Expression Across Clusters (DotPlot)

<img src="05_DotPlot_Identity.png" width="600">

*Expression of the same six lineage markers summarized by cluster. Dot size = percent of cells expressing the gene; color = average expression. This quantitative view was the primary evidence used to assign cluster identities.*

---

## 5. Annotated Cell Populations

<img src="06__DimPlot_umap.png" width="600">

*Final annotated UMAP showing major cell populations identified in the tumor sample: Epithelial/Tumor, T cell, Macrophage, Fibroblast, Stromal, and Immune (other). Two clusters remained unresolved.*

---

## 6. Tumor Cell Checkpoint / Antigen Presentation Markers

<img src="09_pdl1_hla_violin.png" width="700">

*CD274 (PD-L1), HLA-A, HLA-B, and B2M expression in the Epithelial/Tumor cluster. HLA-A, HLA-B, and B2M were broadly and robustly expressed (antigen presentation intact). CD274 showed heterogeneous expression: 42.6% of cells with any detectable signal, 5.75% showing robust expression.*

---

## 7. T Cell Exhaustion Score

<img src="07_Exhaustion_Score_Summary.png" width="500">

*Distribution of a composite exhaustion score (AddModuleScore, based on PDCD1, LAG3, HAVCR2, CTLA4, TOX, ENTPD1) across all T cells. Most cells cluster near baseline; a distinct minority subpopulation shows notably elevated scores, consistent with a genuine exhausted subset rather than a uniform shift.*

---

## 8. Exhaustion Score Validation

<img src="08_Exhaustion_Tcells.png" width="600">

*Composite exhaustion score plotted against raw PDCD1 (PD-1) expression per T cell (r = 0.43, moderate positive correlation). Confirms the composite score captures real exhaustion-associated signal rather than noise. Independently, 29.2% of T cells showed detectable PDCD1 expression (28.5% at a stricter threshold).*

---

## Key Finding

A defined subset of tumor cells express PD-L1 (small strongly-expressing core within a wider low-level-expressing population), while a substantial, consistently robust fraction of tumor-infiltrating T cells express PD-1 and show elevated exhaustion scores. Antigen presentation machinery remains intact. This pattern is consistent with active PD-1/PD-L1-mediated immune checkpoint escape.

**Limitations:** Single tumor sample, no biological replicates; descriptive/exploratory rather than statistically confirmatory. See [full procedure documentation](cervical_cancer_full_project_steps.md) for complete methodology and code.
