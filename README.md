# LncGx
Shiny web application for exploring the response of the human non-coding genome to glucocorticoids (GCs). The data and figures for this application have yet to be published so only the code is included in this repository.

# Data
Total RNA-seq was performed in 9 primary human cell types (B cells, CD4+ T cells, Endothelial cells, Fibroblasts, Monocytes, Myoblasts, Neutrophils, Osteoblasts, Preadipocytes) obtained from 4 unrelated healthy volunteers. Cells were treated in vitro with the GC methylprednisolone, or with vehicle, and sampled 2 and 6 hours after treatment.

# Features
1. ***Intuitive data visualization***: Easily visualize the degree of differential expression (DE) of long non-coding RNA (lncRNA) genes in nine primary human cell types on a Manhattan-style plot based on the following user-defined parameters:
- Cell type
- Timepoint (2h or 6h after GC treatment)
- Direction of change (upregulation, downregulation, not signficant)
- Padj. value
- Log2 fold change
2. ***User-friendly plot interaction***: Allows users to click on or brush across genes of interest on the plot. Upon interaction, a downloadable table including relevant info. about the clicked/brushed genes is output.
3. ***Display and query lncRNA genes of interest***: Users may upload a .txt file or manually enter genes of interest to display on plot and retrieve info. about. The genes are highlighted on the plot and a table containing relevant info. about the genes is output. 
