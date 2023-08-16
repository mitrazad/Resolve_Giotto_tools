# Resolve_Giotto_tools

The ResolveSeurat_to_Giotto (RStG) function converts a Seurat object from Resolve into a Giotto object for further analysis of spatial transcriptomic data.

## 1) ResolveSeurat_to_Giotto function

The ResolveSeurat_to_Giotto (RStG) function is based on the [seuratToGiotto function](https://github.com/RubD/Giotto/blob/suite/R/interoperability.R).

The RStG function was tested against an official Resolve tutorial from the Giotto suite webpage. See [test example](https://github.com/mitrazad/Resolve_Giotto_tools/blob/main/RStG_tutorial.R).

The Resolve tutorial from the Giotto suite webpage: [Resolve Bioscience Breast Cancer Subcellular](https://giottosuite.readthedocs.io/en/latest/subsections/datasets/resolve_bc_210928.html).

The RStG function uses a Seurat object created with [ResolveSeurat](https://codebase.helmholtz.cloud/resolve_tools/resolve-analysis) and converts it into a Giotto object. In the Resolve tutorial from Giotto suite, the Giotto object is created starting from subcellular polygon (e.g. cell) and point (e.g. transcript) information.

The script currently requires the [Giotto suite](https://github.com/drieslab/Giotto_site_suite) version.

Note: 
+ The coordinate systems for cells are different between Seurat and Giotto (in Seurat X and Y are inverted).
