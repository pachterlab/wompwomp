# wompwomp

## Note: This is the version of wompwomp that is reflected in version 1 of our preprint on arXiv [here](https://doi.org/10.48550/arXiv.2509.03761).
Since the release of this preprint, wompwomp on the main branch has undergone some changes in function naming and functionality - namely, we have focused wompwomp on the sorting/grouping functions, and have removed the end-to-end plot_alluvial function in favor of working more seamlessly with ggalluvial and the tidyverse directly. The following functions have changed names from this branch to the main branch:

- data_sort --\> sort_to_uncross
- data_color --\> get_lode_clusters
- make_stratum_color_list --\> lode_cluster_pal
- data_preprocess --\> prep_for_lodes
- determine_crossing_edges --\> compute_crossing_objective
- plot_alluvial --\> *deleted* (see https://github.com/pachterlab/biowomp).

In order to work with the latest vesion of wompwomp, see the main branch. In order to reproduce figures, see https://github.com/pachterlab/ROP_2025. 

wompwomp solves the **W**eighted (permutation) **O**ptimization of **M**ultiple **P**artitions-**W**eighted (label) **O**ptimization of **M**ultiple **P**artitions (W<sub>P</sub>OMP--W<sub>L</sub>OMP) problem. Sort k-partite graphs with node order, layer order, and node grouping optimized with a heuristic to (nearly) minimize edge crossings. Useful for improving visualizations with alluvial plots — see [biowomp](https://github.com/pachterlab/biowomp) for end-to-end alluvial plot generation with vignettes.

![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_before_after.png)

## wompwomp functions/commands

![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/schematic.png)

## Installation:

### R - Requires system [R](https://www.r-project.org/) to be installed
CRAN (not yet released on CRAN - please install from GitHub)
```         
install.packages("wompwomp")
```

GitHub
```
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("pachterlab/wompwomp")
```

The first time any command is run on the command line, a prompt will appear asking to install any missing R dependencies.


## Usage

The I/O for each of wompwomp's functions is as follows:

1.  data_preprocess: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
2.  data_sort: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
3.  determine_crossing_edges: dataframe, csv, or tibble (grouped or ungrouped) --\> list

The input table can have one of two formats:
1. Ungrouped: columns specified by graphing_columns, where each row corresponds to a separate entity
2. Grouped: columns specified by graphing_columns and column_weights, where each row corresponds to a combination of graphing_columns, and column_weights specified the number of items in this combination

Read our preprint on arXiv [here](https://doi.org/10.48550/arXiv.2509.03761).
See examples and vignettes [here](https://github.com/pachterlab/biowomp).
