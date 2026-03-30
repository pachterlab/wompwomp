# wompwomp

Sort k-partite graphs with node order, layer order, and node grouping optimized with a heuristic to (nearly) minimize edge crossings. Useful for improving visualizations with alluvial plots — see [biowomp](https://github.com/pachterlab/biowomp) for end-to-end alluvial plot generation with vignettes. Thank you to Cory Brunson for utilizing this algorithm in [ggalluvial](https://github.com/corybrunson/ggalluvial).

wompwomp solves the **W**eighted (permutation) **O**ptimization of **M**ultiple **P**artitions-**W**eighted (label) **O**ptimization of **M**ultiple **P**artitions (W<sub>P</sub>OMP--W<sub>L</sub>OMP) problem.

<!-- ![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_before_after.png) -->
<img src="https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_before_after.png" width="600">

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
2.  sort_to_uncross: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
3.  determine_crossing_edges: dataframe, csv, or tibble (grouped or ungrouped) --\> list

The input table can have one of two formats:
1. Ungrouped: columns specified by graphing_columns, where each row corresponds to a separate entity
2. Grouped: columns specified by graphing_columns and column_weights, where each row corresponds to a combination of graphing_columns, and column_weights specified the number of items in this combination

Read our preprint on arXiv [here](https://doi.org/10.48550/arXiv.2509.03761).
See examples and vignettes [here](https://github.com/pachterlab/biowomp).
