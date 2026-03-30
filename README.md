# wompwomp

Sort k-partite graphs with node order, layer order, and node grouping optimized with a heuristic to (nearly) minimize edge crossings. Useful for improving visualizations with alluvial plots — see [biowomp](https://github.com/pachterlab/biowomp) for end-to-end alluvial plot generation with vignettes. Thank you to Cory Brunson for utilizing this algorithm in [ggalluvial](https://github.com/corybrunson/ggalluvial).

wompwomp solves the **W**eighted (permutation) **O**ptimization of **M**ultiple **P**artitions-**W**eighted (label) **O**ptimization of **M**ultiple **P**artitions (W<sub>P</sub>OMP--W<sub>L</sub>OMP) problem.

![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_before_after.png)

## wompwomp functions/commands

![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/schematic.png)

## Installation:

### CRAN - Requires system [R](https://www.r-project.org/) to be installed
```         
install.packages("wompwomp")
```

### GitHub
```
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("pachterlab/wompwomp")
```

The first time any command is run on the command line, a prompt will appear asking to install any missing R dependencies.


## Usage

The I/O for each of wompwomp's functions is as follows:

1.  prep_for_lodes: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
2.  sort_to_uncross: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
3.  compute_crossing_objective: dataframe, csv, or tibble (grouped or ungrouped) --\> list

The input table can have one of two formats:
1. Ungrouped: columns specified by graphing_columns, where each row corresponds to a separate entity
2. Grouped: columns specified by graphing_columns and column_weights, where each row corresponds to a combination of graphing_columns, and column_weights specified the number of items in this combination

Read our preprint on arXiv [here](https://doi.org/10.48550/arXiv.2509.03761). Note that some functionalities and names have been changed since the preprint. To download the version of wompwomp reflected in the preprint, please see the "arxiv_v1" branch of this repository.
Reproduce figures from our preprint at https://github.com/pachterlab/ROP_2025.