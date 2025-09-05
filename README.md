# wompwomp

Make alluvial plots with node order and colors optimized to minimize edge crossings with wompwomp!

wompwomp solves the **W**eighted (permutation) **O**ptimization of **M**ultiple **P**artitions--**W**eighted (label) **O**ptimization of **M**ultiple **P**artitions (W<sub>P</sub>OMP--W<sub>L</sub>OMP) problem.

<!-- ![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_before_after.png) -->
<img src="https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_before_after.png" width="600">

## wompwomp functions/commands

![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/schematic.png)

## Installation:

### R - Requires system [R](https://www.r-project.org/) to be installed
Bioconductor (not yet released on Bioconductor - please install from GitHub)
```         
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("wompwomp")
wompwomp::setup_python_env()
```

GitHub
```
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("pachterlab/wompwomp")
wompwomp::setup_python_env()
```

### Command line - Does not require system R to be installed if using conda.

```         
git clone https://github.com/pachterlab/wompwomp
cd wompwomp
conda env create -f environment.yml  # or to avoid conda: Rscript install.R
conda activate wompwomp_env  # skip if used install.R above
remotes::install_local(".")  # or use --dev flag in commands
```

The first time any command is run on the command line, a prompt will appear asking to install any missing R dependencies.

While Python is not strictly required for use of the package, it is required for some options, including default package options (i.e., NeighborNet algorithm for sorting_algorithm == "neighbornet" or column_sorting_algorithm == "neighbornet", Leiden clustering for coloring_algorithm == "advanced", fenwick tree optimization for objective calculation).

## Usage

The I/O for each of wompwomp's functions is as follows:

1.  plot_alluvial: dataframe, csv, or tibble (grouped or ungrouped) --\> plot
2.  data_preprocess: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
3.  data_sort: dataframe, csv, or tibble (grouped or ungrouped) --\> dataframe (grouped)
4.  plot_alluvial_internal: dataframe, csv, or tibble (grouped) --\> plot
5.  determine_crossing_edges: dataframe, csv, or tibble (grouped or ungrouped) --\> list
6.  determine_weighted_layer_free_objective: dataframe, csv, or tibble (grouped or ungrouped) --\> integer

The input table can have one of two formats:
1. Ungrouped: columns specified by column1 and column2, where each row corresponds to a separate entity
2. Grouped: columns specified by column1, column2, and column_weights, where each row corresponds to a combination of column1 and column2, and column_weights specified the number of items in this combination

## Examples in R

Ungrouped input

```         
library("wompwomp")
df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
head(df)
#>   method1    method2
#> 1   1   1
#> 2   1   3
#> 3   1   2
#> 4   1   1
#> 5   2   1
#> 6   2   2

p <- plot_alluvial(df)
p
```

Grouped input

```         
set.seed(42)
raw_df <- data.frame(
    method1 = sample(1:3, 100, TRUE),
    method2 = sample(1:3, 100, TRUE)
)

# Aggregate by combination
df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "weight"))
head(df)

#>   method1    method2     weight
#> 1    1   1   13  
#> 2    1   2   15  
#> 3    1   3   12  
#> 4    2   1   12  
#> 5    2   2   17  
#> 6    2   3   10  

p <- plot_alluvial(df, column_weights = "weight")
p
```

## Examples in Command Line:

``` bash
./exec/wompwomp plot_alluvial --df mydata.csv --graphing_columns column1 column2
```

For help on any command, run `./exec/wompwomp COMMAND --help`

Notes about command line usage:
- all parameter values should be space-separted
ex. ./exec/wompwomp plot_alluvial --df data.csv, NOT --df=data.csv
- all parameters that take a single argument have identical names between R and command line, with the value immediately following the argument
ex. plot_alluvial(df=data.csv), ./exec/wompwomp plot_alluvial --df data.csv
- all parameters that take a vector/list of arguments have identical names between R and command line, with the values immediately following the argument, all separated by spaced
ex. plot_alluvial(graphing_columns=c("tissue", "cluster")), ./exec/wompwomp plot_alluvial --graphing_columns tissue cluster
- all boolean parameters are passed with the flag without any following arguments; boolean parameters that default to FALSE have identical names between R and command line, while boolean parameters that default to TRUE have "disable_" prepended to the name in the command line
ex. (note that the defaults for include_group_sizes=FALSE and include_axis_titles=TRUE): plot_alluvial(include_group_sizes=TRUE, include_axis_titles=FALSE), ./exec/wompwomp plot_alluvial --include_group_sizes --disable_include_axis_titles

## See a full tutorial in our introductory vignette [wompwomp-intro.Rmd](vignettes/wompwomp-intro.Rmd)

Read our preprint on arXiv [here](https://doi.org/10.48550/arXiv.2509.03761).
