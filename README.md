# wompwomp
Make alluvial plots with node order and colors optimized to minimize edge crossings with wompwomp! wompwomp solves the **W**eighted (permutation) **O**ptimization of **M**ultiple **P**artitions--**W**eighted (label) **O**ptimization of **M**ultiple **P**artitions (W<sub>P</sub>OMP–W<sub>L</sub>OMP) problem. 

## Before wompwomp
![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_unsorted.png)

## After wompwomp
![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/wompwomp_NN.png)

## wompwomp functions/commands
![alt text](https://github.com/pachterlab/wompwomp/blob/main/figures/schematic.png)


## Installation:
### R - Requires [R](https://www.r-project.org/) to be installed
```
remotes::install_github("pachterlab/wompwomp")
wompwomp::setup_python_env()
```

### Command line - Does not require R to be installed if using conda and environment.yml.
```
git clone https://github.com/pachterlab/wompwomp
cd wompwomp
conda env create -f environment.yml && conda activate wompwomp_env
```
As an alternative to conda: `Rscript install.R`

The first time any command is run on the command line, a prompt will appear asking to install any missing R dependencies.


## Usage
The wompwomp library has four functions: plot_alluvial, greedy_wolf, determine_crossing_edges, and determine_weighted_layer_free_objective. plot_alluvial plots data as an alluvial plot using ggalluvial as a framework, running maximum weighted matching to maximize color concordance and optionally running a greedy_wolf sorting algorithm implementation of the weighted one-layer free problem. greedy_wolf runs the sorting algorithm and returns the dataframe object without plotting. determine_crossing_edges determines the edges of a graph structure that cross. determine_weighted_layer_free_objective returns the sum of products of crossing edge weights from the output of determine_crossing_edges.

The I/O for each function is as follows:

1. plot_alluvial: dataframe, csv, or tibble (grouped or ungrouped) --> plot
1. greedy_wolf: dataframe, csv, or tibble (grouped or ungrouped) --> grouped dataframe/csv (sorted)
1. determine_crossing_edges: grouped dataframe/csv --> list of lists
1. determine_weighted_layer_free_objective: list of lists (from determine_crossing_edges) --> integer

The input table can have one of two formats: 
1) Ungrouped: columns specified by column1 and column2, where each row corresponds to a separate entity
2) Grouped: columns specified by column1, column2, and column_weights, where each row corresponds to a combination of column1 and column2, and column_weights specified the number of items in this combination

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
#> 1	1	1	13	
#> 2	1	2	15	
#> 3	1	3	12	
#> 4	2	1	12	
#> 5	2	2	17	
#> 6	2	3	10	

p <- plot_alluvial(df, column_weights = "weight")
p
```



## Examples in Command Line:
```bash
./exec/wompwomp plot_alluvial --df mydata.csv --column1 method1 --column2 method2
./exec/wompwomp greedy_wolf --df mydata.csv --column1 method1 --column2 method2
./exec/wompwomp determine_crossing_edges --df tmp_files/clus_df_gather_new2.csv --column1 tissue --column2 leiden --output_df_path tmp_files/crossing.csv
./exec/wompwomp determine_weighted_layer_free_objective --df tmp_files/crossing.csv
```


For help on any command, run 
`./exec/wompwomp plot_alluvial COMMAND --help`

## See a full tutorial in our introductory vignette [wompwomp-intro.Rmd](vignettes/wompwomp-intro.Rmd)
