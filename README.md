# alluvialmatch
Make alluvial plots with order and colors optimized to minimize cluster cross-over (i.e., the product of weighted edge overlap)

## Before alluvialmatch
![alt text](https://github.com/pachterlab/varseek/blob/main/figures/ggalluvial.pdf?raw=true)

## After alluvialmatch
![alt text](https://github.com/pachterlab/varseek/blob/main/figures/alluvialmatch.pdf?raw=true)

## Installation: 
`remotes::install_github("pachterlab/alluvialmatch")`

## Usage
The alluvialmatch library has a single function, plot_alluvial, which takes in a data frame, CSV, or tibble as input; and returns a ggplot2 object as output. The input table can have one of two formats: 
1) Ungrouped: columns specified by column1 and column2, where each row corresponds to a separate entity
2) Grouped: columns specified by column1, column2, and column_weights, where each row corresponds to a combination of column1 and column2, and column_weights specified the number of items in this combination

Ungrouped input
```
library("alluvialmatch")
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
