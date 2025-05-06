# alluvialmatch
Make alluvial plots with order and colors optimized to minimize cluster cross-over

## Installation: 
`remotes::install_github("pachterlab/alluvialmatch")`

## Usage
The alluvialmatch library has a single function, plot_alluvial, which takes in a dataframe, CSV, or tibble as input; and returns a ggplot2 object as output.
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
