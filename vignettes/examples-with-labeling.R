## -----------------------------------------------------------------------------
reticulate::use_condaenv("wompwomp_env", required = TRUE) # !!! erase
devtools::load_all() # !!! erase

#library(wompwomp)  #!!! uncomment
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggforce)
library(igraph)
library(tibble)
library(tidyr)
library(reticulate)

set.seed(42)

## ----setup-env, eval = FALSE--------------------------------------------------
#  wompwomp::setup_python_env()

## ----check-python, include = FALSE--------------------------------------------
if (!reticulate::py_module_available("splitspy")) {
    message("❌ 'splitspy' not found. Run wompwomp::setup_python_env() and restart R.")
    knitr::knit_exit()
}

## -----------------------------------------------------------------------------
df <- data.frame(
    tissue = c(
        "BRAIN", "BRAIN", "BRAIN",
        "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH",
        "HEART", "HEART", "HEART", "HEART", "HEART", "HEART", "HEART",
        "T CELL", "T CELL",
        "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL"
    ),
    cluster = c(
        1, 1, 2,
        1, 2, 2, 2, 2, 2,
        1, 3, 3, 3, 3, 3, 3,
        4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4
    ),
    sex = c(
        "male", "female", "male",
        "female", "male", "female", "female", "male", "female",
        "male", "female", "male", "female", "male", "female", "male",
        "female", "male",
        "female", "male", "female", "male", "female", "male", "female", "female", "male"
    ),
    special = c(
        1, 1, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        1, 0,
        0, 0, 1, 1, 0, 1, 0, 0, 1
    )
)

# write.csv(df, file = "vignette_intro_df_ungrouped.csv", row.names = FALSE, quote = FALSE)
head(df)

## -----------------------------------------------------------------------------
p1 <- plot_alluvial(df, graphing_columns = c("tissue", "cluster", "sex"), sorting_algorithm = "none",
                    coloring_algorithm = "none")
p1

## -----------------------------------------------------------------------------
p2 <- plot_alluvial(df, graphing_columns = c("tissue", "cluster", "sex"), sorting_algorithm = "none",
                    coloring_algorithm = "none",
                    add_legend = TRUE)
p2

## -----------------------------------------------------------------------------
p3 <- plot_alluvial(df, graphing_columns = c("tissue", "cluster", "sex"), sorting_algorithm = "none",
                    coloring_algorithm = "none",
                    add_legend = TRUE,
                    text_size = 6, auto_adjust_text = FALSE, )
p3

## -----------------------------------------------------------------------------
p4 <- plot_alluvial(df, graphing_columns = c("tissue", "cluster", "sex"), sorting_algorithm = "none",
                    coloring_algorithm = "none",
                    add_legend = TRUE,
                    min_text = 10, auto_adjust_text = TRUE, )
p4

## -----------------------------------------------------------------------------
p5 <- plot_alluvial(df, graphing_columns = c("tissue", "cluster", "sex"), sorting_algorithm = "none",
                    coloring_algorithm = "none",
                    add_legend = TRUE,
                    min_text = 10, auto_adjust_text = TRUE, 
                    keep_y_labels = TRUE, axis_text_size = 10, include_axis_titles = TRUE,
                    include_group_sizes = TRUE)
p5

## -----------------------------------------------------------------------------
sessioninfo::session_info()

