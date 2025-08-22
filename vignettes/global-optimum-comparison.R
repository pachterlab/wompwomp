## -----------------------------------------------------------------------------
reticulate::use_condaenv("wompwomp_env", required = TRUE) # !!! erase
devtools::load_all() # !!! erase

# library(wompwomp)  #!!! uncomment
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
        1, 1, 1,
        2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3,
        4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5
    ),
    cluster = c(
        2, 2, 1,
        2, 1, 1, 1, 1, 1,
        2, 3, 3, 3, 3, 3, 3,
        4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 5
    )
)
column1 <- "tissue"
column2 <- "cluster"

head(df)

## -----------------------------------------------------------------------------
clus_df_gather_sorted <- data_sort(df, column1 = "tissue", column2 = "cluster", sorting_algorithm = "greedy_wolf", fixed_column = "tissue")
crossing_edges_output <- determine_crossing_edges(clus_df_gather_sorted, column1 = "col1_int", column2 = "col2_int")
crossing_edges <- crossing_edges_output$crossing_edges_df
wompwomp_objective_greedy_wolf <- determine_weighted_layer_free_objective(crossing_edges)
print(wompwomp_objective_greedy_wolf)

## -----------------------------------------------------------------------------
p1 <- plot_alluvial(clus_df_gather_sorted, column1 = "tissue", column2 = "cluster", column_weights = "value", sorting_algorithm = "greedy_wolf", fixed_column = "tissue", color_bands = TRUE)
p1

## -----------------------------------------------------------------------------
# Sample: assume clus_df_gather_sorted has factor 'cluster' with n levels
n <- nlevels(clus_df_gather_sorted$cluster)

# Get all n! permutations of labels 1 to n
perms <- gtools::permutations(n, n)

# Original cluster levels
orig_levels <- levels(clus_df_gather_sorted$col2_int)

# Loop over permutations
objective_minimum <- Inf
objectives <- numeric(nrow(perms))
for (i in 1:nrow(perms)) {
    new_labels <- perms[i, ] # as.character(perms[i, ])
    right_map <- setNames(new_labels, orig_levels)

    # Create a relabeled copy
    clus_permuted <- clus_df_gather_sorted %>%
        mutate(
            col2_int = right_map[as.character(col2_int)]
        )

    map_dict <- setNames(new_labels, orig_levels)
    clus_permuted_objective <- determine_crossing_edges(clus_permuted, column1 = "col1_int", column2 = "col2_int", return_weighted_layer_free_objective = TRUE)

    objectives[i] <- clus_permuted_objective

    if (clus_permuted_objective < objective_minimum) {
        clus_permuted_best <- clus_df_gather_sorted %>%
            mutate(cluster = factor(cluster, levels = orig_levels, labels = new_labels))
    }
}

## -----------------------------------------------------------------------------
clus_permuted_best

## -----------------------------------------------------------------------------
p_best_wolf <- plot_alluvial(clus_permuted_best, column1 = "col1_int", column2 = "col2_int", column_weights = "value", sorting_algorithm = "none", color_bands = TRUE)
p_best_wolf

## -----------------------------------------------------------------------------
df_plot <- data.frame(
    perm_id = seq_along(objectives),
    objective = sort(objectives)
)

ggplot(df_plot, aes(x = perm_id, y = objective)) +
    geom_point() +
    geom_hline(yintercept = wompwomp_objective_greedy_wolf, linetype = "dashed", color = "blue") +
    annotate("text",
        x = Inf, y = wompwomp_objective_greedy_wolf, label = "wompwomp objective",
        hjust = 1.1, vjust = -0.5, color = "blue", size = 3
    ) +
    labs(x = "Permutation index", y = "Objective value") +
    theme_minimal()

## -----------------------------------------------------------------------------
clus_df_gather_sorted <- data_sort(df, column1 = "tissue", column2 = "cluster", sorting_algorithm = "greedy_wblf", random_initializations = 10)
crossing_edges_output <- determine_crossing_edges(clus_df_gather_sorted, column1 = "col1_int", column2 = "col2_int")
crossing_edges <- crossing_edges_output$crossing_edges_df
wompwomp_objective_greedy_wblf <- determine_weighted_layer_free_objective(crossing_edges)
print(wompwomp_objective_greedy_wblf)

## -----------------------------------------------------------------------------
p2 <- plot_alluvial(clus_df_gather_sorted, column1 = "tissue", column2 = "cluster", column_weights = "value", sorting_algorithm = "greedy_wblf", color_bands = TRUE)
p2

