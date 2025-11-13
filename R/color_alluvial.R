#' wompwomp: Coloring alluvia
#'
#' Coloring alluvia
#' @docType package
#' @name wompwomp
#'
#' @importFrom igraph V cluster_louvain cluster_leiden E
#' @importFrom dplyr add_count
#' @importFrom magrittr %>%
#' @importFrom utils read.csv write.csv

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name"
))
default_colors <- c(
    "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685",
    "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2",
    "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666",
    "#3D3D3D"
)

#' Color alluvia
#'
#' Colors a dataframe with the algorithm specified by \code{coloring_algorithm}.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting).
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param color_list Logical named list or vector of colors to override default group colors.
#' @param match_colors Optional match colors
#' @param color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
#' @param coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{graphing_columns}.
#' @param coloring_algorithm_advanced_option Character. If \code{coloring_algorithm == 'advanced'}, then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
#' @param resolution Numeric If \code{coloring_algorithm == 'advanced'}, then choose resolution for the graph clustering algorithm.
#' @param cutoff Numeric If \code{coloring_algorithm != 'none' and coloring_algorithm != 'advanced'}, sets the cutoff for color matching, below which a new color will be assigned.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param load_df Internal flag; not recommended to modify.
#'
#' @return A data frame where each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} ('value' if \code{column_weights} == NULL) represents the number of entities in that combination of groupings. There will be an additional column 'node_color'.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_preprocess(df, graphing_columns = c("method1", "method2"))
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_preprocess(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value"
#' )
#'
#' @export
data_color <- function(df, graphing_columns, column_weights = NULL, color_list = NULL, match_colors = TRUE,
                       color_val = NULL, coloring_algorithm = "advanced", coloring_algorithm_advanced_option = "leiden", 
                       resolution = 1, cutoff = .5, preprocess_data = TRUE, output_df_path = NULL, verbose = FALSE, print_params = FALSE, load_df = TRUE) {
    if (print_params) print_function_params()
    lowercase_args(c("coloring_algorithm", "coloring_algorithm_advanced_option"))
    
    #* Type Checking Start
    valid_algorithms <- c("advanced", "left", "right", "none", "random", graphing_columns)
    if (!(coloring_algorithm %in% valid_algorithms)) {
        stop(sprintf(
            "Invalid sorting_algorithm: '%s'. Must be one of: %s",
            coloring_algorithm, paste(valid_algorithms, collapse = ", ")
        ))
    }
    
    if (is.vector(graphing_columns) && length(graphing_columns) < 2) {
        stop("graphing_columns must have at least 2 entries.")
    }
    
    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }
    
    if ((is.character(df)) && (load_df)) {
        if (verbose) message("Loading in data")
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
        load_df <- FALSE
    }
    #* Type Checking End
    
    
    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, 
                                          column_weights = column_weights, load_df = load_df, 
                                          do_gather_set_data = FALSE, do_add_int_columns = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather <- df
    }
    clus_df_gather <- data_sort(clus_df_gather, graphing_columns = graphing_columns, 
                                sorting_algorithm="none", column_weights = column_weights,
                                optimize_column_order = FALSE)
    #!!!! do coloring here
    if (match_colors) {
        unused_colors <- default_colors
        first <- TRUE
        if (coloring_algorithm == "left") {
            n <- 1
            for (col_group in graphing_columns) {
                num_levels <- length(levels(clus_df_gather[[col_group]]))
                if (first) {
                    temp_df <- data.frame(name = factor(1:num_levels))
                    temp_df[[paste0("col", n, "_int_colors")]] <- unused_colors[1:num_levels]
                    names(temp_df) <- c(paste0("col", n, "_int"), paste0("col", n, "_int_colors"))
                    clus_df_gather_color <- dplyr::left_join(
                        clus_df_gather,
                        temp_df,
                        by = paste0("col", n, "_int")
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0("col", n, "_int_colors")]])]
                    first <- FALSE
                } else {
                    clus_df_gather_color <- find_group2_colors(clus_df_gather_color, unused_colors, 
                                                      group1_name = paste0("col", n - 1, "_int"), group2_name = paste0("col", n, "_int"),
                                                      cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0("col", n, "_int_colors")]])]
                }
                n <- n + 1
            }
        } else if (coloring_algorithm == "right") {
            n <- length(graphing_columns)
            for (col_group in rev(graphing_columns)) {
                num_levels <- length(levels(clus_df_gather[[col_group]]))
                if (first) {
                    temp_df <- data.frame(name = factor(1:num_levels))
                    temp_df[[paste0("col", n, "_int_colors")]] <- unused_colors[1:num_levels]
                    names(temp_df) <- c(paste0("col", n, "_int"), paste0("col", n, "_int_colors"))
                    
                    clus_df_gather_color <- dplyr::left_join(
                        clus_df_gather,
                        temp_df,
                        by = paste0("col", n, "_int")
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% df[[paste0("col", n, "_int_colors")]])]
                    first <- FALSE
                } else {
                    clus_df_gather_color <- find_group2_colors(clus_df_gather_color, unused_colors, 
                                             group1_name = paste0("col", n + 1, "_int"), group2_name = paste0("col", n, "_int"),
                                             cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0("col", n, "_int_colors")]])]
                }
                n <- n - 1
            }
        } else if (coloring_algorithm == "advanced") {
            # check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("igraph", "leidenalg"), additional_message = "do not set coloring_algorithm to 'advanced'")
            clus_df_gather_color <- find_colors_advanced(clus_df_gather, unused_colors, 
                                                         graphing_columns, coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, 
                                                         resolution = resolution)
        } else {
            col_group <- coloring_algorithm
            num_levels <- length(levels(factor(df[[col_group]])))
            
            ref_group_n <- which(coloring_algorithm == graphing_columns)[[1]]
            temp_df <- data.frame(name = factor(1:num_levels))
            temp_df[[paste0("col", ref_group_n, "_int_colors")]] <- unused_colors[1:num_levels]
            names(temp_df) <- c(paste0("col", ref_group_n, "_int"), paste0("col", ref_group_n, "_int_colors"))
            clus_df_gather_color <- dplyr::left_join(
                clus_df_gather,
                temp_df,
                by = paste0("col", ref_group_n, "_int")
            )
            unused_colors <- unused_colors[!(unused_colors %in% df[[paste0("col", ref_group_n, "_int_colors")]])]
            
            for (col_group in graphing_columns) {
                if (!(col_group == coloring_algorithm)) {
                    col_group_n <- which(col_group == graphing_columns)[[1]]
                    clus_df_gather_color <- find_group2_colors(clus_df_gather_color, unused_colors, 
                                                               group1_name = paste0("col", ref_group_n, "_int"), group2_name = paste0("col", col_group_n, "_int"),
                                                               cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0("col", col_group_n, "_int_colors")]])]
                }
            }
        }
    } else {
        clus_df_gather_color <- clus_df_gather
        }
    
    # Save if desired
    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving sorted dataframe to=%s", output_df_path))
        write.csv(clus_df_gather_color, file = output_df_path, row.names = FALSE)
    }
    
    return (clus_df_gather_color)
}


find_group2_colors <- function(clus_df_gather, unused_colors,
                               group1_name = "col1_int", group2_name = "col2_int",
                               cutoff = .5) {
    num_levels <- length(levels(clus_df_gather[[group2_name]]))
    
    clus_df_ungrouped <- clus_df_gather[, c(group1_name, group2_name, "value",
                                            paste0(group1_name, '_colors'))]
    clus_df_filtered <- clus_df_ungrouped %>%
        dplyr::add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
        select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n, !!rlang::sym(paste0(group1_name, '_colors')))
    clus_df_filtered <- dplyr::distinct(clus_df_filtered)
    colnames(clus_df_filtered) <- c(group1_name, group2_name, "value", paste0(group1_name, '_colors'))
    
    
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group1_name)) %>%
        mutate(group1_size = sum(value))
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group2_name)) %>%
        mutate(group2_size = sum(value))
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group1_name)) %>%
        mutate(weight = value / group2_size)
    
    parent_df <- clus_df_filtered %>%
        group_by(!!rlang::sym(group2_name)) %>%
        dplyr::filter(weight == max(weight)) # %>% select(!!rlang::sym(group1_name))
    parent_df <- parent_df[parent_df$weight > cutoff,]

    parent_df <- parent_df[,c(group1_name, group2_name, paste0(group1_name, '_colors'))]
    colnames(parent_df) <- c(group1_name, group2_name, paste0(group2_name, '_colors'))
    
    final_df <- dplyr::left_join(
        clus_df_gather,
        parent_df,
        by = c(group1_name, group2_name)
        )
    
    need_colors <- is.na(final_df[[paste0(group2_name, '_colors')]])
    final_df[[paste0(group2_name, '_colors')]][need_colors] <- unused_colors[1:sum(need_colors)]
    
    return(final_df)
}

find_colors_advanced <- function(clus_df_gather, ditto_colors, graphing_columns, coloring_algorithm_advanced_option = "leiden", resolution = 1) {
    column_int_names <- c()
    for (col_int in seq_along(graphing_columns)) {
        int_name <- paste0("col", col_int, "_int")
        column_int_names <- c(column_int_names, int_name)
    }
    clus_df_ungrouped <- clus_df_gather[, c(column_int_names, "value")]
    
    first <- TRUE
    compared <- c()
    for (group1_name in column_int_names) {
        for (group2_name in column_int_names) {
            if (!(group1_name == group2_name)) {
                comp1 <- paste0(group1_name, group2_name)
                comp2 <- paste0(group2_name, group1_name)
                if (!(comp1 %in% compared | comp2 %in% compared)) {
                    if (first) {
                        clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        clus_df_filtered <- clus_df_filtered %>%
                            dplyr::add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
                            select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        clus_df_filtered <- dplyr::distinct(clus_df_filtered)
                        colnames(clus_df_filtered) <- c("group1", "group2", "value")
                        
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group1_size = sum(value))
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group2_size = sum(value))
                        
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(weight = value) # /group2_size)
                        
                        clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), clus_df_filtered[["group1"]])
                        clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), clus_df_filtered[["group2"]])
                        
                        first <- FALSE
                    } else {
                        temp_clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            dplyr::add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
                            select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        temp_clus_df_filtered <- dplyr::distinct(temp_clus_df_filtered)
                        colnames(temp_clus_df_filtered) <- c("group1", "group2", "value")
                        
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group1_size = sum(value))
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group2) %>%
                            mutate(group2_size = sum(value))
                        
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(weight = value) # /group2_size)
                        
                        temp_clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), temp_clus_df_filtered[["group1"]])
                        temp_clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), temp_clus_df_filtered[["group2"]])
                        
                        clus_df_filtered <- rbind(clus_df_filtered, temp_clus_df_filtered)
                    }
                }
                # compared <- c(compared, comp1, comp2)
            }
        }
    }
    
    clus_df_extra_filtered <- clus_df_filtered[, c("group1", "group2", "value")]
    g <- igraph::graph_from_data_frame(d = clus_df_extra_filtered, directed = FALSE)
    if (coloring_algorithm_advanced_option == "louvain") {
        partition <- igraph::cluster_louvain(g, weights = igraph::E(g)$value, resolution = resolution)
    } else if (coloring_algorithm_advanced_option == "leiden") {
        partition <- igraph::cluster_leiden(g, weights = igraph::E(g)$value, resolution = resolution)
    } else {
        stop(sprintf("coloring_algorithm_advanced_option '%s' is not recognized. Please choose from 'leiden' (default) or 'louvain'.", coloring_algorithm_advanced_option))
    }
    
    clus_df_leiden <- data.frame(group_name = partition$names, leiden = partition$membership)
    clus_df_leiden <- clus_df_leiden %>% tidyr::separate_wider_delim(group_name, names = c("col", "trash", "group"), delim = "_")
    clus_df_leiden <- clus_df_leiden %>% tidyr::separate_wider_delim(col, names = c("trash2", "col"), delim = "col")
    clus_df_leiden <- clus_df_leiden[, c("col", "group", "leiden")]
    clus_df_leiden[["colors"]] <- unlist(Map(function(x) ditto_colors[x], clus_df_leiden$leiden))
    
    final_df <- clus_df_gather
    final_colors <- c()
    for (col_int in unique(clus_df_leiden$col)) {
        intcol <- paste0('col', col_int, '_int')
        temp_df <- clus_df_leiden[clus_df_leiden$col == col_int, ]
        colnames(temp_df) <- c('col', intcol, 'leiden', paste0(intcol, '_colors'))
        final_df <- dplyr::left_join(
            final_df,
            temp_df[, c(intcol,  paste0(intcol, '_colors'))],
            by = intcol
        )
    }
    
    return(final_df)
}
