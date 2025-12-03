#' wompwomp: Coloring alluvia
#'
#' Coloring alluvia
#' @docType package
#' @name wompwomp
#'
#' @importFrom igraph V cluster_louvain cluster_leiden E
#' @importFrom dplyr add_count mutate select group_by
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

#' Control Options for `data_color()`
#'
#' Creates a list of control parameters that modify the behavior of
#' `data_color()`. These options allow fine-grained control over the coloring
#' algorithm without cluttering the main function interface.
#'
#' @param color_list Optional named list or vector mapping values in the
#'   graphing columns to colors. Overrides default palette.
#' @param color_val Optional named list: names correspond to data values, and
#'   each element is a color to assign directly.
#' @param coloring_algorithm_advanced_option Character. When
#'   `method == "advanced"`, selects the clustering algorithm.
#'   Options: `"leiden"` (default) or `"louvain"`.
#' @param cutoff Numeric. If using a non-advanced algorithm, similarity score
#'   threshold for deciding whether to reuse a color or assign a new one.
#' @param preprocess_data Logical. If TRUE (default), preprocess data with
#'   `data_preprocess()`.
#' @param print_params Logical. If TRUE, prints parameters during execution.
#' @param load_df Internal flag; not recommended to modify.
#'
#' @return A named list of control parameters for use in `data_color()`.
#'
#' @examples
#' opts <- data_color_options(method = "advanced",
#'                            coloring_algorithm_advanced_option = "leiden")
#' df_colored <- data_color(data, cols, options = opts)
#'
#' @export
data_color_options <- function(
        color_list = NULL,
        color_val = NULL,
        coloring_algorithm_advanced_option = c("leiden", "louvain"),
        cutoff = 0.5,
        preprocess_data = TRUE,
        print_params = FALSE,
        load_df = TRUE
) {
    
    coloring_algorithm_advanced_option <- match.arg(coloring_algorithm_advanced_option)
    
    list(
        color_list = color_list,
        color_val = color_val,
        coloring_algorithm_advanced_option = coloring_algorithm_advanced_option,
        cutoff = cutoff,
        preprocess_data = preprocess_data,
        print_params = print_params,
        load_df = load_df
    )
}

#' Color alluvia
#'
#' Colors a dataframe with the algorithm specified by \code{method}.
#'
#' @param data A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) wt == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two cols).
#' (2) wt != NULL: Each row represents a combination of groupings, each column from \code{cols} represents a grouping, and the column \code{wt} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{cols}, one \code{wt}).
#' @param cols Character vector. Vector of column names from \code{data} to be used in graphing (i.e., alluvial plotting).
#' @param wt Optional character. Column name from \code{data} that contains the weights of each combination of groupings if \code{data} is in format (2) (see above).
#' @param method Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{cols}.
#' @param resolution Numeric If \code{method == 'advanced'}, then choose resolution for the graph clustering algorithm.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param options Additional arguments. See data_color_options
#'
#' @return A data frame where each row represents a combination of groupings, each column from \code{cols} represents a grouping, and the column \code{wt} ('value' if \code{wt} == NULL) represents the number of entities in that combination of groupings. There will be an additional column 'node_color'.
#'
#' @examples
#' # Example 1: data format 1
#' data <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_preprocess(data, cols = c("method1", "method2"))
#'
#' # Example 2: data format 2
#' data <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_preprocess(
#'     clus_df_gather,
#'     cols = c("method1", "method2"),
#'     wt = "value"
#' )
#'
#' @export
data_color <- function(data, cols, wt = NULL, method = "advanced", resolution = 1, verbose = FALSE, options = data_color_options()) {
    list2env(options, envir = environment())  # adds my options to be variables here
    if (print_params) print_function_params()
    # lowercase_args(c("method", "coloring_algorithm_advanced_option"))
    
    #* Type Checking Start
    valid_algorithms <- c("advanced", "left", "right", "none", "random", cols)
    if (!(method %in% valid_algorithms)) {
        stop(sprintf(
            "Invalid sorting_algorithm: '%s'. Must be one of: %s",
            method, paste(valid_algorithms, collapse = ", ")
        ))
    }
    
    if (is.vector(cols) && length(cols) < 2) {
        stop("cols must have at least 2 entries.")
    }
    
    if (!is.null(cols) && any(!cols %in% colnames(data))) {
        stop("Some cols are not present in the dataframe.")
    }
    
    if ((is.character(data)) && (load_df)) {
        if (verbose) message("Loading in data")
        data <- load_in_df(data = data, cols = cols, wt = wt)
        load_df <- FALSE
    }
    #* Type Checking End
    
    
    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- data_preprocess(data = data, cols = cols, 
                                          wt = wt, load_df = load_df, 
                                          do_gather_set_data = FALSE, do_add_int_columns = FALSE)
        for (col in cols) {
            clus_df_gather[[col]] <- factor(clus_df_gather[[col]])
        }
        if (is.null(wt)) {
            wt <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather <- data
    }
    unused_colors <- default_colors
    first <- TRUE
    if (method == "left") {
        for (col_group in cols) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            if (first) {
                temp_df <- data.frame(name = levels(clus_df_gather[[col_group]]))
                temp_df[[paste0(col_group, "_colors")]] <- unused_colors[1:num_levels]
                names(temp_df) <- c(col_group, paste0(col_group, "_colors"))
                clus_df_gather_color <- dplyr::left_join(
                    clus_df_gather,
                    temp_df,
                    by = col_group
                )
                unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0(col_group, "_colors")]])]
                old_col_group <- col_group
                first <- FALSE
            } else {
                clus_df_gather_color <- find_group2_colors(clus_df_gather_color, unused_colors, 
                                                           group1_name = old_col_group, group2_name = col_group,
                                                           cutoff = cutoff
                )
                old_col_group <- col_group
                unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[col_group]])]
            }
        }
    } else if (method == "right") {
        for (col_group in cols) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            if (first) {
                temp_df <- data.frame(name = levels(clus_df_gather[[col_group]]))
                temp_df[[paste0(col_group, "_colors")]] <- unused_colors[1:num_levels]
                names(temp_df) <- c(col_group, paste0(col_group, "_colors"))
                clus_df_gather_color <- dplyr::left_join(
                    clus_df_gather,
                    temp_df,
                    by = col_group
                )
                unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0(col_group, "_colors")]])]
                old_col_group <- col_group
                first <- FALSE
            } else {
                clus_df_gather_color <- find_group2_colors(clus_df_gather_color, unused_colors, 
                                                           group1_name = old_col_group, group2_name = col_group,
                                                           cutoff = cutoff
                )
                old_col_group <- col_group
                unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[col_group]])]
            }
        }
    } else if (method == "advanced") {
        # check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("igraph", "leidenalg"), additional_message = "do not set method to 'advanced'")
        clus_df_gather_color <- find_colors_advanced(clus_df_gather, cols, unused_colors, 
                                                     coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, 
                                                     resolution = resolution)
        return (clus_df_gather_color)
    } else {
        ref_group <- method
        num_levels <- length(levels(factor(data[[col_group]])))
        
        temp_df <- data.frame(name = levels(clus_df_gather[[col_group]]))
        temp_df[[paste0(ref_group, "_colors")]] <- unused_colors[1:num_levels]
        names(temp_df) <- c(ref_group, paste0(ref_group, "_colors"))
        clus_df_gather_color <- dplyr::left_join(
            clus_df_gather,
            temp_df,
            by = ref_group
        )
        unused_colors <- unused_colors[!(unused_colors %in% data[[paste0(ref_group, "_colors")]])]
        
        for (col_group in cols) {
            if (!(col_group == method)) {
                clus_df_gather_color <- find_group2_colors(clus_df_gather_color, unused_colors, 
                                                           group1_name = ref_group, group2_name = col_group,
                                                           cutoff = cutoff
                )
                unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0(col_group, "_colors")]])]
            }
        }
    }
    clus_df_gather_color <- clus_df_gather_color %>% dplyr::select(-!!rlang::sym(wt)) 
    
    final_df <- data.frame(axis=c(),
                           value=c(),
                           color=c())
    for (col in cols) {
        temp_df <- clus_df_gather_color[, c(col, paste0(col, '_colors'))]
        names(temp_df) <- c('value', 'color')
        temp_df[['axis']] <- col
        final_df <- rbind(final_df, temp_df)
    }
    
    final_df <- split(unique(final_df), unique(final_df)[['axis']])
    final_list <- lapply(final_df, function(subdf) {
        setNames(as.list(subdf[['color']]), subdf[['value']])
    })
    return (final_list)
}