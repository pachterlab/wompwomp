#' wompwomp: Cluster-matching alluvial plots
#'
#' @name wompwomp-imports
#' @rdname wompwomp
#' @importFrom igraph V cluster_louvain cluster_leiden E
#' @importFrom dplyr add_count mutate select group_by
#' @importFrom utils read.csv write.csv
#' @importFrom tidyselect eval_select

utils::globalVariables(c(
    ".data", ":=", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name"
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
#' @param method_advanced_option Character. When
#'   `method == "advanced"`, selects the clustering algorithm.
#'   Options: `"leiden"` (default) or `"louvain"`.
#' @param cutoff Numeric. If using a non-advanced algorithm, similarity score
#'   threshold for deciding whether to reuse a color or assign a new one.
#' @param preprocess_data Logical. If TRUE (default), preprocess data with
#'   `data_preprocess()`.
#' @param print_params Logical. If TRUE, prints parameters during execution.
#'
#' @return A named list of control parameters for use in `data_color()`.
#'
#' @examples
#' 
#' data <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' opts <- data_color_options(cutoff = 0.3, method_advanced_option = "louvain")
#' mapping <- data_color(data = data, cols = c('method1', 'method2'), options = opts)
#'
#' @export
data_color_options <- function(
        method_advanced_option = c("leiden", "louvain"),
        cutoff = 0.5,
        preprocess_data = TRUE,
        print_params = FALSE
) {
    
    method_advanced_option <- match.arg(method_advanced_option)
    
    list(
        method_advanced_option = method_advanced_option,
        cutoff = cutoff,
        preprocess_data = preprocess_data,
        print_params = print_params
    )
}

data_color_internal <- function(data, cols, wt = NULL, method = "advanced", resolution = 1, verbose = FALSE, options = NULL) {
    default_opt <- data_color_options()
    if (!is.null(options)) {
        if (!is.list(options)) stop("`options` must be a list.")
        for (nm in names(default_opt)) {
            if (!nm %in% names(options)) {
                val <- ifelse(!is.null(default_opt[[nm]]), default_opt[[nm]], "NULLTMP")
                options[[nm]] <- val
            }
        }
    } else {
        options <- default_opt
        for (nm in names(options)) {
            if (is.null(options[[nm]])) {
                options[[nm]] <- "NULLTMP"
            }
        }
    }
    for (nm in names(options)) {
        if (is.null(options[[nm]]) || options[[nm]] == "NULLTMP") {
            val <- NULL
        } else {
            val <- options[[nm]]
        }
        assign(nm, val, envir = environment())
    }
    
    if (print_params) print_function_params()
    # lowercase_args(c("method", "method_advanced_option"))

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
    #* Type Checking End
    
    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- data_preprocess(data = data, cols = cols, wt = wt, 
                                          do_gather_set_data = FALSE, do_add_int_columns = FALSE)
        if (is.null(wt) || length(wt) == 0) {
            wt <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather <- data
    }
    
    for (col in cols) {
        if (!is.factor(clus_df_gather[[col]])) {
            clus_df_gather[[col]] <- as.factor(clus_df_gather[[col]])
        }
    }
    
    unused_colors <- default_colors
    first <- TRUE
    if (method == "left") {
        for (col_group in cols) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            if (first) {
                temp_df <- data.frame(name = levels(clus_df_gather[[col_group]]))
                temp_df[[paste0(col_group, "_colors")]] <- 1:num_levels
                names(temp_df) <- c(col_group, paste0(col_group, "_colors"))
                clus_df_gather_color <- dplyr::left_join(
                    clus_df_gather,
                    temp_df,
                    by = col_group
                )
                #unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0(col_group, "_colors")]])]
                max_level <- num_levels
                old_col_group <- col_group
                first <- FALSE
            } else {
                clus_df_gather_color <- find_group2_colors(clus_df_gather_color, max_level, 
                                                           group1_name = old_col_group, group2_name = col_group,
                                                           cutoff = cutoff
                )
                old_col_group <- col_group
                max_level <- max(max_level, max(clus_df_gather_color[[paste0(col_group, "_colors")]]))
                #unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[col_group]])]
            }
        }
    } else if (method == "right") {
        for (col_group in cols) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            if (first) {
                temp_df <- data.frame(name = levels(clus_df_gather[[col_group]]))
                temp_df[[paste0(col_group, "_colors")]] <- 1:num_levels
                names(temp_df) <- c(col_group, paste0(col_group, "_colors"))
                clus_df_gather_color <- dplyr::left_join(
                    clus_df_gather,
                    temp_df,
                    by = col_group
                )
                #unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0(col_group, "_colors")]])]
                max_level <- num_levels
                old_col_group <- col_group
                first <- FALSE
            } else {
                clus_df_gather_color <- find_group2_colors(clus_df_gather_color, max_level, 
                                                           group1_name = old_col_group, group2_name = col_group,
                                                           cutoff = cutoff
                )
                old_col_group <- col_group
                max_level <- max(max_level, max(clus_df_gather_color[[paste0(col_group, "_colors")]]))
                #unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[col_group]])]
            }
        }
    } else if (method == "advanced") {
        # check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("igraph", "leidenalg"), additional_message = "do not set method to 'advanced'")
        clus_df_gather_color <- find_colors_advanced(clus_df_gather, cols, unused_colors, 
                                                     method_advanced_option = method_advanced_option, 
                                                     resolution = resolution)
        return (clus_df_gather_color)
    } else {
        ref_group <- method
        num_levels <- length(levels(factor(data[[col_group]])))
        
        temp_df <- data.frame(name = levels(clus_df_gather[[col_group]]))
        temp_df[[paste0(ref_group, "_colors")]] <- 1:num_levels
        names(temp_df) <- c(ref_group, paste0(ref_group, "_colors"))
        clus_df_gather_color <- dplyr::left_join(
            clus_df_gather,
            temp_df,
            by = ref_group
        )
        max_level <- num_levels
        
        for (col_group in cols) {
            if (!(col_group == method)) {
                clus_df_gather_color <- find_group2_colors(clus_df_gather_color, max_level, 
                                                           group1_name = ref_group, group2_name = col_group,
                                                           cutoff = cutoff
                )
                max_level <- max(max_level, max(clus_df_gather_color[[paste0(col_group, "_colors")]]))
                #unused_colors <- unused_colors[!(unused_colors %in% clus_df_gather_color[[paste0(col_group, "_colors")]])]
            }
        }
    }
    
    clus_df_gather_color <- clus_df_gather_color |> dplyr::select(-!!rlang::sym(wt)) 
    
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

find_group2_colors <- function(clus_df_gather, max_level,
                               group1_name, group2_name,
                               cutoff = .5) {
    num_levels <- length(levels(clus_df_gather[[group2_name]]))
    
    clus_df_ungrouped <- clus_df_gather[, c(group1_name, group2_name, "value",
                                            paste0(group1_name, '_colors'))]
    clus_df_filtered <- clus_df_ungrouped |>
        dplyr::add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) |>
        dplyr::select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n, !!rlang::sym(paste0(group1_name, '_colors')))
    clus_df_filtered <- dplyr::distinct(clus_df_filtered)
    colnames(clus_df_filtered) <- c(group1_name, group2_name, "value", paste0(group1_name, '_colors'))
    
    
    clus_df_filtered <- clus_df_filtered |>
        dplyr::group_by(!!rlang::sym(group1_name)) |>
        dplyr::mutate(group1_size = sum(value))
    clus_df_filtered <- clus_df_filtered |>
        dplyr::group_by(!!rlang::sym(group2_name)) |>
        dplyr::mutate(group2_size = sum(value))
    clus_df_filtered <- clus_df_filtered |>
        dplyr::group_by(!!rlang::sym(group1_name)) |>
        dplyr::mutate(weight = value / group2_size)
    
    parent_df <- clus_df_filtered |>
        dplyr::group_by(!!rlang::sym(group2_name)) |>
        dplyr::filter(weight == max(weight)) # |> select(!!rlang::sym(group1_name))
    parent_df <- parent_df[parent_df$weight > cutoff,]
    
    parent_df <- parent_df[,c(group1_name, group2_name, paste0(group1_name, '_colors'))]
    colnames(parent_df) <- c(group1_name, group2_name, paste0(group2_name, '_colors'))
    final_df <- dplyr::left_join(
        unique(clus_df_gather[, c(group2_name)]),
        parent_df[, c(group2_name, paste0(group2_name, '_colors'))],
        by = c(group2_name)
    )
    need_colors <- is.na(final_df[[paste0(group2_name, '_colors')]])
    
    final_df[[paste0(group2_name, '_colors')]][need_colors] <- (max_level+1):(max_level+sum(need_colors))#unused_colors[1:sum(need_colors)]
    
    final_df <- dplyr::left_join(
        clus_df_gather,
        final_df,
        by = c(group2_name)
    )
    
    return(final_df)
}

find_colors_advanced <- function(clus_df_gather, graphing_columns, ditto_colors = NULL, method_advanced_option = "leiden", resolution = 1) {
    if (is.null(ditto_colors)) {
        ditto_colors <- default_colors
    }
    clus_df_ungrouped <- clus_df_gather[, c(graphing_columns, "value")]
    
    first <- TRUE
    compared <- c()
    for (group1_name in graphing_columns) {
        for (group2_name in graphing_columns) {
            if (!(group1_name == group2_name)) {
                comp1 <- paste0(group1_name, group2_name)
                comp2 <- paste0(group2_name, group1_name)
                if (!(comp1 %in% compared | comp2 %in% compared)) {
                    if (first) {
                        clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        clus_df_filtered <- clus_df_filtered |>
                            dplyr::add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) |>
                            dplyr::select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        clus_df_filtered <- dplyr::distinct(clus_df_filtered)
                        colnames(clus_df_filtered) <- c("group1", "group2", "value")
                        
                        clus_df_filtered <- clus_df_filtered |>
                            dplyr::group_by(group1) |>
                            dplyr::mutate(group1_size = sum(value))
                        clus_df_filtered <- clus_df_filtered |>
                            dplyr::group_by(group1) |>
                            dplyr::mutate(group2_size = sum(value))
                        
                        clus_df_filtered <- clus_df_filtered |>
                            dplyr::group_by(group1) |>
                            dplyr::mutate(weight = value)
                        
                        clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), clus_df_filtered[["group1"]])
                        clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), clus_df_filtered[["group2"]])
                        
                        first <- FALSE
                    } else {
                        temp_clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        temp_clus_df_filtered <- temp_clus_df_filtered |>
                            dplyr::add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) |>
                            dplyr::select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        temp_clus_df_filtered <- dplyr::distinct(temp_clus_df_filtered)
                        colnames(temp_clus_df_filtered) <- c("group1", "group2", "value")
                        
                        temp_clus_df_filtered <- temp_clus_df_filtered |>
                            dplyr::group_by(group1) |>
                            dplyr::mutate(group1_size = sum(value))
                        temp_clus_df_filtered <- temp_clus_df_filtered |>
                            dplyr::group_by(group2) |>
                            dplyr::mutate(group2_size = sum(value))
                        
                        temp_clus_df_filtered <- temp_clus_df_filtered |>
                            dplyr::group_by(group1) |>
                            dplyr::mutate(weight = value) 
                        
                        temp_clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), temp_clus_df_filtered[["group1"]])
                        temp_clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), temp_clus_df_filtered[["group2"]])
                        
                        clus_df_filtered <- rbind(clus_df_filtered, temp_clus_df_filtered)
                    }
                }
            }
        }
    }
    
    clus_df_extra_filtered <- clus_df_filtered[, c("group1", "group2", "value")]
    g <- igraph::graph_from_data_frame(d = clus_df_extra_filtered, directed = FALSE)
    if (method_advanced_option == "louvain") {
        partition <- igraph::cluster_louvain(g, weights = igraph::E(g)$value, resolution = resolution)
    } else if (method_advanced_option == "leiden") {
        partition <- igraph::cluster_leiden(g, weights = igraph::E(g)$value, resolution = resolution)
    } else {
        stop(sprintf("method_advanced_option '%s' is not recognized. Please choose from 'leiden' (default) or 'louvain'.", method_advanced_option))
    }
    
    clus_df_leiden <- data.frame(group_name = partition$names, leiden = partition$membership)
    clus_df_leiden <- clus_df_leiden |> tidyr::separate_wider_delim(group_name, names = c("axis", 'value'), delim = "_")
    
    #clus_df_leiden[["leiden"]] <- unlist(Map(function(x) ditto_colors[x], clus_df_leiden$leiden))
    
    final_df <- split(clus_df_leiden, clus_df_leiden[['axis']])
    final_list <- lapply(final_df, function(subdf) {
        setNames(as.list(subdf[['leiden']]), subdf[['value']])
    })
    
    # numeric --> integer
    final_list <- lapply(final_list, function(inner) {
        lapply(inner, as.integer)
    })
    
    return(final_list)
}


convert_mapping_to_colors <- function(mapping, default_colors) {
    idx <- as.integer(mapping)
    idx_mod <- ((idx - 1) %% length(default_colors)) + 1
    out <- default_colors[idx_mod]
    names(out) <- names(mapping)
    out
}

#' Make stratum color list
#'
#' Convert color mapping made by data_color into list for scale_fill_manual
#'
#' @param data A data frame with columns corresponding to axes, with factors indicating sorting. (eg run through data_sort)
#' @param cols Character vector. Vector of column names from \code{data} to be used in graphing (i.e., alluvial plotting).
#' @param mapping List. Output from data_color.
#' @param color_palette Optional named list or vector mapping values in the graphing columns to colors. Overrides default palette.
#'
#' @return A vector of colors.
#'
#' @examples
#' # Example 1
#' data <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' cols = c("method1", "method2")
#' clus_df_gather <- data_sort(data, cols = cols, method = "tsp")
#' color_mapping <- data_color(data = clus_df_gather, cols = cols)
#' color_list <- make_stratum_color_list(data = clus_df_gather, cols = cols, mapping = color_mapping)
#'
#' @export
make_stratum_color_list <- function(data, cols, mapping, color_palette = NULL) {
    cols_tmp <- substitute(cols)
    if (is.call(cols_tmp) && cols_tmp[[1]] == 'c') {
        items <- as.list(cols_tmp)[-1]
        cols <- as.list(sapply(items, function(x) rlang::as_string(x)))
    }

    # 1. Collect all factor levels across all columns
    vals <- unique(unlist(lapply(cols, function(col) levels(data[[col]]))))
    
    # 2. Initialize named color vector
    flat_colors <- setNames(rep(NA_character_, length(vals)), vals)
    
    # 3. Fill colors from each sub-mapping
    for (col in cols) {
        if (!col %in% names(mapping)) {
            stop(sprintf("Column '%s' not found in mapping", col))
        }
        for (k in names(mapping[[col]])) {
            if (k %in% names(flat_colors)) {
                flat_colors[k] <- mapping[[col]][[k]]
            }
        }
    }
    
    # 4. Check for missing colors
    if (any(is.na(flat_colors))) {
        missing_levels <- names(flat_colors)[is.na(flat_colors)]
        stop(
            "Missing colors for: ",
            paste(missing_levels, collapse = ", "),
            "\nAdd them to stratum_to_color_mapping."
        )
    }
    
    if (is.null(color_palette)) {
        color_palette <- default_colors
    }
    flat_colors <- convert_mapping_to_colors(flat_colors, color_palette)
    
    return(flat_colors)
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
#' set.seed(429144)
#' data <- data.frame(
#'   method1 = factor(LETTERS[sample(1:3, 100, TRUE)]),
#'   method2 = factor(LETTERS[27 - sample(1:3, 100, TRUE)])
#' )
#' lapply(data, levels)
#' cluster_mapping <- data |> 
#'   data_color(cols = c(method1, method2))
#'
#' @export
data_color <- function(data, cols, wt = NULL, method = "advanced", resolution = 1, verbose = FALSE, options = NULL) {
    cols_expr <- rlang::enquo(cols)
    
    # if (missing(wt)) {
    #     col_names <- names(
    #         tidyselect::eval_select(cols_expr, data)
    #     )
    #     data <- data_preprocess(data = data, cols = col_names, do_gather_set_data = FALSE, do_add_int_columns = TRUE)
    #     wt <- "value" # is set during data_preprocess
    # }
    
    wt_expr <- rlang::enquo(wt)
    cols_pos <- tidyselect::eval_select(cols_expr, data = data)
    wt_pos <- tidyselect::eval_select(wt_expr, data = data)
    res <- rlang::set_names(
        data[c(cols_pos, wt_pos)],
        c(names(cols_pos), names(wt_pos))
    )
    
    map <- data_color_internal(res, cols = names(cols_pos), wt = names(wt_pos), method = method, resolution = resolution, verbose = verbose, options = options)
    # make_stratum_color_list(data = data, cols = names(cols_pos), mapping = map)
}
