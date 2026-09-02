#' wompwomp alluvial plotting
#'
#' plot_alluvial() and its helpers: build an alluvial diagram from a data frame
#' with prep_for_lodes() + sort_to_uncross(), colour it, and render with ggalluvial.
#' (Merged from the former biowomp package.)
#'
#' @importFrom dplyr mutate select group_by summarise arrange desc ungroup slice n pull filter bind_rows across matches all_of add_count distinct rename
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual labs after_stat annotate theme_void theme element_text rel ggsave guides scale_color_manual scale_x_continuous element_blank coord_flip
#' @importFrom ggalluvial geom_alluvium geom_stratum stat_stratum stat_alluvium
#' @importFrom ggforce gather_set_data
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv write.csv combn
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @name plot_alluvial-internals
NULL

StatStratum <- ggalluvial::StatStratum # avoid the error Can't find stat called "stratum" - and make sure to do stat = StatStratum instead of stat = "stratum"

# `default_colors` is defined in color_alluvial.R.

load_in_df <- function(df, graphing_columns = NULL, column_weights = NULL) {
    if (is.character(df) && grepl("\\.csv$", df)) {
        df <- read.csv(df) # load in CSV as dataframe
    } else if (tibble::is_tibble(df)) {
        df <- as.data.frame(df) # convert tibble to dataframe
    } else if (!is.data.frame(df)) {
        stop("Input must be a data frame, tibble, or CSV file path.")
    }
    
    if (!is.null(column_weights)) {
        if (!(column_weights %in% colnames(df))) {
            stop(sprintf("column_weights '%s' is not a column in the dataframe.", column_weights))
        }
        df <- tidyr::uncount(df, weights = !!rlang::sym(column_weights))
    }
    
    if (!(is.null(graphing_columns))) {
        for (col in graphing_columns) {
            if (!(col %in% colnames(df))) {
                stop(sprintf("column '%s' is not a column in the dataframe.", col))
            }
            # convert to factor
            if (!is.factor(df[[col]])) {
                df[[col]] <- as.factor(df[[col]])
            }
        }
    }
    
    return(df)
}





#' Generate an Alluvial Plot with Minimal Cluster Cross-over
#'
#' Creates a two-axis alluvial plot to visualize the relationship between two categorical groupings (e.g., cluster assignments from different methods),
#' minimizing crossover and aligning clusters based on agreement.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Optional character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param column1 Optional character. Can be used along with \code{column2} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column2 Optional character. Can be used along with \code{column1} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from: 'neighbornet', 'tsp', 'greedy_wolf', 'greedy_wblf', 'random', 'none'. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'tsp' performs Traveling Salesman Problem solver from the TSP package. 'greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when \code{graphing_columns} has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_metric Character. Metric to use for determining column order. Options are "edge_crossing" (default) or "ARI". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_algorithm Character. Algorithm to use for determining column order. Options are "tsp" (default) or "neighbornet". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weighted Logical. Weighted objective
#' @param cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_wolf'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'}.
#' @param color_boxes Logical. Whether to color the strata/boxes (representing groups).
#' @param color_bands Logical. Whether to color the alluvia/edges (connecting the strata).
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_band_list Optional named list or vector of colors to override default band colors.
#' @param color_band_column Optional Character. Which column to use for coloring bands.
#' @param color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
#' @param color_band_boundary Logical. Whether or not to color boundaries between bands
#' @param coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{graphing_columns}.
#' @param coloring_algorithm_advanced_option Character. If \code{coloring_algorithm == 'advanced'}, then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
#' @param resolution Numeric If \code{coloring_algorithm == 'advanced'}, then choose resolution for the graph clustering algorithm. Affects coloring of both bands and boxes.
#' @param cutoff Numeric If \code{coloring_algorithm != 'none' and coloring_algorithm != 'advanced'}, sets the cutoff for color matching, below which a new color will be assigned.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param output_plot_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, then will not be saved.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{prep_for_lodes} function.
#' @param default_sorting Character. Default column sorting in prep_for_lodes if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param box_width Numeric between 0 and 1. Box width
#' @param text_width Numeric between 0 and 1. Text width
#' @param min_text Integer greater than 0. Min text
#' @param text_size Integer greater than 0. Text size (works whether auto_adjust_text is TRUE or FALSE).
#' @param auto_adjust_text Logical. Whether to automatically adjust text size to fit in box.
#' @param axis_text_size Integer greater than 0. Axis text size
#' @param axis_text_vjust Integer. Axis text vjust
#' @param save_height Integer greater than 0. Save height, in inches
#' @param save_width Integer greater than 0. Save width, in inches
#' @param dpi Integer greater than 0. DPI for \code{output_plot_path}, if \code{output_plot_path} is a raster image or \code{rasterise_alluvia} is TRUE
#' @param rasterise_alluvia Logical. Whether to rasterize the alluvia if \code{output_plot_path} is a PDF. Can save space if DPI low enough
#' @param keep_y_labels Keep y labels
#' @param keep_x_labels Keep x labels
#' @param box_line_width Box line width
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param add_legend Logical. If TRUE, will generate a legend of the colors of boxes and alluvial
#' @param legend_loc Character. Location of legend. Only applies if \code{add_legened == TRUE}. Choices are 'right' (default), 'left', 'bottom', 'top'
#' @param flip_xy Logical. Flip x and y (rotate plot 90 degrees).
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' p <- plot_alluvial(df,
#'     graphing_columns = c("method1", "method2"),
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp",
#'     coloring_algorithm = "right"
#' )
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' p <- plot_alluvial(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value",
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp",
#'     coloring_algorithm = "right"
#' )
#'
#' @export
plot_alluvial <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL,
                          column_weights = NULL, sorting_algorithm = "tsp",
                          optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE,
                          matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6,
                          weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6,
                          weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing",
                          column_sorting_algorithm = "tsp", weighted = TRUE, cycle_start_positions = NULL, fixed_column = NULL,
                          random_initializations = 1, color_boxes = TRUE, color_bands = FALSE,
                          color_list = NULL, color_band_list = NULL, color_band_column = NULL, color_val = NULL,
                          color_band_boundary = FALSE, coloring_algorithm = "advanced", coloring_algorithm_advanced_option = "leiden", resolution = 1, cutoff = .5,
                          alluvial_alpha = 0.5, include_labels_in_boxes = TRUE, include_axis_titles = TRUE, include_group_sizes = FALSE,
                          output_plot_path = NULL, output_df_path = NULL, preprocess_data = TRUE,
                          default_sorting = "alphabetical", box_width = 1 / 3, text_width = 1 / 4, min_text = 4, text_size = 14,
                          auto_adjust_text = TRUE, axis_text_size = 2, axis_text_vjust = 0, save_height = 6, save_width = 6, dpi = 300, rasterise_alluvia = FALSE,
                          keep_y_labels = FALSE, keep_x_labels = TRUE, 
                          box_line_width = 1, verbose = FALSE, print_params = FALSE,
                          add_legend = FALSE, legend_loc = "right", flip_xy=FALSE) {
    if (print_params) print_function_params()
    lowercase_args(c("sorting_algorithm", "column_sorting_metric", "column_sorting_algorithm", "coloring_algorithm", "coloring_algorithm_advanced_option", "default_sorting", "legend_loc"))
    
    #* Type Checking Start
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }
    
    if (is.vector(graphing_columns) && length(graphing_columns) < 2) {
        stop("graphing_columns must have at least 2 entries.")
    }
    
    if (preprocess_data) {
        if (verbose) message("Loading in data")
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
        # load_in_df() expands weighted rows, so the data is now raw per-row and
        # the old weight column is gone; prep_for_lodes() re-derives "value".
        column_weights <- NULL
    }

    if (nrow(df) == 0) {
        stop("df has no rows.")
    }
    
    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }
    
    if (ncol(df) < 2) {
        stop("Dataframe must have at least 2 columns when column_weights is NULL.")
    } else if (ncol(df) > 2) {
        if (is.null(graphing_columns) && is.null(column1) && is.null(column2)) {
            stop("graphing_columns must be specified when dataframe has more than 2 columns and column_weights is NULL.")
        }
    } else { # length 2
        if (is.null(column1) && !is.null(column2)) {
            column1 <- setdiff(colnames(df), column2)
        } else if (is.null(column2) && !is.null(column1)) {
            column2 <- setdiff(colnames(df), column1)
        } else if (is.null(column1) && is.null(column2)) {
            column1 <- colnames(df)[1]
            column2 <- colnames(df)[2]
        }
    }
    
    # if someone specifies column1/2, then use it
    if (length(graphing_columns) == 2) {
        column1 <- graphing_columns[1]
        column2 <- graphing_columns[2]
    }
    
    if (is.null(graphing_columns)) {
        graphing_columns <- c(column1, column2)
    }
    
    if (is.null(fixed_column)) {
        fixed_column <- column1
    } else if ((is.integer(fixed_column) || (is.double(fixed_column)))) {
        if (fixed_column > length(colnames(df))) {
            stop(sprintf("fixed_column index '%s' is not a column in the dataframe.", fixed_column))
        } else {
            fixed_column <- colnames(df)[fixed_column]
        }
    } else if (!(fixed_column %in% colnames(df))) {
        stop(sprintf("fixed_column '%s' is not a column in the dataframe.", fixed_column))
    }
    
    if (!is.null(color_band_column)) {
        color_bands <- TRUE
    }
    #* Type Checking End
    
    if (sorting_algorithm == "greedy_wolf") {
        default_sorting <- "fixed"
    }
    # Preprocess
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather_unsorted <- prep_for_lodes(data = df, cols = graphing_columns, wt = column_weights, color_band_column = color_band_column, default_sorting = default_sorting, do_gather_set_data = FALSE, do_add_int_columns = TRUE)
        if (is.null(column_weights)) {
            column_weights <- "value" # is set during prep_for_lodes
        }
    } else {
        clus_df_gather_unsorted <- df
    }

    if (verbose) {
        compute_alluvial_statistics(clus_df_gather = clus_df_gather_unsorted, cols = graphing_columns, wt = column_weights)
    }

    # Sort. The tuning parameters that used to be direct arguments now live in
    # sort_to_uncross_options(); optimize_column_order is expressed as
    # column_method == "none"; and sort_to_uncross() now returns the sorted data
    # frame directly rather than a list.
    if (verbose) message(sprintf("Sorting data with sorting_algorithm=%s", sorting_algorithm))
    column_method <- if (isTRUE(optimize_column_order)) column_sorting_algorithm else "none"
    df <- sort_to_uncross(
        data = clus_df_gather_unsorted,
        cols = graphing_columns,
        wt = column_weights,
        method = sorting_algorithm,
        column_method = column_method,
        weight_scalar = weight_scalar,
        fixed_column = fixed_column,
        verbose = verbose,
        options = sort_to_uncross_options(
            optimize_column_order_per_cycle = optimize_column_order_per_cycle,
            matrix_initialization_value = matrix_initialization_value,
            same_side_matrix_initialization_value = same_side_matrix_initialization_value,
            matrix_initialization_value_column_order = matrix_initialization_value_column_order,
            weight_scalar_column_order = weight_scalar_column_order,
            column_metric = column_sorting_metric,
            weighted_metric = weighted,
            cycle_start_positions = cycle_start_positions,
            random_initializations = random_initializations,
            preprocess_data = FALSE
        )
    )
    if (!is.null(output_df_path)) utils::write.csv(df, output_df_path, row.names = FALSE)
    # sort_to_uncross() returns the graphing columns first, in the sorted layer
    # order, as factors whose level order encodes the sorted stratum order, and
    # drops the col*_int helper columns. Recover the layer order and rebuild the
    # integer position columns that the plotting / colouring code below expects.
    graphing_columns <- intersect(names(df), graphing_columns)
    for (i in seq_along(graphing_columns)) {
        gcol <- df[[graphing_columns[i]]]
        if (!is.factor(gcol)) gcol <- factor(gcol)
        df[[graphing_columns[i]]] <- gcol
        df[[paste0("col", i, "_int")]] <- as.integer(gcol)
    }
    # Clean grouped copy (one row per block combination) for the colouring step,
    # before gather_set_data() reshapes df for plotting.
    df_for_colors <- df

    # Plot
    if (verbose) message("Plotting data")
    #* beginning of the old plot_alluvial_internal function
    lowercase_args(c("coloring_algorithm", "coloring_algorithm_advanced_option", "legend_loc"))
    
    if (verbose) compute_alluvial_statistics(clus_df_gather = df, cols = graphing_columns, wt = column_weights)
    
    geom_alluvium <- if (rasterise_alluvia) {
        function(...) ggrastr::rasterise(ggalluvial::geom_alluvium(...), dpi = dpi)
    } else {
        ggalluvial::geom_alluvium
    }
    
    if (column_weights != "value") {
        df <- df %>% dplyr::rename(value = !!sym(column_weights))
        column_weights <- "value"
    }
    
    df <- ggforce::gather_set_data(df, 1:2)
    if (!is.numeric(df$x)) {
        df$x <- match(as.character(df$x), graphing_columns)
    } # weird Docker issue
    df <- df[df$x == 1, ]
    
    if (!is.null(color_list)) {
        ditto_colors <- color_list
    } else {
        ditto_colors <- default_colors
    }
    
    # remove user-defined colors from available color list
    if (!(is.null(color_val))) {
        # convert named list into named vector
        if (is.list(color_val)) {
            color_val <- unlist(color_val)
        }
        ditto_colors <- ditto_colors[!(ditto_colors %in% color_val)]
    }
    
    match_colors <- (coloring_algorithm != "none") # we want to match color if coloring_algorithm is not none
    
    if (!(coloring_algorithm %in% c("none", "left", "right", "advanced", graphing_columns))) {
        stop("Invalid coloring_algorithm. Options are 'none', 'left', 'right', 'advanced', or any value in graphing_columns.")
    }
    
    # warning if coloring_algorithm is not none but color_boxes is False
    if (!color_boxes && !color_bands && match_colors) {
        if (verbose) message("Warning: color_boxes and color_bands are False but coloring_algorithm is specified. boxes will not be colored.")
    }
    
    # One color per (column, block) position. All matching methods -- "advanced"
    # (community detection on the parent-score graph), "left"/"right", and a
    # named reference column -- go through get_lode_clusters(); within-column
    # colour reuse is allowed, so a group a method over-splits shows the same
    # colour in two adjacent boxes.
    if (match_colors) {
        lode_map <- get_lode_clusters(
            df_for_colors, cols = graphing_columns, wt = column_weights,
            method = coloring_algorithm, resolution = resolution,
            options = get_lode_clusters_options(
                method_advanced_option = coloring_algorithm_advanced_option,
                cutoff = cutoff, preprocess_data = FALSE
            )
        )
        per_axis <- lode_cluster_pal(
            df_for_colors, cols = graphing_columns, mapping = lode_map,
            color_palette = ditto_colors, per_axis = TRUE
        )
        final_colors <- unlist(lapply(graphing_columns, function(cn) {
            rev(stats::setNames(per_axis[[cn]], seq_along(per_axis[[cn]])))
        }), use.names = TRUE)
    } else {
        remaining_colors <- ditto_colors
        final_colors <- c()
        for (col_group in graphing_columns) {
            num_levels <- length(levels(df[[col_group]]))
            if (length(remaining_colors) < num_levels) {
                if (verbose) message("Warning: Some colors will be recycled.")
                remaining_colors <- ditto_colors
            }
            old_colors <- remaining_colors[1:num_levels]
            final_colors <- c(final_colors, rev(old_colors))
            remaining_colors <- remaining_colors[(1 + num_levels):length(remaining_colors)]
        }
    }
    
    if (!(is.null(color_val))) {
        final_value_order <- c()
        for (col_int in seq_along(graphing_columns)) {
            int_name <- paste0("col", col_int, "_int")
            group_name <- graphing_columns[[col_int]]
            
            curr_label <- as.character(unique(df[order(df[[int_name]]), ][[group_name]]))
            
            final_value_order <- c(final_value_order, rev(curr_label))
        }
        names(final_colors) <- final_value_order
        
        for (box_val in names(color_val)) {
            if (box_val %in% names(final_colors)) {
                box_val_color <- color_val[names(color_val) == box_val][1]
                # find where value is in final colors
                val_match <- which(box_val == names(final_colors))
                
                # JR added
                to_change <- c()
                for (vm in val_match) {
                    matched_color <- final_colors[vm]
                    maybe_to_change <- which(matched_color == final_colors)
                    
                    for (matched in maybe_to_change) {
                        testing <- names(final_colors)[matched]
                        if (!(testing %in% names(color_val)) || (testing == box_val)) {
                            to_change <- c(to_change, matched)
                        }
                    }
                }
                # JR added
                
                # JR commented out
                # matched_color <- final_colors[val_match[1]]
                # maybe_to_change <- which(matched_color == final_colors)
                # to_change <- c()
                # for (matched in maybe_to_change) {
                #     testing <- names(final_colors)[matched]
                #     if (!(testing %in% names(color_val)) | (testing == box_val)) {
                #         to_change <- c(to_change, matched)
                #     }
                # }
                # JR commented out
                
                final_colors[as.integer(to_change)] <- box_val_color
            }
        }
    }
    
    remaining_colors <- ditto_colors[!(ditto_colors %in% final_colors)]
    
    final_colors_legend <- final_colors
    # generate label names
    final_label_names <- c()
    for (col_int in seq_along(graphing_columns)) {
        int_name <- paste0("col", col_int, "_int")
        group_name <- graphing_columns[[col_int]]
        
        curr_label <- as.character(unique(df[order(df[[int_name]]), ][[group_name]]))
        
        final_label_names <- c(final_label_names, rev(curr_label))
    }
    names(final_colors_legend) <- final_label_names
    
    # if color_bands, add to named list
    
    if (!is.null(color_band_list)) {
        final_colors_legend <- c(color_band_list, final_colors_legend)
    } else {
        if (!is.null(color_band_column)) {
            if (!(color_band_column %in% graphing_columns)) {
                num_levels <- length(unique(df[[color_band_column]]))
                color_band_list <- remaining_colors[1:num_levels]
                names(color_band_list) <- unique(df[[color_band_column]])
                final_colors_legend <- c(color_band_list, final_colors_legend)
            }
        }
    }
    # remove duplicate names
    final_colors_legend <- final_colors_legend[!duplicated(names(final_colors_legend))]
    # remove duplicate dims
    temp_df <- df # [1:as.integer(dim(df)[1]/2),1:dim(df)[2]]
    
    # uncomment to attempt mapping
    p <- ggplot(data = temp_df, aes(y = value), )
    for (x in seq_along(graphing_columns)) {
        p$mapping[[paste0("axis", x)]] <- sym(paste0("col", x, "_int"))
    }
    
    if (color_bands) {
        if (is.null(color_band_column)) {
            color_band_column <- graphing_columns[1]
        }
        if (color_band_boundary) {
            p <- p +
                geom_alluvium(aes(fill = !!sym(color_band_column), color = !!sym(color_band_column)),
                              alpha = alluvial_alpha
                )
        } else {
            p <- p +
                geom_alluvium(aes(fill = !!sym(color_band_column)), alpha = alluvial_alpha)
        }
        p <- p + scale_fill_manual(values = final_colors_legend)
    } else {
        if (color_band_boundary) {
            p <- p + geom_alluvium(color = "grey2", alpha = alluvial_alpha)
        } else {
            p <- p + geom_alluvium(alpha = alluvial_alpha)
        }
    }
    
    if (color_band_boundary) {
        if (add_legend) {
            p <- p + scale_color_manual(values = final_colors_legend, guide = "none")
        } else {
            p <- p + scale_color_manual(values = color_band_list, guide = "none")
        }
    }
    
    
    if (color_boxes) {
        if (add_legend) {
            p <- p + geom_stratum(width = box_width, aes(fill = after_stat(!!sym("final_label_names"))), linewidth = box_line_width) + scale_fill_manual(values = final_colors_legend)
        } else {
            p <- p + geom_stratum(width = box_width, fill = final_colors, linewidth = box_line_width)
        }
    } else {
        p <- p + geom_stratum(width = box_width, linewidth = box_line_width)
    }
    
    
    if (!(include_labels_in_boxes == FALSE)) {
        if (auto_adjust_text) {
            if (!requireNamespace("ggfittext", quietly = TRUE)) {
                stop("auto_adjust_text = TRUE needs the 'ggfittext' package; install it or set auto_adjust_text = FALSE.")
            }
            p <- p +
                ggfittext::geom_fit_text(
                    reflow = TRUE, stat = StatStratum, width = text_width, min.size = min_text, size = text_size,
                    aes(label = after_stat(final_label_names))
                )
        } else {
            p <- p +
                geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)), size = text_size)
        }
    }
    
    top_y <- 0
    for (test_x in unique(df$x)) {
        curr_y <- df %>%
            filter(x == test_x) %>%
            group_by(y) %>%
            summarise(total = sum(value), .groups = "drop") %>%
            arrange(desc(total)) %>%
            mutate(cum_y = cumsum(total)) %>%
            pull(cum_y) %>%
            max()
        top_y <- max(curr_y, top_y)
    } # top_y1 and top_y2 are probably the same
    
    if (include_axis_titles) {
        p <- p +
            scale_x_continuous(
                breaks = 1:length(graphing_columns), labels = graphing_columns,
                position = "top"
            )
    }
    
    if (include_group_sizes) {
        offset_below <- top_y * 0.075
        x <- 1
        for (col_group in graphing_columns) {
            p <- p +
                annotate("text", x = x, y = -offset_below, label = length(levels(df[[col_group]])), hjust = 0.5, size = 5) # Adjust x, y for Scanpy
            x <- x + 1
        }
    }
    
    if (flip_xy) {
        p <- p + ggplot2::coord_flip()
    }
    
    p <- p +
        theme_void() +
        theme(
            text = element_text(family = "sans"))
    
    if (add_legend) {
        p <- p + theme(legend.position = legend_loc) + labs(fill = "")
    } else {
        p <- p + theme(legend.position = "none")
    }
    
    if (keep_x_labels) {
        p <- p + theme(axis.text.x = element_text(size = axis_text_size, vjust = axis_text_vjust)) # vjust adjusts the vertical position of column titles)
    } 
    if (keep_y_labels) {
        p <- p + theme(axis.text.y = element_text(size = axis_text_size, vjust = axis_text_vjust))
    }
    
    
    if (!is.null(output_plot_path)) {
        if (verbose) message(sprintf("Saving plot to: %s", output_plot_path))
        ggsave(output_plot_path,
               plot = p,
               height = save_height, width = save_width, dpi = dpi, bg = "white"
        )
    }
    
    return(p)
}
