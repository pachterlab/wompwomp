#' wompwomp: Cluster-matching alluvial plots
#'
#' @name wompwomp-imports
#' @rdname wompwomp
#' @importFrom rlang sym .data
#' @importFrom tidyselect eval_select
#' @importFrom vctrs vec_cast

utils::globalVariables(c(
    ".data", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "total", "cum_y", "best_cluster_agreement"
))

print_function_params <- function(display_df = FALSE) {
    f <- sys.function(sys.parent())
    call <- sys.call(sys.parent())
    defaults <- as.list(formals(f))
    call_args <- as.list(call)[-1]
    
    # Evaluate args, but wrap in list() so NULL survives
    eval_args <- lapply(call_args, function(arg) {
        if (is.symbol(arg) && as.character(arg) == "NULL") {
            list(NULL)  # wrapper preserves NULL
        } else {
            list(eval(arg, envir = parent.frame()))
        }
    })
    
    formal_names <- names(defaults)
    unnamed <- which(!nzchar(names(call_args)))
    names(call_args)[unnamed] <- formal_names[unnamed]
    
    # Flatten one level
    names(eval_args) <- names(call_args)
    eval_args <- lapply(eval_args, `[[`, 1)
    
    # Manual merge (defaults first, then override)
    all_args <- defaults
    for (nm in names(eval_args)) {
        all_args[nm] <- list(eval_args[[nm]])  # assign inside list()
    }
    
    # Print
    for (nm in names(all_args)) {
        val <- all_args[[nm]]
        if (is.null(val)) {
            message(nm, " = NULL")
        } else if (is.data.frame(val) && !display_df) {
            message(nm, " = dataframe")
        } else if (length(val) == 1) {
            message(nm, " = ", val)
        } else {
            message(nm, " =")
            for (j in seq_along(val)) {
                elt_name <- names(val)[j]
                if (!is.null(elt_name) && nzchar(elt_name)) {
                    message("  ", elt_name, " : ", val[[j]])
                } else {
                    message("  [", j, "] : ", val[[j]])
                }
            }
        }
    }
}

lowercase_args <- function(arg_names) {
    for (nm in arg_names) {
        val <- get(nm, envir = parent.frame())
        if (is.character(val)) {
            assign(nm, tolower(val), envir = parent.frame())
        }
    }
}


generalized_reorder <- function(clus_df_gather, clus_df_gather_sorted, cols) {
    # 1. Reorder factors of each column in cols according to ascending order of integers in coli_int, where i represents the index of the column in cols
    ordering_columns <- paste0("col", seq_along(cols), "_int")
    missing_cols <- setdiff(ordering_columns, colnames(clus_df_gather_sorted))
    if (length(missing_cols) > 0) {
        stop(sprintf("The following ordering_columns are missing from clus_df_gather_sorted: %s", paste(missing_cols, collapse = ", ")))
    }
    
    for (i in seq_along(cols)) {
        grp_col <- cols[i]
        ord_col <- ordering_columns[i]
        
        order_vec <- clus_df_gather_sorted |>
            arrange(.data[[ord_col]]) |>
            pull(.data[[grp_col]]) |>
            unique()
        
        clus_df_gather[[grp_col]] <- factor(clus_df_gather[[grp_col]], levels = order_vec)
    }
    
    # 2. Reorder columns so cols come first, in that order
    other_cols <- setdiff(colnames(clus_df_gather), cols)
    clus_df_gather <- clus_df_gather[, c(cols, other_cols)]
    clus_df_gather <- clus_df_gather[, !(names(clus_df_gather) %in% ordering_columns), drop = FALSE]
    
    clus_df_gather
}

generalized_make_int_columns <- function(clus_df_gather, cols) {
    # Check that all columns exist
    missing_cols <- setdiff(cols, colnames(clus_df_gather))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "The following cols do not exist in clus_df_gather: %s",
            paste(missing_cols, collapse = ", ")
        ))
    }
    
    # Check that each is a factor
    not_factors <- cols[!sapply(clus_df_gather[cols], is.factor)]
    if (length(not_factors) > 0) {
        stop(sprintf(
            "The following columns are not factors (factor ordering required): %s",
            paste(not_factors, collapse = ", ")
        ))
    }
    
    # Construct integer column names
    ordering_columns <- paste0("col", seq_along(cols), "_int")
    
    # Construct integer columns
    for (i in seq_along(cols)) {
        grp_col <- cols[i]
        ord_col <- ordering_columns[i]
        
        # Mapping: factor level order → integer (1, 2, ..., n)
        lvl <- levels(clus_df_gather[[grp_col]])
        map <- setNames(seq_along(lvl), lvl)  # named vector; e.g. "B"→1, "A"→2, etc
        
        # Assign integer values based on factor order
        clus_df_gather[[ord_col]] <- unname(map[clus_df_gather[[grp_col]]])
    }
    
    # Return with new integer columns added
    clus_df_gather
}

# check_col <- function(data, col) {
#     col_quo <- rlang::enquo(col)
#     
#     # Try resolving selection
#     col_sel <- tidyselect::eval_select(col_quo, data)
#     
#     # No column returned → fail
#     if (length(col_sel) == 0) {
#         stop("Column does not exist.")
#     }
#     
#     # Only allow one column
#     if (length(col_sel) > 1) {
#         stop("Column must refer to exactly one column.")
#     }
#     
#     # Extract name
#     col_name <- names(col_sel)
#     
#     message("Column exists: ", col_name)
# }

get_col_name <- function(data, col) {
    col_quo <- as_name_selection(rlang::enquo(col))
    
    # Try resolving selection.
    # If eval_select() fails (unknown column), return NULL.
    sel <- tryCatch(
        tidyselect::eval_select(col_quo, data),
        error = function(e) return(NULL)
    )
    
    # 0 columns selected → treat as NULL
    if (is.null(sel) || length(sel) == 0) {
        return(NULL)
    }
    
    # >1 column selected → invalid
    if (length(sel) > 1) {
        stop("`col` must refer to exactly one column.")
    }
    
    names(sel)
}

check_col <- function(data, col) {
    col_quo <- as_name_selection(rlang::enquo(col))
    
    # Try resolving tidyselect. If it fails, return FALSE.
    sel <- tryCatch(
        tidyselect::eval_select(col_quo, data),
        error = function(e) return(FALSE)
    )
    
    # 0 columns selected → FALSE
    if (length(sel) == 0) {
        return(FALSE)
    }
    
    # >1 column selected → ambiguous → FALSE
    if (length(sel) > 1) {
        return(FALSE)
    }
    
    # Valid: exactly one column selected
    TRUE
}

# `cols` and `wt` are documented as character vectors but are captured with
# `enquo()` so that bare column names also work. A quosure whose expression is
# a *symbol* that resolves to a character vector is a variable holding column
# names -- the way a programmatic caller (e.g. ggalluvial) passes them -- and
# tidyselect deprecated selecting by such an "external vector". Rewrite those
# to `all_of(<value>)`, which selects by value; leave every other expression
# (bare column names, `c(...)`, `starts_with(...)`, ...) untouched.
as_name_selection <- function(quo) {
    if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
        return(quo)
    }
    expr <- rlang::quo_get_expr(quo)
    if (!rlang::is_symbol(expr)) {
        return(quo)
    }
    # A symbol that does not resolve outside the data is a bare column name, so
    # it must be left for tidyselect to evaluate against the data. Distinguish
    # that from a variable that resolves to NULL.
    unresolved <- new.env()
    val <- tryCatch(rlang::eval_tidy(quo), error = function(e) unresolved)
    if (identical(val, unresolved)) {
        return(quo)
    }
    env <- rlang::quo_get_env(quo)
    if (is.null(val)) {
        # a variable holding no column at all
        return(rlang::new_quosure(NULL, env))
    }
    if (!is.character(val)) {
        return(quo)
    }
    rlang::new_quosure(rlang::expr(tidyselect::all_of(!!val)), env)
}

# `missing_value` is the label given to the stratum of missing (`NA`) entries:
# a single string, or `NA` (where `allow_na`) to keep them as missing values.
check_missing_value <- function(missing_value, allow_na = TRUE) {
    if (length(missing_value) != 1 || !(is.character(missing_value) || is.na(missing_value))) {
        stop("`missing_value` must be a single string", if (allow_na) " or NA" else "", ".")
    }
    if (is.na(missing_value) && !allow_na) {
        stop("`missing_value` must be a single string, not NA.")
    }
    invisible(missing_value)
}

# The string that stands in for missing entries while the data are processed.
# The sorting and coloring code identifies strata by name, so it cannot carry
# real `NA`s: with `missing_value = NA` they travel under a placeholder that no
# entry of `cols` uses, which restore_missing_as_na() turns back into `NA`.
resolve_missing_label <- function(data, cols, missing_value) {
    if (!is.na(missing_value)) {
        return(as.character(missing_value))
    }
    used <- unique(unlist(lapply(intersect(cols, colnames(data)), function(col) as.character(data[[col]]))))
    used <- used[!is.na(used)]
    utils::tail(make.unique(c(used, "NA")), 1)
}

# Replaces the missing entries of each column of `cols` with `missing_label`,
# returning the columns as character vectors. as.character() comes first so
# that an explicit `NA` factor level counts as missing, which is.na() on the
# factor itself would not report.
fill_missing_values <- function(data, cols, missing_label) {
    merged <- character(0)
    for (col in cols) {
        values <- as.character(data[[col]])
        is_missing <- is.na(values)
        if (any(is_missing) && missing_label %in% values) {
            merged <- c(merged, col)
        }
        values[is_missing] <- missing_label
        data[[col]] <- values
    }
    if (length(merged) > 0) {
        warning(sprintf("Missing (NA) entries are labeled '%s', which is already a value of column(s) %s, so the two are merged into one stratum. Set `missing_value` to a label the data do not use (or to NA) to keep them apart.", missing_label, paste0("'", merged, "'", collapse = ", ")), call. = FALSE)
    }
    data
}

# Inverse of fill_missing_values() for `missing_value = NA`: the placeholder
# becomes an explicit `NA` factor level (what forcats::fct_na_value_to_level()
# produces), holding the position in the level order that the placeholder had.
restore_missing_as_na <- function(data, cols, missing_label) {
    for (col in cols) {
        values <- data[[col]]
        lv <- if (is.factor(values)) levels(values) else unique(as.character(values))
        if (!(missing_label %in% lv)) next
        values <- as.character(values)
        values[values == missing_label] <- NA
        lv[lv == missing_label] <- NA
        if (is.factor(data[[col]])) {
            data[[col]] <- factor(values, levels = lv, exclude = NULL)
        } else {
            data[[col]] <- values
        }
    }
    data
}
