# Missing (`NA`) entries of `cols` are kept as a stratum of their own, labeled
# by `missing_value`. The default label "Missing" can already be a value of the
# data (as in ggalluvial's `vaccinations`), in which case the two would merge.

# `a` has NAs and the value "Missing"; `b` has NAs only; `c` has "Missing" only
missing_counts <- function() {
    data.frame(
        a = c("x", "x", "y", "Missing", NA, NA, "y", "Missing"),
        b = c("p", "q", NA, "p", "q", "p", "q", NA),
        c = c("u", "v", "u", "Missing", "v", "u", "v", "u"),
        value = c(5, 3, 4, 2, 6, 1, 7, 2),
        stringsAsFactors = FALSE
    )
}
cols <- c("a", "b", "c")

# weight per stratum of `col`, with the missing stratum under the name "<NA>"
stratum_weights <- function(data, col) {
    key <- as.character(data[[col]])
    key[is.na(key)] <- "<NA>"
    weights <- tapply(data$value, key, sum)
    stats::setNames(as.numeric(weights), names(weights))
}

test_that("the default label warns when it merges with an existing value", {
    data <- missing_counts()
    expect_warning(
        sorted <- sort_to_uncross(data, cols = cols, wt = "value", method = "greedy"),
        "already a value of column\\(s\\) 'a'"
    )
    expect_false(anyNA(sorted[cols]))
    expect_equal(unname(stratum_weights(sorted, "a")["Missing"]), 2 + 6 + 1 + 2)

    # no warning without a collision: NAs but no "Missing", and vice versa
    expect_no_warning(sort_to_uncross(data, cols = c("b", "c"), wt = "value", method = "greedy"))
})

test_that("a custom missing_value keeps missing entries apart", {
    data <- missing_counts()
    for (method in c("neighbornet", "tsp", "greedy", "barycenter", "median", "none", "random")) {
        expect_no_warning(sorted <- sort_to_uncross(
            data, cols = cols, wt = "value", method = method,
            options = sort_to_uncross_options(missing_value = "(missing)")
        ))
        expect_false(anyNA(sorted[cols]))
        expect_setequal(levels(sorted$a), c("x", "y", "Missing", "(missing)"))
        expect_setequal(levels(sorted$b), c("p", "q", "(missing)"))
        expect_setequal(levels(sorted$c), c("u", "v", "Missing"))
        weights <- stratum_weights(sorted, "a")
        expect_equal(unname(weights[c("Missing", "(missing)")]), c(2 + 2, 6 + 1))
    }
})

test_that("missing_value = NA returns an explicit NA level", {
    data <- missing_counts()
    # "NA" as a value must not be confused with the missing entries
    data$b[data$b == "p"] <- "NA"
    for (method in c("neighbornet", "tsp", "greedy", "barycenter", "median", "none", "random")) {
        expect_no_warning(sorted <- sort_to_uncross(
            data, cols = cols, wt = "value", method = method,
            options = list(missing_value = NA)
        ))
        for (col in cols) {
            expect_true(is.factor(sorted[[col]]))
            # every row is assigned a level, so nothing is dropped downstream
            expect_false(anyNA(as.integer(sorted[[col]])))
            expect_equal(stratum_weights(sorted, col), stratum_weights(data, col))
        }
        expect_setequal(levels(sorted$a), c("x", "y", "Missing", NA))
        expect_setequal(levels(sorted$b), c("NA", "q", NA))
        expect_setequal(levels(sorted$c), c("u", "v", "Missing"))
    }
})

test_that("missing_value = NA matches the sort under a label", {
    data <- missing_counts()
    labeled <- sort_to_uncross(data, cols = cols, wt = "value", method = "greedy", column_method = "none",
                               options = list(missing_value = "NA"))
    as_na <- sort_to_uncross(data, cols = cols, wt = "value", method = "greedy", column_method = "none",
                             options = list(missing_value = NA))
    for (col in cols) {
        lv <- levels(labeled[[col]])
        lv[lv == "NA"] <- NA
        expect_identical(levels(as_na[[col]]), lv)
    }
})

test_that("missing_value works on uncounted data", {
    data <- missing_counts()
    data <- data[rep(seq_len(nrow(data)), data$value), cols]
    expect_no_warning(sorted <- sort_to_uncross(data, cols = cols, method = "greedy",
                                                options = list(missing_value = NA)))
    expect_equal(sum(sorted$value), nrow(data))
    expect_equal(sum(sorted$value[is.na(as.character(sorted$a))]), 6 + 1)
    expect_equal(sum(sorted$value[as.character(sorted$a) %in% "Missing"]), 2 + 2)
})

test_that("prep_for_lodes() labels missing entries", {
    data <- missing_counts()
    expect_warning(prep_for_lodes(data, cols = cols, wt = "value"), "already a value")
    prepped <- prep_for_lodes(data, cols = cols, wt = "value", missing_value = "none")
    expect_equal(sum(prepped$a == "none"), 2)
    expect_equal(sum(prepped$a == "Missing"), 2)

    # an explicit NA level (e.g. the output of missing_value = NA) counts as missing
    data$b <- addNA(factor(data$b))
    prepped <- prep_for_lodes(data, cols = cols, wt = "value", missing_value = "none")
    expect_equal(sum(prepped$b == "none"), 2)

    prepped <- prep_for_lodes(data, cols = cols, wt = "value", missing_value = NA, do_add_int_columns = TRUE)
    expect_equal(sum(is.na(as.character(prepped$a))), 2)
    expect_true(anyNA(levels(prepped$a)))
    expect_false(anyNA(prepped$col1_int))

    expect_error(prep_for_lodes(data, cols = cols, wt = "value", missing_value = c("a", "b")), "missing_value")
    expect_error(sort_to_uncross_options(missing_value = 1), "missing_value")
    expect_error(get_lode_clusters_options(missing_value = NA), "missing_value")
})

test_that("coloring names the missing stratum by missing_value", {
    data <- missing_counts()
    expect_warning(get_lode_clusters(data, cols = cols, wt = "value"), "already a value")
    expect_no_warning(mapping <- get_lode_clusters(
        data, cols = cols, wt = "value",
        options = get_lode_clusters_options(missing_value = "(missing)")
    ))
    expect_setequal(names(mapping$a), c("x", "y", "Missing", "(missing)"))
    expect_setequal(names(mapping$b), c("p", "q", "(missing)"))

    # the NA level of a sorted data frame is colored under the same name, and
    # lode_cluster_pal() leaves it to the scale's `na.value`
    sorted <- sort_to_uncross(data, cols = cols, wt = "value", method = "greedy", options = list(missing_value = NA))
    mapping_sorted <- get_lode_clusters(sorted, cols = cols, wt = "value", options = list(missing_value = "(missing)"))
    expect_setequal(names(mapping_sorted$a), names(mapping$a))
    pal <- lode_cluster_pal(sorted, cols = cols, mapping = mapping_sorted)
    expect_false(anyNA(names(pal)))
    expect_false(anyNA(pal))
    per_axis <- lode_cluster_pal(sorted, cols = cols, mapping = mapping_sorted, per_axis = TRUE)
    expect_setequal(names(per_axis$a), c("x", "y", "Missing"))

    expect_no_warning(agreement <- compute_color_agreement(data, cols = cols, mapping = mapping, wt = "value", missing_value = "(missing)"))
    from_sorted <- compute_color_agreement(sorted, cols = cols, mapping = mapping, wt = "value", missing_value = "(missing)")
    expect_equal(from_sorted$color_agreement, agreement$color_agreement)
})

test_that("plot_alluvial() labels the missing stratum", {
    data <- missing_counts()
    expect_warning(plot_alluvial(data, graphing_columns = cols, column_weights = "value", sorting_algorithm = "greedy"), "already a value")
    expect_no_warning(p <- plot_alluvial(data, graphing_columns = cols, column_weights = "value", sorting_algorithm = "greedy", missing_value = "(missing)"))
    expect_s3_class(p, "ggplot")
    expect_error(plot_alluvial(data, graphing_columns = cols, column_weights = "value", missing_value = NA), "missing_value")
})
