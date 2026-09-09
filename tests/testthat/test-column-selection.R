# `cols` and `wt` are documented as character vectors, but they are captured
# with `enquo()` so that bare column names also work. A programmatic caller
# (e.g. ggalluvial's `sort_strata`/`color_strata` parameters) has to pass them
# as variables holding character vectors, which tidyselect deprecated selecting
# on. `as_name_selection()` rewrites those to `all_of(<value>)`.

toy_counts <- function() {
    set.seed(42)
    raw_df <- data.frame(
        method1 = sample(c("a1", "a2", "a3"), 200, TRUE),
        method2 = sample(c("b1", "b2", "b3", "b4"), 200, TRUE)
    )
    as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
}

test_that("as_name_selection() rewrites only symbols bound to character", {
    data <- data.frame(method1 = 1, method2 = 2, other = 3)
    cols <- c("method1", "method2")
    quo <- as_name_selection(rlang::quo(cols))
    expect_true(rlang::is_call(rlang::quo_get_expr(quo)))
    expect_no_warning(sel <- tidyselect::eval_select(quo, data))
    expect_identical(names(sel), c("method1", "method2"))

    # a literal, a tidyselect call and a bare column name are left alone
    for (e in list(rlang::quo("method1"),
                   rlang::quo(c(method1, method2)),
                   rlang::quo(method1),
                   rlang::quo(dplyr::starts_with("method")))) {
        expect_identical(rlang::quo_get_expr(as_name_selection(e)),
                         rlang::quo_get_expr(e))
    }

    # a symbol bound to something other than column names is left alone
    n <- 3
    expect_identical(rlang::quo_get_expr(as_name_selection(rlang::quo(n))),
                     rlang::quo_get_expr(rlang::quo(n)))

    # a variable holding no column becomes an explicit NULL selection, which
    # `eval_select()` reads as "no columns" rather than an external vector
    empty <- NULL
    expect_null(rlang::quo_get_expr(as_name_selection(rlang::quo(empty))))
    expect_no_warning(sel <- tidyselect::eval_select(
        as_name_selection(rlang::quo(empty)), data
    ))
    expect_length(sel, 0L)

    expect_identical(as_name_selection(rlang::quo(NULL)), rlang::quo(NULL))
})

test_that("get_lode_clusters() accepts cols and wt as character variables", {
    data <- toy_counts()
    cols <- c("method1", "method2")
    wt <- "value"

    expect_no_warning(from_vars <- get_lode_clusters(data, cols = cols, wt = wt))
    from_literal <- get_lode_clusters(data, cols = c("method1", "method2"),
                                      wt = "value")
    from_symbols <- get_lode_clusters(data, cols = c(method1, method2),
                                      wt = value)

    expect_identical(from_vars, from_literal)
    expect_identical(from_vars, from_symbols)
    expect_named(from_vars, cols)
})

test_that("get_lode_clusters() still accepts a NULL weight", {
    data <- data.frame(
        method1 = c("a1", "a1", "a2", "a2"),
        method2 = c("b1", "b2", "b1", "b2")
    )
    cols <- c("method1", "method2")
    expect_no_warning(mapping <- get_lode_clusters(data, cols = cols))
    expect_named(mapping, cols)
})

test_that("sort_to_uncross() accepts cols and wt as character variables", {
    data <- toy_counts()
    cols <- c("method1", "method2")
    wt <- "value"

    expect_no_warning(from_vars <- sort_to_uncross(
        data, cols = cols, wt = wt, method = "tsp", column_method = "none"
    ))
    from_literal <- sort_to_uncross(
        data, cols = c("method1", "method2"), wt = "value",
        method = "tsp", column_method = "none"
    )
    expect_identical(lapply(from_vars[cols], levels),
                     lapply(from_literal[cols], levels))
})

test_that("lode_cluster_pal() accepts cols as a character variable", {
    data <- toy_counts()
    cols <- c("method1", "method2")
    sorted <- sort_to_uncross(data, cols = cols, wt = "value", method = "tsp",
                              column_method = "none")
    mapping <- get_lode_clusters(sorted, cols = cols, wt = "value")
    pal <- lode_cluster_pal(sorted, cols = cols, mapping = mapping)
    expect_true(all(unlist(lapply(cols, function(c) levels(sorted[[c]]))) %in%
                        names(pal)))
})
