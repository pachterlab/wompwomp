test_that("plot_alluvial returns a ggplot object", {
    set.seed(42)
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    p <- plot_alluvial(df)
    expect_s3_class(p, "ggplot")
})

test_that("plot_alluvial returns a ggplot object", {
    set.seed(42)
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    p <- plot_alluvial(df)
    expect_s3_class(p, "ggplot")
})

test_that("plot_alluvial returns an error with 0 columns", {
    df <- data.frame()
    expect_error(plot_alluvial(df), "at least two columns")
})

test_that("plot_alluvial returns an error with 1 column", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df), "at least two columns")
})

test_that("plot_alluvial returns an error with 2 columns, where column1 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, column1 = "bad_col", column2 = "method2"), "column1.*not a column")
})

test_that("plot_alluvial returns an error with 2 columns, where column2 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, column1 = "method1", column2 = "bad_col"), "column2.*not a column")
})

test_that("plot_alluvial returns an error with 3 columns, where column1 is NULL", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE),
                     method2 = sample(1:3, 100, TRUE),
                     method3 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, column1 = NULL, column2 = "method2"), "Dataframe has more than two columns")
})

test_that("plot_alluvial returns an error with 3 columns, where column2 is NULL", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE),
                     method2 = sample(1:3, 100, TRUE),
                     method3 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, column1 = "method1", column2 = NULL), "Dataframe has more than two columns")
})

test_that("plot_alluvial returns an error with 3 columns, where column1 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE),
                     method2 = sample(1:3, 100, TRUE),
                     method3 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, column1 = "bad_col", column2 = "method2"), "column1.*not a column")
})

test_that("plot_alluvial returns an error with 3 columns, where column2 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE),
                     method2 = sample(1:3, 100, TRUE),
                     method3 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, column1 = "method1", column2 = "bad_col"), "column2.*not a column")
})

test_that("plot_alluvial works with a df with 2 columns, both column1 and column2 NULL", {
    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_s3_class(plot_alluvial(df), "ggplot")
})

test_that("plot_alluvial works with 2 columns, column1 specified, column2 NULL", {
    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_s3_class(plot_alluvial(df, column1 = "method1"), "ggplot")
})

test_that("plot_alluvial works with 2 columns, column2 specified, column1 NULL", {
    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_s3_class(plot_alluvial(df, column2 = "method2"), "ggplot")
})

test_that("VDIFFR - plot_alluvial works with 2 columns, column2 specified, column1 NULL", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, column2 = "method2")
    vdiffr::expect_doppelganger("basic alluvial plot", p)
})

test_that("plot_alluvial works with 4 columns, column1 and column2 specified", {
    set.seed(42)
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        meta = rnorm(100),
        sample_id = paste0("sample_", 1:100)
    )
    expect_s3_class(plot_alluvial(df, column1 = "method1", column2 = "method2"), "ggplot")
})

test_that("plot_alluvial works with a pre-aggregated df with weight column", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "weight"))

    expect_s3_class(plot_alluvial(df, column1 = "method1", column2 = "method2", column_weights = "weight"), "ggplot")
})

test_that("VDIFFR - plot_alluvial works with a pre-aggregated df with weight column", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "weight"))

    p <- plot_alluvial(df, column1 = "method1", column2 = "method2", column_weights = "weight")

    vdiffr::expect_doppelganger("basic alluvial plot weighted", p)
})

