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

test_that("plot_alluvial returns an error with 0 rows", {
    df <- data.frame()
    expect_error(plot_alluvial(df), "no rows")
})

test_that("plot_alluvial returns an error with 1 column", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df), "at least 2 columns")
})

test_that("plot_alluvial returns an error with 2 columns, where column1 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, graphing_columns = c("bad_col", "method2")), "column 'bad_col' is not a column in the dataframe")
})

test_that("plot_alluvial returns an error with 2 columns, where column2 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, graphing_columns = c("method1", "bad_col")), "column 'bad_col' is not a column in the dataframe")
})

test_that("plot_alluvial returns an error with 3 columns, where graphing_columns is NULL", {
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        method3 = sample(1:3, 100, TRUE)
    )
    expect_error(plot_alluvial(df, graphing_columns = NULL), "dataframe has more than 2 columns")
})

test_that("plot_alluvial returns an error with 3 columns, where column2 is not in df", {
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        method3 = sample(1:3, 100, TRUE)
    )
    expect_error(plot_alluvial(df, graphing_columns = c("method1", "bad_col")), "column 'bad_col' is not a column in the dataframe")
})

test_that("plot_alluvial works with a df with 2 columns, both column1 and column2 NULL", {
    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_s3_class(plot_alluvial(df), "ggplot")
})

test_that("plot_alluvial works with 4 columns, column1 and column2 specified", {
    set.seed(42)
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        meta = rnorm(100),
        sample_id = paste0("sample_", 1:100)
    )
    expect_s3_class(plot_alluvial(df, graphing_columns = c("method1", "method2")), "ggplot")
})

test_that("plot_alluvial works with a pre-aggregated df with weight column", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))

    expect_s3_class(plot_alluvial(df, graphing_columns = c("method1", "method2"), column_weights = "value"), "ggplot")
})

test_that("data_sort works with unsorted algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))

    unsorted_df <- data_sort(df, column1 = "method1", column2 = "method2", column_weights = "value", sorting_algorithm = "none")
    # unsorted_df <- dplyr::ungroup(unsorted_df)

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "unsorted_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(unsorted_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)

    expect_equal(unsorted_df, ground_truth_df)
})

test_that("data_sort works with greedy_wolf algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    greedy_wolf_df <- data_sort(df, column1 = "method1", column2 = "method2", column_weights = "value", sorting_algorithm = "greedy_wolf")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "greedy_wolf_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(greedy_wolf_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)

    expect_equal(greedy_wolf_df, ground_truth_df)
})

test_that("data_sort works with greedy_wblf algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))

    greedy_wblf_df <- data_sort(df, column1 = "method1", column2 = "method2", column_weights = "value", sorting_algorithm = "greedy_wblf")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "greedy_wblf_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(greedy_wblf_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)

    expect_equal(greedy_wblf_df, ground_truth_df)
})


test_that("data_sort works with neighbornet algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))

    neighbornet_df <- data_sort(df, column1 = "method1", column2 = "method2", column_weights = "value", sorting_algorithm = "neighbornet")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "neighbornet_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(neighbornet_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)

    expect_equal(neighbornet_df, ground_truth_df)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "none")
    vdiffr::expect_doppelganger("basic_alluvial_plot_unsorted", p)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, greedy_wolf", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "greedy_wolf")
    vdiffr::expect_doppelganger("basic_alluvial_plot_WOLF", p)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, greedy_wblf", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "greedy_wblf")
    vdiffr::expect_doppelganger("basic_alluvial_plot_WBLF", p)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, greedy_NN", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "neighbornet")
    vdiffr::expect_doppelganger("basic_alluvial_plot_NN", p)
})


make_more_neighbornet_2_layer_df <- function() {
    df <- data.frame(
        tissue = c(
            1, 1, 1,
            2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3,
            4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5
        ),
        cluster = c(
            6, 6, 7,
            6, 7, 7, 7, 7, 7,
            6, 8, 8, 8, 8, 8, 8,
            8, 8,
            8, 8, 8, 8, 8, 8, 8, 8, 8
        )
    )
    column1 <- "tissue"
    column2 <- "cluster"
    graphing_columns <- c(column1, column2)

    list(
        df = df,
        graphing_columns = graphing_columns
    )
}

test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 2 layers, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_2_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "none", color_bands = TRUE)

    vdiffr::expect_doppelganger("more_neighbornet_2layer_unsorted", p)
})


test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 2 layers, neighbornet, optimize_column_order_per_cycle FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_2_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order_per_cycle = FALSE)

    vdiffr::expect_doppelganger("more_neighbornet_2layer_neighbornet_oFALSE", p)
})

test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 2 layers, neighbornet, optimize_column_order_per_cycle TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_2_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order_per_cycle = TRUE)

    vdiffr::expect_doppelganger("more_neighbornet_2layer_neighbornet_oTRUE", p)
})




make_more_neighbornet_3_layer_df <- function() {
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
            "male", "male", "male", "male", "male", "male", "male", "male", "male"
        )
    )
    graphing_columns <- c("tissue", "cluster", "sex")

    list(
        df = df,
        graphing_columns = graphing_columns
    )
}

make_more_neighbornet_3_layer_df_with_2_identical_layers <- function() {
    df <- data.frame(
        tissue = c(
            "BRAIN", "BRAIN", "BRAIN",
            "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH",
            "HEART", "HEART", "HEART", "HEART", "HEART", "HEART", "HEART",
            "T CELL", "T CELL",
            "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL"
        ),
        sex = c(
            "male", "female", "male",
            "female", "male", "female", "female", "male", "female",
            "male", "female", "male", "female", "male", "female", "male",
            "female", "male",
            "female", "male", "female", "male", "female", "male", "female", "female", "male"
        ),
        cluster = c(
            "BRAIN", "BRAIN", "BRAIN",
            "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH",
            "HEART", "HEART", "HEART", "HEART", "HEART", "HEART", "HEART",
            "T CELL", "T CELL",
            "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL"
        )
    )
    graphing_columns <- c("tissue", "cluster", "sex")

    list(
        df = df,
        graphing_columns = graphing_columns
    )
}

test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "none", color_bands = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_unsorted", p)
})



test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers, neighbornet, optimize_column_order FALSE, optimize_column_order_per_cycle FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order = FALSE, optimize_column_order_per_cycle = FALSE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_unsorted_oFALSE_ocFALSE", p)
})


test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers, neighbornet, optimize_column_order FALSE, optimize_column_order_per_cycle TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order = FALSE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_unsorted_oFALSE_ocTRUE", p)
})

test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers, neighbornet, optimize_column_order TRUE, optimize_column_order_per_cycle FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_unsorted_oTRUE_ocFALSE", p)
})

test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers, neighbornet, optimize_column_order TRUE, optimize_column_order_per_cycle TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order = TRUE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_unsorted_ooTRUE_ocTRUE", p)
})


test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers with 2 identical layers, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "none", color_bands = TRUE)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_2ident", p)
})

test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers with 2 identical layers, neighbornet, optimize_column_order FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order = FALSE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_with_2ident_oFALSE", p)
})


test_that("VDIFFR - plot_alluvial, more_neighbornet.Rmd, 3 layers with 2 identical layers, neighbornet, optimize_column_order TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_neighbornet_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "neighbornet", color_bands = TRUE, optimize_column_order = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_neighbornet_3layer_with_2_ident_oTRUE", p)
})


test_that("Objective calculation, more_neighbornet.Rmd, 3 layers, unsorted", {
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "none", optimize_column_order = FALSE)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE)

    testthat::expect_equal(num, 225)
})

test_that("Objective calculation, more_neighbornet.Rmd, 3 layers, neighbornet, optimize_column_order FALSE", {
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "neighbornet", optimize_column_order = FALSE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE)

    testthat::expect_equal(num, 30)
})

test_that("Objective calculation, more_neighbornet.Rmd, 3 layers, neighbornet, optimize_column_order TRUE", {
    input <- make_more_neighbornet_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "neighbornet", optimize_column_order = TRUE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE)

    testthat::expect_equal(num, 29)
})



test_that("Objective calculation, more_neighbornet.Rmd, 3 layers with 2 identical layers, unsorted", {
    input <- make_more_neighbornet_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "none", optimize_column_order = FALSE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE)

    testthat::expect_equal(num, 74)
})

test_that("Objective calculation, more_neighbornet.Rmd, 3 layers with 2 identical layers, neighbornet, optimize_column_order FALSE", {
    input <- make_more_neighbornet_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "neighbornet", optimize_column_order = FALSE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE)

    testthat::expect_equal(num, 48)
})

test_that("Objective calculation, more_neighbornet.Rmd, 3 layers with 2 identical layers, neighbornet, optimize_column_order TRUE", {
    input <- make_more_neighbornet_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "neighbornet", optimize_column_order = TRUE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE)

    testthat::expect_equal(num, 48)
})
