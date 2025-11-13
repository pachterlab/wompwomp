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


test_that("data_sort works with tsp algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))

    tsp_df <- data_sort(df, column1 = "method1", column2 = "method2", column_weights = "value", sorting_algorithm = "tsp")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "tsp_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(tsp_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)

    expect_equal(tsp_df, ground_truth_df)
})


make_more_tsp_2_layer_df <- function() {
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

make_more_tsp_3_layer_df <- function() {
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

make_more_tsp_3_layer_df_with_2_identical_layers <- function() {
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

test_that("Objective calculation, more_tsp.Rmd, 3 layers, unsorted", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "none", optimize_column_order = FALSE)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns)$output_objective

    testthat::expect_equal(num, 316)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers, tsp, optimize_column_order FALSE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "tsp", optimize_column_order = FALSE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns)$output_objective

    testthat::expect_equal(num, 153)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers, tsp, optimize_column_order TRUE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "tsp", optimize_column_order = TRUE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns)$output_objective

    testthat::expect_equal(num, 153)
})



test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, unsorted", {
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "none", optimize_column_order = FALSE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns)$output_objective

    testthat::expect_equal(num, 74)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order FALSE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "tsp", optimize_column_order = FALSE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns)$output_objective

    testthat::expect_equal(num, 95)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order TRUE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns)

    clus_df_gather_sorted <- data_sort(clus_df_gather, graphing_columns = graphing_columns, column_weights = "value", sorting_algorithm = "tsp", optimize_column_order = TRUE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns)$output_objective

    testthat::expect_equal(num, 95)
})
