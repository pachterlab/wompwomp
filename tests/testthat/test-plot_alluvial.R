test_that("data_sort works with unsorted algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")

    unsorted_df <- data_sort(data, cols = cols, wt = "value", method = "none")
    # unsorted_df <- dplyr::ungroup()(unsorted_df)

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "unsorted_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(unsorted_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(unsorted_df), as.data.frame(ground_truth_df))
})

test_that("data_sort works with greedy_wolf algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    greedy_wolf_df <- data_sort(data, cols = cols, wt = "value", method = "greedy_wolf")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "greedy_wolf_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(greedy_wolf_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(greedy_wolf_df), as.data.frame(ground_truth_df))
})

test_that("data_sort works with greedy_wblf algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    greedy_wblf_df <- data_sort(data, cols = cols, wt = "value", method = "greedy_wblf")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "greedy_wblf_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(greedy_wblf_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(greedy_wblf_df), as.data.frame(ground_truth_df))
})


test_that("data_sort works with tsp algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    tsp_df <- data_sort(data, cols = cols, wt = "value", method = "tsp")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "tsp_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(tsp_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(tsp_df), as.data.frame(ground_truth_df))
})


make_more_tsp_2_layer_df <- function() {
    data <- data.frame(
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
    cols <- c(column1, column2)

    list(
        data = data,
        cols = cols
    )
}

make_more_tsp_3_layer_df <- function() {
    data <- data.frame(
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
    cols <- c("tissue", "cluster", "sex")

    list(
        data = data,
        cols = cols
    )
}

make_more_tsp_3_layer_df_with_2_identical_layers <- function() {
    data <- data.frame(
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
    cols <- c("tissue", "cluster", "sex")

    list(
        data = data,
        cols = cols
    )
}

test_that("Objective calculation, more_tsp.Rmd, 3 layers, unsorted", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- data_preprocess(data = data, cols = cols)

    clus_df_gather_sorted <- data_sort(clus_df_gather, cols = cols, wt = "value", method = "none", column_method = "none")

    num <- determine_crossing_edges(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 316)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers, tsp, optimize_column_order FALSE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- data_preprocess(data = data, cols = cols)

    clus_df_gather_sorted <- data_sort(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "none", weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 153)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers, tsp, optimize_column_order TRUE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- data_preprocess(data = data, cols = cols)

    clus_df_gather_sorted <- data_sort(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "tsp", options = list(optimize_column_order_per_cycle = TRUE), weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 153)
})



test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, unsorted", {
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- data_preprocess(data = data, cols = cols)

    clus_df_gather_sorted <- data_sort(clus_df_gather, cols = cols, wt = "value", method = "none", column_method = "none", weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 74)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order FALSE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- data_preprocess(data = data, cols = cols)

    clus_df_gather_sorted <- data_sort(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "none", weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 95)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order TRUE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- data_preprocess(data = data, cols = cols)

    clus_df_gather_sorted <- data_sort(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "tsp", options = list(optimize_column_order_per_cycle = TRUE), weight_scalar = 1)

    num <- determine_crossing_edges(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 95)
})

test_that("data_color correctly handles multiple factor columns", {
    set.seed(429144)
    
    data <- data.frame(
        method1 = factor(LETTERS[sample(1:3, 100, TRUE)]),
        method2 = factor(LETTERS[27 - sample(1:3, 100, TRUE)])
    )
    
    # sanity check input
    expect_equal(lapply(data, levels), list(
        method1 = c("A", "B", "C"),
        method2 = c("X", "Y", "Z")
    ))
    
    cluster_mapping <- data |>
        data_color(cols = c(method1, method2), method = "left")
    
    # ---- expectations (adjust to actual return type) ----
    
    expect_true(!is.null(cluster_mapping))
    expect_true(length(cluster_mapping) > 0)
    
    # If it's a data.frame
    if (is.data.frame(cluster_mapping)) {
        expect_true(all(c("method1", "method2") %in% names(cluster_mapping)))
    }
    
    # If it's a named list
    if (is.list(cluster_mapping)) {
        expect_true(all(c("method1", "method2") %in% names(cluster_mapping)))
    }
    
    expected <- list(
        method1 = list(
            A = 1L,
            B = 2L,
            C = 3L
        ),
        method2 = list(
            X = 4L,
            Y = 5L,
            Z = 6L
        )
    )
    
    expect_identical(cluster_mapping, expected)
})

test_that("data_color correctly handles multiple factor columns with string column names", {
    set.seed(429144)
    
    data <- data.frame(
        method1 = factor(LETTERS[sample(1:3, 100, TRUE)]),
        method2 = factor(LETTERS[27 - sample(1:3, 100, TRUE)])
    )
    
    # sanity check input
    expect_equal(lapply(data, levels), list(
        method1 = c("A", "B", "C"),
        method2 = c("X", "Y", "Z")
    ))
    
    cluster_mapping <- data |>
        data_color(cols = c("method1", "method2"), method = "left")
    
    # ---- expectations (adjust to actual return type) ----
    
    expect_true(!is.null(cluster_mapping))
    expect_true(length(cluster_mapping) > 0)
    
    # If it's a data.frame
    if (is.data.frame(cluster_mapping)) {
        expect_true(all(c("method1", "method2") %in% names(cluster_mapping)))
    }
    
    # If it's a named list
    if (is.list(cluster_mapping)) {
        expect_true(all(c("method1", "method2") %in% names(cluster_mapping)))
    }
    
    expected <- list(
        method1 = list(
            A = 1L,
            B = 2L,
            C = 3L
        ),
        method2 = list(
            X = 4L,
            Y = 5L,
            Z = 6L
        )
    )
    
    expect_identical(cluster_mapping, expected)
})

test_that("data_color correctly handles multiple factor columns with method advanced", {
    set.seed(429144)
    
    data <- data.frame(
        method1 = factor(LETTERS[sample(1:3, 100, TRUE)]),
        method2 = factor(LETTERS[27 - sample(1:3, 100, TRUE)])
    )
    
    # sanity check input
    expect_equal(lapply(data, levels), list(
        method1 = c("A", "B", "C"),
        method2 = c("X", "Y", "Z")
    ))
    
    cluster_mapping <- data |>
        data_color(cols = c("method1", "method2"), method = "advanced", resolution=10)
    
    # ---- expectations (adjust to actual return type) ----
    
    expect_true(!is.null(cluster_mapping))
    expect_true(length(cluster_mapping) > 0)
    
    # If it's a data.frame
    if (is.data.frame(cluster_mapping)) {
        expect_true(all(c("method1", "method2") %in% names(cluster_mapping)))
    }
    
    # If it's a named list
    if (is.list(cluster_mapping)) {
        expect_true(all(c("method1", "method2") %in% names(cluster_mapping)))
    }
    
    expected <- list(
        method1 = list(
            A = 1L,
            B = 1L,
            C = 2L
        ),
        method2 = list(
            X = 2L,
            Y = 1L,
            Z = 1L
        )
    )
    
    expect_identical(cluster_mapping, expected)
})