test_that("sort_to_uncross works with unsorted algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")

    unsorted_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "none")
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

test_that("sort_to_uncross works with greedy algorithm and a fixed column", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    greedy_wolf_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "greedy", fixed_column = "method1")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "greedy_wolf_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(greedy_wolf_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(greedy_wolf_df), as.data.frame(ground_truth_df))
})

test_that("sort_to_uncross works with greedy algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    greedy_wblf_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "greedy")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "greedy_wblf_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(greedy_wblf_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(greedy_wblf_df), as.data.frame(ground_truth_df))
})


test_that("sort_to_uncross works with barycenter algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    barycenter_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "barycenter")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "barycenter_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(barycenter_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(barycenter_df), as.data.frame(ground_truth_df))
})

test_that("sort_to_uncross works with median algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    median_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "median")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "median_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(median_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(median_df), as.data.frame(ground_truth_df))
})

test_that("sort_to_uncross works with barycenter algorithm and a fixed column", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    barycenter_one_sided_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "barycenter", fixed_column = "method1")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "barycenter_one_sided_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(barycenter_one_sided_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(barycenter_one_sided_df), as.data.frame(ground_truth_df))
})

test_that("sort_to_uncross works with median algorithm and a fixed column", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    median_one_sided_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "median", fixed_column = "method1")

    ground_truth_df_path <- normalizePath(testthat::test_path("ground_truth", "median_one_sided_df.rds"))

    if (!file.exists(ground_truth_df_path)) {
        saveRDS(median_one_sided_df, file = ground_truth_df_path)
    }

    ground_truth_df <- readRDS(ground_truth_df_path)
    ground_truth_df <- ground_truth_df[, c(cols, "value"), drop = FALSE]
    ground_truth_df <- ground_truth_df |> dplyr::ungroup()

    expect_equal(as.data.frame(median_one_sided_df), as.data.frame(ground_truth_df))
})

sweep_test_df <- function(n_cols, seed = 42) {
    set.seed(seed)
    latent <- sample(1:4, 400, TRUE)
    raw_df <- as.data.frame(lapply(setNames(seq_len(n_cols), paste0("method", seq_len(n_cols))), function(i) {
        noisy <- ifelse(runif(400) < 0.8, latent, sample(1:4, 400, TRUE))
        # scramble the labels so the default alphabetical order is far from the clustered one
        sample(LETTERS[1:4])[noisy]
    }))
    list(
        data = as.data.frame(dplyr::count(raw_df, dplyr::across(dplyr::everything()), name = "value")),
        cols = paste0("method", seq_len(n_cols))
    )
}

test_that("fixed_column keeps the order of every fixed column, for every method and number of axes", {
    for (n_cols in 2:4) {
        input <- sweep_test_df(n_cols)
        cols <- input$cols
        unsorted_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = "none", column_method = "none")
        fixed_sets <- list(cols[1], cols[n_cols], cols[c(1, n_cols)])
        if (n_cols > 2) fixed_sets <- c(fixed_sets, list(cols[2]))
        for (m in c("greedy", "barycenter", "median", "neighbornet", "tsp", "random")) {
            for (fixed in fixed_sets) {
                set.seed(1)
                sorted_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = m, column_method = "none", fixed_column = fixed)
                for (col in fixed) {
                    expect_equal(levels(sorted_df[[col]]), levels(unsorted_df[[col]]), info = sprintf("method %s, %d axes, fixed %s", m, n_cols, paste(fixed, collapse = "+")))
                }
            }
        }
    }
})

test_that("fixed_column accepts positions in cols", {
    input <- sweep_test_df(3)
    by_name <- sort_to_uncross(input$data, cols = input$cols, wt = "value", method = "barycenter", column_method = "none", fixed_column = c("method1", "method3"))
    by_position <- sort_to_uncross(input$data, cols = input$cols, wt = "value", method = "barycenter", column_method = "none", fixed_column = c(1, 3))
    expect_identical(by_name, by_position)
})

test_that("fixed_column must name or index cols", {
    input <- sweep_test_df(3)
    expect_error(sort_to_uncross(input$data, cols = input$cols, wt = "value", method = "greedy", fixed_column = "value"), "not in cols")
    expect_error(sort_to_uncross(input$data, cols = input$cols, wt = "value", method = "greedy", fixed_column = 4), "not positions in cols")
})

test_that("fixing every column leaves the strata unsorted", {
    input <- sweep_test_df(3)
    unsorted_df <- sort_to_uncross(input$data, cols = input$cols, wt = "value", method = "none", column_method = "none")
    for (m in c("greedy", "barycenter", "median")) {
        sorted_df <- sort_to_uncross(input$data, cols = input$cols, wt = "value", method = m, column_method = "none", fixed_column = input$cols)
        expect_identical(lapply(sorted_df[input$cols], levels), lapply(unsorted_df[input$cols], levels))
    }
})

test_that("greedy/barycenter/median sort any number of axes", {
    for (n_cols in 3:4) {
        input <- sweep_test_df(n_cols)
        cols <- input$cols
        unsorted_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = "none", column_method = "none")
        unsorted_objective <- compute_crossing_objective(unsorted_df, cols = cols, wt = "value")$output_objective
        for (m in c("greedy", "barycenter", "median")) {
            for (fixed in list(NULL, cols[2])) {
                sorted_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = m, column_method = "none", fixed_column = fixed)
                expect_identical(names(sorted_df)[seq_along(cols)], cols)
                objective <- compute_crossing_objective(sorted_df, cols = cols, wt = "value")$output_objective
                expect_lt(objective, unsorted_objective)
            }
        }
    }
})

test_that("greedy/barycenter/median honour column_method with more than two axes", {
    input <- sweep_test_df(2)
    data <- input$data
    data$copy <- data$method1
    cols <- c("method1", "method2", "copy")
    for (m in c("greedy", "barycenter", "median")) {
        set.seed(42)
        sorted_df <- sort_to_uncross(data, cols = cols, wt = "value", method = m, column_method = "tsp")
        axis_order <- names(sorted_df)[1:3]
        expect_setequal(axis_order, cols)
        # the two identical layers belong next to each other
        expect_equal(abs(diff(match(c("method1", "copy"), axis_order))), 1)
        expect_equal(compute_crossing_objective(sorted_df, cols = axis_order, wt = "value")$output_objective,
                     compute_crossing_objective(sorted_df, cols = c("method2", "method1"), wt = "value")$output_objective)
    }
})

test_that("sweep passes reproduce the two-axis schedules and propagate away from fixed axes", {
    expect_equal(sweep_passes(2, integer(0)), list(c(reordered = 2, stable = 1), c(reordered = 1, stable = 2)))
    expect_equal(sweep_passes(2, 1), list(c(reordered = 2, stable = 1)))
    expect_equal(sweep_passes(2, 2), list(c(reordered = 1, stable = 2)))
    expect_equal(sweep_passes(4, 2), list(c(reordered = 3, stable = 2), c(reordered = 4, stable = 3), c(reordered = 1, stable = 2)))
})

test_that("random_initializations never does worse than a single initialization", {
    input <- sweep_test_df(3)
    for (m in c("greedy", "barycenter", "median")) {
        single <- sort_to_uncross(input$data, cols = input$cols, wt = "value", method = m, column_method = "none")
        set.seed(7)
        several <- sort_to_uncross(input$data, cols = input$cols, wt = "value", method = m, column_method = "none", options = list(random_initializations = 5))
        expect_lte(
            compute_crossing_objective(several, cols = input$cols, wt = "value")$output_objective,
            compute_crossing_objective(single, cols = input$cols, wt = "value")$output_objective
        )
    }
})

test_that("deprecated method names warn and map onto the new ones", {
    input <- sweep_test_df(2)
    cols <- input$cols
    deprecated <- list(
        greedy_wolf = list(method = "greedy", fixed_column = "method1"),
        greedy_wblf = list(method = "greedy", fixed_column = NULL),
        barycenter_one_sided = list(method = "barycenter", fixed_column = "method1"),
        median_one_sided = list(method = "median", fixed_column = "method1")
    )
    for (old in names(deprecated)) {
        expect_warning(old_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = old), "deprecated")
        new_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = deprecated[[old]]$method, fixed_column = deprecated[[old]]$fixed_column)
        expect_identical(old_df, new_df)
    }
    expect_warning(old_df <- sort_to_uncross(input$data, cols = cols, wt = "value", method = "greedy_wolf", fixed_column = "method2"), "deprecated")
    expect_identical(old_df, sort_to_uncross(input$data, cols = cols, wt = "value", method = "greedy", fixed_column = "method2"))
})

test_that("sort_to_uncross works with tsp algorithm", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    data <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))
    cols = c("method1", "method2")
    tsp_df <- sort_to_uncross(data, cols = cols, wt = "value", method = "tsp")

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

    clus_df_gather <- prep_for_lodes(data = data, cols = cols)

    clus_df_gather_sorted <- sort_to_uncross(clus_df_gather, cols = cols, wt = "value", method = "none", column_method = "none")

    num <- compute_crossing_objective(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 225)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers, tsp, optimize_column_order FALSE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- prep_for_lodes(data = data, cols = cols)

    clus_df_gather_sorted <- sort_to_uncross(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "none", alpha = 1e6, options = list(weight_scalar = 1))

    num <- compute_crossing_objective(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 44) # was 57 before the nearest-right/nearest-left within-stratum ordering
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers, tsp, optimize_column_order TRUE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- prep_for_lodes(data = data, cols = cols)

    clus_df_gather_sorted <- sort_to_uncross(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "tsp", alpha = 1e6, options = list(optimize_column_order_per_cycle = TRUE, weight_scalar = 1))

    num <- compute_crossing_objective(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 44) # was 57 before the nearest-right/nearest-left within-stratum ordering
})



test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, unsorted", {
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- prep_for_lodes(data = data, cols = cols)

    clus_df_gather_sorted <- sort_to_uncross(clus_df_gather, cols = cols, wt = "value", method = "none", column_method = "none", alpha = 1e6, options = list(weight_scalar = 1))

    num <- compute_crossing_objective(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 74)
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order FALSE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- prep_for_lodes(data = data, cols = cols)

    clus_df_gather_sorted <- sort_to_uncross(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "none", alpha = 1e6, options = list(weight_scalar = 1))

    num <- compute_crossing_objective(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 50) # was 56 before the nearest-right/nearest-left within-stratum ordering
})

test_that("Objective calculation, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order TRUE", {
    set.seed(42)
    
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    data <- input$data
    cols <- input$cols

    clus_df_gather <- prep_for_lodes(data = data, cols = cols)

    clus_df_gather_sorted <- sort_to_uncross(clus_df_gather, cols = cols, wt = "value", method = "tsp", column_method = "tsp", alpha = 1e6, options = list(optimize_column_order_per_cycle = TRUE, weight_scalar = 1))

    num <- compute_crossing_objective(clus_df_gather_sorted, cols = cols)$output_objective

    testthat::expect_equal(num, 50) # was 56 before the nearest-right/nearest-left within-stratum ordering
})

test_that("get_lode_clusters correctly handles multiple factor columns", {
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
        get_lode_clusters(cols = c(method1, method2), method = "left")
    
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

test_that("get_lode_clusters correctly handles multiple factor columns with string column names", {
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
        get_lode_clusters(cols = c("method1", "method2"), method = "left")
    
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

test_that("get_lode_clusters correctly handles multiple factor columns with method advanced", {
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
    
    # `resolution` is a modularity resolution (find_colors_advanced() passes
    # objective_function = "modularity" to igraph::cluster_leiden(), whose own
    # default is CPM). At resolution 1 the two layers' blocks merge into the
    # two communities the data actually contains; raising it splits every block
    # into its own community.
    cluster_mapping <- data |>
        get_lode_clusters(cols = c("method1", "method2"), method = "advanced", resolution = 1)
    
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

    # A high resolution gives every block its own colour.
    cluster_mapping_fine <- data |>
        get_lode_clusters(cols = c("method1", "method2"), method = "advanced", resolution = 10)
    expect_identical(
        cluster_mapping_fine,
        list(
            method1 = list(A = 1L, B = 2L, C = 3L),
            method2 = list(X = 4L, Y = 5L, Z = 6L)
        )
    )

    # M: total weight of observations coloured the same on both sides.
    agreement <- compute_color_agreement(
        data, cols = c("method1", "method2"), mapping = cluster_mapping
    )
    expect_equal(agreement$color_agreement, 63)
    expect_equal(nrow(agreement$per_pair), 1L)
    # Every block its own colour => nothing agrees across layers.
    expect_equal(
        compute_color_agreement(
            data, cols = c("method1", "method2"), mapping = cluster_mapping_fine
        )$color_agreement,
        0
    )
})