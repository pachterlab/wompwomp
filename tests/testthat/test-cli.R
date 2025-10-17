skip("Skipping this test file during development or check")

# dev mode for testing, not for compiling
dev <- !("wompwomp" %in% installed.packages())  # !requireNamespace("wompwomp", quietly = TRUE)

cli_cmd_path <- system.file("exec", "wompwomp", package = "wompwomp")
if (cli_cmd_path == "") {  # file not installed
    cli_cmd_path <- file.path(here::here(), "exec", "wompwomp")  # cli_cmd_path <- testthat::test_path("..", "..", "exec", "wompwomp")
}

type_checking_files <- function(output_path, truth_path, check = TRUE) {
    # check for file path type
    if (!is.character(output_path) || length(output_path) != 1) {
        stop("`output_path` must be a character string representing a file path.")
    }
    if (!is.character(truth_path) || length(truth_path) != 1) {
        stop("`truth_path` must be a character string representing a file path.")
    }

    # Check that output file exists
    if (!file.exists(output_path)) {
        stop(sprintf("Output file not found at: %s", output_path))
    }

    # If ground truth doesn't exist, initialize it
    if (!file.exists(truth_path)) {
        dir.create(dirname(truth_path), recursive = TRUE, showWarnings = FALSE)
        file.copy(output_path, truth_path)
        message(sprintf("Ground truth image initialized at: %s", truth_path))
        return(invisible(TRUE)) # No comparison needed on first run
    }

    if (!check) {
        return(invisible(TRUE))
    } # Allow skipping comparison
}

compare_images <- function(output_path, truth_path, tolerance_fraction = 0.01, check = TRUE) {
    # tolerance_fraction of 0.01 means up to 1% difference allowable
    type_checking_files(output_path = output_path, truth_path = truth_path, check = check)

    # Load and compare images
    img_new <- magick::image_read(output_path)
    img_truth <- magick::image_read(truth_path)

    # Ensure same size before comparing
    stopifnot(identical(
        magick::image_info(img_new)[, c("width", "height")],
        magick::image_info(img_truth)[, c("width", "height")]
    ))

    # Total pixel count
    dims <- magick::image_info(img_new)
    total_pixels <- dims$width * dims$height

    # Compare images using Absolute Error metric
    diff <- magick::image_compare(img_new, img_truth, metric = "AE")
    diff_val <- as.numeric(attr(diff, "distortion"))

    tolerance <- tolerance_fraction * total_pixels

    testthat::expect_true(
        diff_val <= tolerance,
        info = sprintf("Images differ (distortion = %f)", diff_val)
    )

    invisible(TRUE)
}

compare_csvs <- function(output_path, truth_path, check = TRUE) {
    type_checking_files(output_path = output_path, truth_path = truth_path, check = check)

    df_output <- read.csv(output_path, stringsAsFactors = FALSE)
    df_truth <- read.csv(truth_path, stringsAsFactors = FALSE)

    testthat::expect_equal(df_output, df_truth)

    invisible(TRUE)
}

# ./exec/wompwomp plot_alluvial --df tests/testthat/ground_truth/df_tests_cli.csv --output_plot_path tests/testthat/ground_truth/sorting_none.png --sorting_algorithm none
test_that("CLI plot_alluvial, no sort", {
    # Paths
    command <- "plot_alluvial"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".png")
    # output_path <- file.path(here::here(), "tmp_files", "tmp_image21.png")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_none.png"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_plot_path", output_path,
        "--sorting_algorithm", "none",
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_images(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI plot_alluvial, WOLF left fixed", {
    # Paths
    command <- "plot_alluvial"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".png")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wolf_left.png"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_plot_path", output_path,
        "--sorting_algorithm", "greedy_wolf",
        "--fixed_column", 1,
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_images(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI plot_alluvial, WOLF right fixed", {
    # Paths
    command <- "plot_alluvial"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".png")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wolf_right.png"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_plot_path", output_path,
        "--sorting_algorithm", "greedy_wolf",
        "--fixed_column", 2,
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_images(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI plot_alluvial, WBLF", {
    # Paths
    command <- "plot_alluvial"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".png")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wblf.png"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_plot_path", output_path,
        "--sorting_algorithm", "greedy_wblf",
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_images(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI plot_alluvial, TSP", {
    # Paths
    set.seed(42)
    command <- "plot_alluvial"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".png")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_tsp.png"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_plot_path", output_path,
        "--sorting_algorithm", "tsp",
        "--weight_scalar", "10",
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_images(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI data_sort, WOLF left fixed", {
    # Paths
    command <- "data_sort"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".csv")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wolf_left_df.csv"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_df_path", output_path,
        "--sorting_algorithm", "greedy_wolf",
        "--fixed_column", 1,
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_csvs(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI data_sort, WOLF right fixed", {
    # Paths
    command <- "data_sort"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".csv")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wolf_right_df.csv"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_df_path", output_path,
        "--sorting_algorithm", "greedy_wolf",
        "--fixed_column", 2,
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_csvs(output_path = output_path, truth_path = truth_path, check = TRUE)
})

# ./exec/wompwomp data_sort --df tests/testthat/ground_truth/df_tests_cli.csv --output_df_path tests/testthat/ground_truth/sorting_wblf_df.csv --sorting_algorithm greedy_wblf --graphing_columns tissue leiden
test_that("CLI data_sort, WBLF", {
    # Paths
    command <- "data_sort"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".csv")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wblf_df.csv"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_df_path", output_path,
        "--sorting_algorithm", "greedy_wblf",
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_csvs(output_path = output_path, truth_path = truth_path, check = TRUE)
})

test_that("CLI data_sort, TSP", {
    # Paths
    set.seed(42)
    command <- "data_sort"
    df_path <- normalizePath(testthat::test_path("ground_truth", "df_tests_cli.csv"))
    output_path <- tempfile(fileext = ".csv")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "sorting_tsp_df.csv"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_df_path", output_path,
        "--sorting_algorithm", "tsp",
        "--weight_scalar", "10",
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_csvs(output_path = output_path, truth_path = truth_path, check = TRUE)
})

# ./exec/wompwomp determine_crossing_edges --df tests/testthat/ground_truth/sorting_wblf_df.csv --output_df_path tests/testthat/ground_truth/crossing_wblf_df.csv --column1 col1_int --column2 col2_int
test_that("CLI determine_crossing_edges", {
    # Paths
    command <- "determine_crossing_edges"
    df_path <- normalizePath(testthat::test_path("ground_truth", "sorting_wblf_df.csv"))
    output_path <- tempfile(fileext = ".csv")
    truth_path <- normalizePath(testthat::test_path("ground_truth", "crossing_wblf_df.csv"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path,
        "--output_df_path", output_path,
        "--column1", "col1_int",
        "--column2", "col2_int",
        "--quiet", "TRUE"
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    system2(cli_cmd_path, args)

    compare_csvs(output_path = output_path, truth_path = truth_path, check = TRUE)
})

# ./exec/wompwomp determine_weighted_layer_free_objective --df tests/testthat/ground_truth/crossing_wblf_df.csv
test_that("CLI determine_weighted_layer_free_objective", {
    # Paths
    command <- "determine_weighted_layer_free_objective"
    df_path <- normalizePath(testthat::test_path("ground_truth", "crossing_wblf_df.csv"))

    # Run CLI
    args <- c(
        command,
        "--df", df_path
    )
    if (dev) {
        args <- c(args, "--dev")
    }
    # cat(cli_cmd_path, paste(args, collapse = " "), "\n")
    output <- system2(cli_cmd_path, args, stdout = TRUE)
    num <- as.integer(gsub(".*\\s", "", output))

    testthat::expect_equal(num, 14)
})
