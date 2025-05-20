test_that("CLI plot matches ground truth PNG", {
    devtools::load_all()

    # Paths
    cli_cmd_path <- file.path(here::here(), "exec", "alluvialmatch")
    command <- "plot_alluvial"
    df_path <- file.path(here::here(), "tmp_files", "df_tmp.csv")
    output_path <- tempfile(fileext = ".png")
    # output_path <- file.path(here::here(), "tmp_files", "tmp_image11.png")
    truth_path <- file.path(here::here(), "tests", "testthat", "ground_truth", "sorting_none.png")

    # Run CLI
    browser()
    # cat(cli_cmd_path, c(command, "--df", df_path, "--output_plot_path", output_path))
    plot_alluvial(df=df_path, output_plot_path=output_path, sorting_algorithm="None")
    system2(cli_cmd_path, c(command, "--df", df_path, "--output_plot_path", output_path, "--sorting_algorithm", "None"))

    # Load and compare images
    img_new   <- magick::image_read(output_path)
    img_truth <- magick::image_read(truth_path)

    # Compare (you can also set fuzz if small differences are tolerable)
    diff <- magick::image_compare(img_new, img_truth, metric = "AE")  # Absolute Error
    diff_val <- as.numeric(attr(diff, "distortion"))

    expect_true(diff_val == 0, info = sprintf("Images differ (distortion = %f)", diff_val))
})
