# Skip all tests if Python not available
if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available, skipping reticulate-based tests")
}

# Skip if conda not available
if (inherits(try(reticulate::conda_binary(), silent = TRUE), "try-error")) {
    testthat::skip("Conda not available, skipping reticulate-based tests")
}

# if miniconda not found or environment not found, then prompt setup_python_env
conda_envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
if ((!file.exists(reticulate::conda_binary()) || (!("wompwomp_env" %in% conda_envs)))) {
    if (interactive()) {
        message("Conda env 'wompwomp_env' not found. Run wompwomp::setup_python_env() to create it.")
    }
    testthat::skip("Conda env 'wompwomp_env' not available, skipping reticulate-based tests")
    wompwomp::setup_python_env()
}

# If everything is available, activate
reticulate::use_condaenv("wompwomp_env", required = TRUE)
Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)
