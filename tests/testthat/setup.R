# if miniconda not found or environment not found, then prompt setup_python_env
if ((!file.exists(reticulate::conda_binary()) || (!("alluvialmatch_env" %in% reticulate::conda_list()$name)))) {
    alluvialmatch::setup_python_env()
}
reticulate::use_condaenv("alluvialmatch_env", required = TRUE)
Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)

