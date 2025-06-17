# if miniconda not found or environment not found, then prompt setup_python_env
if ((!file.exists(reticulate::conda_binary()) || (!("wompwomp_env" %in% reticulate::conda_list()$name)))) {
    wompwomp::setup_python_env()
}
reticulate::use_condaenv("wompwomp_env", required = TRUE)
Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)
