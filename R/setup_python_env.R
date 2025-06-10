#' Setup Python environment for this package
#'
#' This function installs Miniconda if needed, creates a dedicated environment, and installs the required Python packages.
#'
#' @param envname Conda environment name
#' @param packages Packages to install
#'
#' @export
setup_python_env <- function(envname = "alluvialmatch_env", packages = c("numpy==1.23.5", "splitspy")) {
    conda_findable <- tryCatch(file.exists(reticulate::conda_binary()), error = function(e) FALSE)
    if (!conda_findable) {
        message("Installing Miniconda...")
        reticulate::install_miniconda()
    }

    envs <- reticulate::conda_list()$name
    if (!(envname %in% envs)) {
        message("Creating conda environment '", envname, "'...")
        reticulate::conda_create(envname, python_version = "3.10")
    }

    reticulate::use_condaenv(envname, required = TRUE)

    # Only install missing Python modules
    for (pkg in packages) {
        # Split "pkg==version" if present
        split <- strsplit(pkg, "==", fixed = TRUE)[[1]]
        pkg_name <- split[1]
        required_version <- if (length(split) > 1) split[2] else NULL

        need_install <- FALSE

        if (!reticulate::py_module_available(pkg_name)) {
            need_install <- TRUE
            message("Module '", pkg_name, "' not found. Will install.")
        } else if (!is.null(required_version)) {
            # Get current installed version
            actual_version <- tryCatch({
                reticulate::import(pkg_name)$`__version__`
            }, error = function(e) NA)

            if (is.na(actual_version) || actual_version != required_version) {
                need_install <- TRUE
                message("Module '", pkg_name, "' has version ", actual_version,
                        " but requires ", required_version, ". Will reinstall.")
            }
        }

        if (need_install) {
            reticulate::conda_install(envname, packages = pkg, pip = TRUE)
        }
    }
}
