#' Setup Python environment for this package
#'
#' This function installs Miniconda if needed, creates a dedicated environment, and installs the required Python packages.
#'
#' @param envname Conda environment name
#' @param packages Packages to install
#' @param use_conda Use conda vs. virtualenv for python environment
#'
#' @return Invisibly returns TRUE if the environment is set up successfully.
#'
#' @examples
#' \dontrun{
#' setup_python_env()
#' }
#'
#' @export
setup_python_env <- function(envname = "wompwomp_env", packages = c(numpy = "numpy==1.23.5", splitspy = "splitspy", pandas = "pandas", scipy = "scipy", leidenalg = "leidenalg", igraph = "python-igraph"), use_conda = TRUE, yes = FALSE) {
    if (use_conda) {
        conda_findable <- tryCatch(file.exists(reticulate::conda_binary()), error = function(e) FALSE)
        if (!conda_findable) {
            if (interactive() || yes) {
                answer <- ifelse(yes, "yes", readline("use_conda is TRUE but Miniconda not found. Would you like to install it now? [y/n]: "))
                if (tolower(answer) %in% c("y", "yes")) {
                    message("Installing Miniconda...")
                    reticulate::install_miniconda()
                } else {
                    stop("Miniconda is required but was not installed. Install with `reticulate::install_miniconda()`, or set use_conda = FALSE.")
                }
            }
        }

        envs <- reticulate::conda_list()$name
        if (!(envname %in% envs)) {
            if (interactive() || yes) {
                answer <- ifelse(yes, "yes", readline(paste0("Conda environment '", envname, "' not found. Create it now? [y/n]: ")))
                if (tolower(answer) %in% c("y", "yes")) {
                    message("Creating conda environment '", envname, "'...")
                    reticulate::conda_create(envname, python_version = "3.10")
                } else {
                    stop("Conda environment '", envname, "' is required but was not created.")
                }
            } else {
                stop(sprintf(
                    "Conda environment '%s' not found, and cannot prompt in non-interactive mode. Create it with: reticulate::conda_create('%s', python_version = \"3.10\")",
                    envname, envname
                ))
            }
        }

        reticulate::use_condaenv(envname, required = TRUE)
    } else {
        env_path <- file.path("~/.virtualenvs", envname)
        if (!dir.exists(path.expand(env_path))) {
            if (interactive() || yes) {
                answer <- ifelse(yes, "yes", readline(paste0("Virtualenv '", envname, "' not found. Create it now? [y/n]: ")))
                if (tolower(answer) %in% c("y", "yes")) {
                    message("Creating virtualenv '", envname, "'...")
                    reticulate::virtualenv_create(envname, python = "python3.10")
                } else {
                    stop("Virtualenv is required but was not created.")
                }
            } else {
                stop(sprintf(
                    "Virtualenv '%s' not found, and cannot prompt in non-interactive mode. Create it with: reticulate::virtualenv_create('%s', python = 'python3.10')",
                    envname, envname
                ))
            }
        }

        reticulate::use_virtualenv(envname, required = TRUE)
    }

    # Only install missing Python modules
    for (import_name in names(packages)) {
        pkg <- packages[[import_name]]
        # Split "pkg==version" if present
        split <- strsplit(pkg, "==", fixed = TRUE)[[1]]
        pypi_name <- split[1]
        required_version <- if (length(split) > 1) split[2] else NULL

        need_install <- FALSE

        if (!reticulate::py_module_available(import_name)) {
            need_install <- TRUE
            message("Module '", pypi_name, "' not found. Will install.")
        } else if (!is.null(required_version)) {
            # Get current installed version
            actual_version <- tryCatch(
                {
                    reticulate::import(import_name)$`__version__`
                },
                error = function(e) NA
            )

            if (is.na(actual_version) || actual_version != required_version) {
                need_install <- TRUE
                message(
                    "Module '", pypi_name, "' has version ", actual_version,
                    " but requires ", required_version, ". Will reinstall."
                )
            }
        }

        if (need_install) {
            if (interactive() || yes) {
                answer <- ifelse(yes, "yes", readline(paste0("Package '", pkg, "' is not installed in conda env '", envname, "'. Install it now with pip? [y/n]: ")))
                if (tolower(answer) %in% c("y", "yes")) {
                    reticulate::conda_install(envname, packages = pkg, pip = TRUE)
                } else {
                    stop("Package '", pkg, "' is required but was not installed.")
                }
            } else {
                stop(
                    "Package '", pkg, "' not found in conda environment '", envname, "', and cannot prompt in non-interactive mode. ",
                    "You can install it manually with:\n",
                    "reticulate::conda_install(\"", envname, "\", packages = \"", pkg, "\", pip = TRUE)"
                )
            }
        }
    }
}

detect_and_setup_python_env <- function(environment = environment, use_conda = use_conda, verbose = FALSE) {
    python_path <- reticulate::py_discover_config(use_environment = environment)$python
    current_env <- basename(dirname(dirname(python_path))) # Works for Conda envs: path typically ends in ".../envs/ENVIRONMENT/bin/python"
    if (current_env != environment) {
        if (verbose) message("Setting up python environment")
        setup_python_env(environment, use_conda = use_conda)
    }
}
