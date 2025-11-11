#' wompwomp: Cluster-matching alluvial plots
#'
#' Utils
#' @docType package
#' @name wompwomp
#'
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%

utils::globalVariables(c(
    ".data", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "total", "cum_y", "best_cluster_agreement"
))

print_function_params <- function(display_df = FALSE) {
    f <- sys.function(sys.parent())
    call <- sys.call(sys.parent())
    defaults <- as.list(formals(f))
    call_args <- as.list(call)[-1]
    
    # Evaluate args, but wrap in list() so NULL survives
    eval_args <- lapply(call_args, function(arg) {
        if (is.symbol(arg) && as.character(arg) == "NULL") {
            list(NULL)  # wrapper preserves NULL
        } else {
            list(eval(arg, envir = parent.frame()))
        }
    })
    
    formal_names <- names(defaults)
    unnamed <- which(!nzchar(names(call_args)))
    names(call_args)[unnamed] <- formal_names[unnamed]
    
    # Flatten one level
    names(eval_args) <- names(call_args)
    eval_args <- lapply(eval_args, `[[`, 1)
    
    # Manual merge (defaults first, then override)
    all_args <- defaults
    for (nm in names(eval_args)) {
        all_args[nm] <- list(eval_args[[nm]])  # assign inside list()
    }
    
    # Print
    for (nm in names(all_args)) {
        val <- all_args[[nm]]
        if (is.null(val)) {
            message(nm, " = NULL")
        } else if (is.data.frame(val) && !display_df) {
            message(nm, " = dataframe")
        } else if (length(val) == 1) {
            message(nm, " = ", val)
        } else {
            message(nm, " =")
            for (j in seq_along(val)) {
                elt_name <- names(val)[j]
                if (!is.null(elt_name) && nzchar(elt_name)) {
                    message("  ", elt_name, " : ", val[[j]])
                } else {
                    message("  [", j, "] : ", val[[j]])
                }
            }
        }
    }
}

lowercase_args <- function(arg_names) {
    for (nm in arg_names) {
        val <- get(nm, envir = parent.frame())
        if (is.character(val)) {
            assign(nm, tolower(val), envir = parent.frame())
        }
    }
}


check_python_setup_with_necessary_packages <- function(necessary_packages_for_this_step = NULL, additional_message = "", environment = "wompwomp_env", use_conda = TRUE) {
    ### make sure that necessary_packages_for_this_step uses the IMPORT package name, not the pypi package name
    
    # Skip check if script was run from command line (including checking from build/check) - this is ok because I set up my python environment in exec/wompwomp now
    if (identical(Sys.getenv("R_SCRIPT_FROM_CLI"), "true")) {
        return(invisible(NULL))
    }
    
    # detect_and_setup_python_env(environment = environment, use_conda = use_conda)  #!!! uncomment later if I want python to be set up upon function call
    
    # can comment out relevant if I call wompwomp::setup_python_env() in here (above)
    if (!reticulate::py_available(initialize = FALSE)) {
        if (is.null(additional_message)) {
            stop("Python environment is not set up.")  #  Please run wompwomp::setup_python_env().
        } else {
            stop(sprintf(
                "Python environment is not set up. %s.",  # Please run wompwomp::setup_python_env(), or %s.
                additional_message
            ))
        }
    }
    if (!is.null(necessary_packages_for_this_step)) {
        for (package in necessary_packages_for_this_step) {
            if (!reticulate::py_module_available(package)) {
                if (is.null(additional_message)) {
                    stop(sprintf(
                        "Python module '%s' is not available.",  # Please run wompwomp::setup_python_env().
                        package
                    ))
                } else {
                    stop(sprintf(
                        "Python module '%s' is not available. %s.",  # Please run wompwomp::setup_python_env(), or %s.
                        package, additional_message
                    ))
                }
            }
        }
    }
}
