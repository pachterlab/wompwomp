required_packages <- c(
    # Imports
    "dplyr",
    "ggplot2",
    "ggforce",
    "ggalluvial",
    "igraph",
    "tibble",
    "rlang",
    "tidyr",
    "magrittr",
    "data.table",
    "purrr",
    "ggfittext",
    "TSP",

    # Suggests
    "testthat",
    "vdiffr",
    "magick",
    "knitr",
    "rmarkdown",
    "gtools",
    "mclust",
    "here",
    "ggrastr",
    "stringr",
    "sessioninfo",
    "remotes"
)

# Optionally filter out already-installed packages
to_install <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
}

# Optional: check version requirement for testthat (>= 3.0.0)
if ("testthat" %in% installed.packages()[,"Package"]) {
    if (packageVersion("testthat") < "3.0.0") {
        install.packages("testthat", repos = "https://cloud.r-project.org")
    }
}

# Optional R version check
if (getRversion() < "3.5.0") {
    stop("R >= 3.5.0 is required.")
}
