# # Default: Python tests not available
# options(wompwomp.python.ok = FALSE)
# 
# # Check if Python is available
# if (reticulate::py_available(initialize = FALSE)) {
#     # Check if conda is available
#     conda_bin <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
#     if (!is.null(conda_bin) && file.exists(conda_bin)) {
#         # Check if the required conda env exists
#         conda_envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
#         if (!"wompwomp_env" %in% conda_envs) {
#             # Try to create/setup env
#             tryCatch(
#                 {
#                     wompwomp::setup_python_env()
#                     conda_envs <- reticulate::conda_list()$name
#                 },
#                 error = function(e) {
#                     message("Failed to create 'wompwomp_env': ", conditionMessage(e))
#                 }
#             )
#         }
#         
#         # If env exists now, activate it
#         if ("wompwomp_env" %in% conda_envs) {
#             tryCatch({
#                 reticulate::use_condaenv("wompwomp_env", required = TRUE)
#                 Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)
#                 options(wompwomp.python.ok = TRUE)
#             }, error = function(e) {
#                 message("Could not activate 'wompwomp_env': ", conditionMessage(e))
#             })
#         }
#     }
# }
# 
# python_test <- function() {
#     testthat::skip_if_not(
#         getOption("wompwomp.python.ok", FALSE),
#         "Python environment not available"
#     )
# }