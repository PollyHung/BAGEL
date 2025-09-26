#' Setup Logging System
#'
#' @param log_level Character, logging level ("DEBUG", "INFO", "WARN", "ERROR")
#' @param log_file Character, optional file path for logging
#' @export
log_bagel <- function(log_level = "INFO", log_file = NULL) {
  .bagel_env <- new.env(parent = emptyenv())
  .bagel_env$log_level <- log_level
  .bagel_env$log_file <- log_file
}