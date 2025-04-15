#' Create Background Files from Breakpoint Data
#'
#' This function processes breakpoint files generated during preprocessing,
#' categorizing them into telomere and centromere segments. It also creates
#' a directory to store the output files if it does not already exist.
#'
#' @param result_dir A string specifying the directory where the breakpoint
#' files are located and where the output background files will be saved.
#' @return None; the function saves two files: one for telomere and one for
#' centromere backgrounds in the specified result directory.
#' @export
background <- function(result_dir) {

  message("Creating Lineage Specific Backgrounds")

  ## Create the result directory if it does not exist
  if (!dir.exists(file.path(result_dir, "backgrounds"))) {
    dir.create(file.path(result_dir, "backgrounds"))
  }

  ## Read in all the preprocess produced Breakpoint files
  files <- list.files(file.path(result_dir, "breakpoints"), full.names = TRUE)
  breakpoints <- lapply(files, function(x) {
    df <- read.delim(x)
    df$amp_del <- ifelse(grepl("_amp", x), "amp", "del")
    df$tel_cent <- ifelse(grepl("_tel", x), "tel", "cent")
    return(df)
  })
  breakpoints <- do.call(rbind, breakpoints)

  ## Split the breakpoint files into telomere and centromere
  tel <- breakpoints[breakpoints$tel_cent == 'tel', c(1:5)]
  cent <- breakpoints[breakpoints$tel_cent == 'cent', c(1:5)]

  ## Save files
  message("writing telomere-bound backgrounds")
  write.table(tel, file.path(result_dir, "backgrounds/background_telomere.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  message("writing centromere-bound backgrounds")
  write.table(cent, file.path(result_dir, "backgrounds/background_centromere.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

}
