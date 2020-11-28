#' @export

runExample <- function() {
  appDir <- system.file("shiny", package = "SpliceSelector")
  if (appDir == "") {
    stop("Could not find shiny directory. Try re-installing `SpliceSelector`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}


