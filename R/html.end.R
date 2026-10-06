.html.end <- function(file) {
  out <- '</body>\n</html>\n'
  cat(out, file=file, append=TRUE)
}
