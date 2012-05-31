#' @examples \dontrun{
#' html.img("index.html", class="myid", src="my.png", width=400)
#' }
html.img <- function(file, ...) {
  att <- unlist(list(...))
  out <- '<img'
  out <- c(out, paste(names(att), dQuote(att), sep="="))
  out <- c(out, '/>')
  out <- paste(out, collapse=" ")
  cat(out, "\n", file=file, append=TRUE)
}
