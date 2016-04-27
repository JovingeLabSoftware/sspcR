#' Trim leading and trailing whitespace
#'
#' Similar to \emph{trim} functions in other languages, \code{trim} removes
#' leading and trailing whitespace from character vectors.
#' @export
#' @examples
#' spacey_string = " hi there  "
#' trim(spacey_string)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
