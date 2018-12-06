#' Remove suffix and path from filename.
#'
#' @param filename string to extract basename from
#' @return Name with suffix and path removed
#' @examples
#' util.basename("bursts/60uM-2017-08-18-16-32/60uM-712.dwt")
#' @export
util.basename <- function(filename) {
    ### Remove the (.dwt|.evt) from filename
    substr(basename(filename), 1, nchar(basename(filename)) - 4) 
}
