#' Read a .xlsx file output from clampfit
#' 
#' Read a .xlsx file output from clampfit. Result is a list of "segments", which is a dataframe extra data. See "segment" for more details. Converts millisecond dwells to seconds.
#'
#' @param filename Filename to read from
#' @param separating_factor In lieu of a known time between segments, seperate with a multple of the longest dwell.
#' @return A list of bursts (possibly a singleton)
#' @examples
#' \dontrun{
#' dwells <- clampfit.read('example.xlsx')
#' }
#' @export
#' @importFrom gdata read.xls
clampfit.read <- function(filename, separating_factor=1000) {

    sc               <- segment.create
    i_read           <- read.xls(filename,sheet=1,header=TRUE) #read in the .xlsx file
    names(i_read)[3] <- 'states'
    names(i_read)[9] <- 'dwells'
    dwells           <- i_read[9] # column 9 are the dwells
    dwells           <- dwells / 1000 # milliseconds to seconds
    max_dwells       <- separating_factor
    max_dwells       <- max(max(dwells)*separating_factor,max_dwells)
    states           <- i_read[3] # column 3 are the conductance levels
    brst             <- list()
    brst[[1]]      <- sc(states,dwells,seg=1,start_time=0,name=util.basename(filename))
    brst             <- bursts.start_times_update(brst,gaps=rep(max_dwell,0))
    return(brst)

}
