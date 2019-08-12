## Copyright (c) 2019 Corrie daCosta, Mathieu Dextraze, Blair Drummond
##
## This file is part of scbursts, a library for analysis and
## sorting of single-channel bursts.
##
## scbursts is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
##
## scbursts is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public
## License along with scbursts; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA


#' Read a .xlsx file output from clampfit
#' 
#' Read a .xlsx file output from clampfit. Result is a list of "segments", which is a dataframe extra data. See "segment" for more details. Converts millisecond dwells to seconds.
#'
#' @param filename Filename to read from
#' @param separating_factor In lieu of a known time between segments, seperate with a multple of the longest dwell.
#' @param header Does the file include a header?
#' @return A list of bursts (possibly a singleton)
#' @examples
#' 
#' infile <- system.file("extdata", "example1_clampfit.xlsx", package = "scbursts")
#' dwells <- clampfit.read(infile)
#' head(dwells)
#' 
#' @export
#' @importFrom readxl read_excel 
#' @importFrom tibble as_tibble
clampfit.read <- function(filename, separating_factor=1000, header=FALSE) {

    i_read           <- read_excel(filename,sheet=1,col_names=header) #read in the .xlsx file
    names(i_read)[3] <- 'states'
    names(i_read)[9] <- 'dwells'
    dwells           <- tibble::as_tibble(i_read[9][,1])$dwells # column 9 are the dwells
    dwells           <- dwells / 1000 # milliseconds to seconds
    max_dwells       <- separating_factor
    max_dwells       <- max(max(dwells)*separating_factor, max_dwells)
    states           <- tibble::as_tibble(i_read[3][,1])$states # column 3 are the conductance levels

    ## Remove leading 0
    if (states[1] == 0) {
        warning("Removing leading 0 dwell.")
        states  <- states[2:length(states)]
        dwells  <- dwells[2:length(dwells)]
    }
    ## Remove trailing 0
    if (states[length(states)] == 0) {
        warning("Removing trailing 0 dwell.")
        states  <- states[1:length(states)-1]
        dwells  <- dwells[1:length(dwells)-1]
    }
        
    brst             <- list()
    brst[[1]]        <- segment.create(states,dwells,seg=1,start_time=0,name=util.basename(filename))
    brst             <- bursts.start_times_update(brst,gaps=rep(max_dwells,0))
    return(brst)

}
