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



#' Read a scan results text file. scan.read returns a 1 segment list
#' Reads in scan results and puts them in the same format as the output
#' of dwt.read. See 'dwt', and 'segment' for more information.
#' 
#' Data is in seconds.
#' @param filename, the file name to read from.
#' @param separating_factor In lieu of a known time between segments, 
#'        seperate with a multple of the longest dwell.
#' @return A list of recording segments from the scan file
#' @examples
#' 
#' infile <- system.file("extdata", "example1_scan.txt", package = "scbursts")
#' record <- scan.read(infile)
#' head(record)
#'
#' @export
#' @importFrom utils read.csv
scan.read <- function(filename,separating_factor=1000){

    consecutives_to_dwells <- function(states,dwells){

        d  = dwells
        s  = states
        sd = list(vector(),vector())
        l  = length
        c  = 1
        i  = 1
        while(i < l(s)){
            if(s[i] == s[i+1]){
                d[i+1] = d[i+1]+d[i]
                if((i+1) == l(s)){
                    sd[[1]] = append(sd[[1]],s[i+1])
                    sd[[2]] = append(sd[[2]],d[i+1])}
                i = i+1}
            else if(s[i]!= s[i+1]){
                    sd[[1]] = append(sd[[1]],s[i])
                    sd[[2]] = append(sd[[2]],d[i])
                    if((i+1) == l(s)){
                        sd[[1]] = append(sd[[1]],s[i+1])
                        sd[[2]] = append(sd[[2]],d[i+1])}
                    i = i+1}
        }
        return(sd)
    }
    
    is.binary <- function(fn,max=1000){
        
        f = file(fn,'rb',raw=TRUE)
        b = readBin(f,'int',max,size=1,signed=FALSE)
        close(f)
        return(max(b)>128)
    }

    ### Warning Message: If file is not in text format:
    if (is.binary(filename)){
        stop('Input file must be a text file.')
    }

    init_read      <- read.csv(filename,sep='\t',header=FALSE)
    dwells         <- init_read[[1]]
    max_dwells     <- separating_factor
    max_dwells     <- max(max(dwells)*separating_factor,max_dwells)
    states         <- init_read[[2]]
    for(i in 1:length(states)){if(states[[i]] != 0){states[[i]] <- 1}}
    sd             <- consecutives_to_dwells(states,dwells)
    brst           <- list()
    states         <- sd[[1]]
    dwells         <- sd[[2]]
    
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
    
    brst[[1]]      <- segment.create(states,dwells,seg=1,start_time=0,name=util.basename(filename))
    brst           <- bursts.start_times_update(brst,gaps=rep(max_dwells,0))
    return(brst)
}



