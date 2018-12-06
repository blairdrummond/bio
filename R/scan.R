#' Read a scan results text file. scan.read returns a 1 segment dataframe.
#' Reads in scan results and puts them in the same format as the output
#' of dwt.read. See 'dwt', and 'segment' for more information.
#' 
#' Data is in seconds.
#' @param filename, the file name to read from.
#' @param separating_factor In lieu of a known time between segments, 
#'        seperate with a multple of the longest dwell.
#' @return A list of bursts.
#' @examples
#' \dontrun{
#' seg <- scan.read('example_scan_output.scntxt')
#'}
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



    sc             <- segment.create
    init_read      <- read.csv(filename,sep='\t',header=FALSE)
    dwells         <- init_read[[1]]
    max_dwells     <- separating_factor
    max_dwells     <- max(max(dwells)*separating_factor,max_dwells)
    states         <- init_read[[2]]
    for(i in 1:length(states)){if(states[[i]] != 0){states[[i]] <- 1}}
    sd             <- consecutives_to_dwells(states,dwells)
    brst           <- list()
    brst[[1]]      <- sc(sd[[1]],sd[[2]],seg=1,start_time=0,name=util.basename(filename))
    brst           <- bursts.start_times_update(brst,gaps=rep(max_dwells,0))
    return(brst)
}



