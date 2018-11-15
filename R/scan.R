#' Read a scan results text file. scan.read returns a 1 segment dataframe.
#' Reads in scan results and puts them in the same format as the output
#' of dwt.read. See 'dwt', and 'segment' for more information.
#' 
#' Data is in seconds.
#' @param fn, the file name to read from.
#' @param separating_factor In lieu of a known time between segments, 
#'        seperate with a multple of the longest dwell.
#' @return A list of bursts.
#' @examples
#' \dontrun{
#' seg <- scan.read('filename.txt')
#'}
#' @export
#' @importFrom utils read.csv
scan.read <- function(fn,separating_factor=1000){

    sc             <- segment.create
    init_read      <- read.csv(fn,sep='\t',header=FALSE)
    dwells         <- init_read[[1]]
    max_dwells     <- separating_factor
    max_dwells     <- max(max(dwells)*separating_factor,max_dwells)
    states         <- init_read[[2]]
    for(i in 1:length(states)){if(states[[i]] != 0){states[[i]] <- 1}}
    sd             <- dwt.consecutives_to_dwells(states,dwells)
    brst           <- list()
    brst[[1]]      <- sc(sd[[1]],sd[[2]],seg=1,start_time=0,name=util.basename(fn))
    brst           <- bursts.start_times_update(brst,gaps=rep(max_dwell,0))
    return(brst)
}


#' Intermediate function that runs inside of scan.read. Scan results come in with 
#' possible consecutive closed or open states. Ex: 0,0,0,1,1,1,0,0,1,1,1,...
#' consecutives_to_dwells will sum together dwell times of consecutive closed or open
#' states to produce a list of dwells with alternating states. Ex: 0,1,0,1.
#'
#' @param states,                                     Ex:0,0,0,1,1,1
#' @param dwells (in seconds, same length as states,) Ex:3,5,2,3,4,1
#' @return A list containing two vectors:         States:0,1  ((0,0,0),(1,1,1))
#'                                                Dwells:10,8 ((3+5+2),(3+4+1))
#' @examples
#' \dontrun{
#' alternating_states <- scan.consecutives_to_dwells(states,dwells)
#'}
#' @export
#' @importFrom utils read.csv
scan.consecutives_to_dwells <- function(states,dwells){

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


