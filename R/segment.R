#' Create a "segment" object
#' 
#' This is a low-level function, mostly for use internally by other functions. There aren't many reasons to use this. Create object containing table data and metadata. The object can be used as a dataframe, and the metadata can be accessed with the functions: segment.seg, segment.start_time, segment.filename
#' 
#' @param states a vector of states in {0,1}
#' @param dwells a vector of dwell durations (same length as states)
#' @param seg The segment number. Defaults to 1
#' @param start_time When the dwells began. Defaults to 0
#' @param name Suffix-less version of the original filename. 60uM.dwt -> '60uM'
#' @param ignore_errors Do not report faulty segments (not many reasons to do this)
#' @return The segment object: A dataframe with extra metadata.
#' @examples
#'
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=0, name="example_segment")
#' 
#' segment.name(my_burst)
#' 
#' @export
segment.create <- function (states, dwells, seg=1, start_time=0, name="burst", ignore_errors=FALSE) {

    data  <- data.frame(states, dwells)

    attr(data, "name") <- name
    attr(data, "seg")  <- seg
    attr(data, "start_time")  <- start_time

    # This is reasonable to include, but will interfere with the check done in bursts.
    if (!ignore_errors && !segment.verify(data))
        warning('Burst seems to have been misrecorded!')
    
    return(data)
    
}




#' Extract segment number from segment.
#'
#' @param segment the segment object
#' @return Segment number (integer)
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=0, name="example_segment")
#' 
#' segment.seg(my_burst)
#' 
#' @export
segment.seg <- function(segment) {attr(segment, "seg")}




#' Extract start_time from segment.
#'
#' @param segment the segment object
#' @return Segment start_time (float)
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.start_time(my_burst)
#' 
#' @export
segment.start_time <- function(segment) {attr(segment, "start_time")}



#' Extract name from segment.
#'
#' @param segment the segment object
#' @return Segment name (string)
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.name(my_burst)
#' 
#' @export
segment.name <- function(segment) {attr(segment, "name")}



#' Get duration of a segment.
#'
#' @param segment the segment object
#' @return the duration
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.duration(my_burst)
#' 
#' @export
segment.duration <- function(segment) {
    sum(segment$dwells)
}



#' Extract number of dwells in segment.
#'
#' @param segment the segment object
#' @return number of dwells
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.count_dwells(my_burst)
#' 
#' @export
segment.count_dwells <- function(segment) {length(segment$states)}



#' Extract number of open dwells.
#'
#' @param segment the segment object
#' @return number of open dwells
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.count_open(my_burst)
#' 
#' @export
segment.count_open<- function(segment) {sum(segment$states == 1)}




#' Extract number of closed dwells.
#'
#' @param segment the segment object
#' @return number of closed dwells
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.count_closed(my_burst)
#' 
#' @export
segment.count_closed<- function(segment) {sum(segment$states == 0)}





#' Extract open dwells.
#'
#' @param segment the segment object
#' @return the open dwells
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' open_dwells <- segment.open_dwells(my_burst)
#' head(open_dwells)
#' 
#' @export
segment.open_dwells <- function(segment) { subset(segment, states == 1)$dwells }



#' Extract closed dwells.
#'
#' @param segment the segment object
#' @return the closed dwells
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' closed_dwells <- segment.closed_dwells(my_burst)
#' head(closed_dwells)
#' 
#' @export
segment.closed_dwells <- function(segment) { subset(segment, states == 0)$dwells }




#' Calculate empirical P(Open) of a segment.
#' 
#' Calculate empirical P(Open) of a segment. NOTE: Assuming that burst starts and ends with 1
#'
#' @param segment The dwells and states table
#' @return The ratio of open time to total time
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#'
#' # P(Open) of this burst
#' segment.popen(my_burst)
#' 
#' @export
segment.popen <- function (segment) {

    open_times <- subset(segment, states == 1, select=dwells)

    total_duration <- sum(segment$dwells)
    
    return (sum(open_times) / total_duration)
}




#' Calculate empirical P(Closed) of a segment.
#' 
#' Calculate empirical P(Closed) of a segment. NOTE: Assuming that burst starts and ends with 1
#'
#' @param segment The dwells and states table
#' @return The ratio of closed time to total time
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' # P(Closed) of this burst
#' segment.pclosed(my_burst)
#' 
#' @export
segment.pclosed <- function (segment) {

    popen <- segment.popen(segment)
    
    return ( 1 - popen )
}


#' Detect misrecorded data.
#' 
#' Segments should have a very specific shape, but recordings can produce errors that make non-sensical segments. In particular, ones contain multiple consecutive openings or closings, or end in closings. This function detects whether a segment satisfies the constraint that the segment states alternate between open and closed, and begin and end with a closing.
#'
#' @param segment The dwells and states table
#' @return True if a valid segment, False otherwise
#' @examples
#' 
#' # It's more likely that you created states or dwells with some function
#' states  <-  c(0,      1,    0,    1,    0,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.0,  1.1,  0.6,  1.1,  0.8,  1.1)
#' my_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="example_segment")
#' 
#' segment.verify(my_burst)
#'
#' # Now, a bad burst with two adjacent open dwells
#' states  <-  c(0,      1,    0,    1,    1,    0,    1,    0,    1)
#' dwells  <-  c(0.1,  1.1,  0.5,  0.2,  1.1,  0.6,  1.1,  0.8,  1.1)
#'
#' # This will issue a warning
#' faulty_burst <- segment.create(states, dwells, seg=1, start_time=3.14159, name="faulty_segment")
#'
#' # This will differentiate good and faulty bursts
#' segment.verify(faulty_burst)
#'
#' # If you have a list of bursts, you can select the good ones with
#' # vbursts <- bursts.select(bursts, segment.verify)
#' 
#' @export
#' @importFrom utils head tail
segment.verify <- function (segment) {

    if (length(segment$states) == 0)
        return(TRUE)
    
    ## Begins with 0
    if (head(segment$states,n=1) == 0)
        return(FALSE)

    ## Ends with 0
    if (tail(segment$states,n=1) == 0)
        return(FALSE)

    ## Contains consecutive 0s or 1s.
    if (any(diff(segment$states) == 0))
        return(FALSE)

    ## Otherwise is OK
    return(TRUE)
}


#' Collapses a segment into dwells with alternating conductance level.
#'
#' Segments may contain consecutive dwells with the same conductance level. 
#' Consecutives_to_dwells sums together all consecutive dwells with the same
#' conductance level. The result is a segment containing dwells that alternate
#' in conductance level (i.e. 1,0,1,0,1,...)
#'
#' @param  segment The dwells and states table
#' @return A modified copy of the original segment 
#' @export
segment.consecutives_to_dwells <- function(segment){

    s2 = copy(segment)
    d  = s2$dwells
    s  = s2$states
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

    s2$dwells = sd[[2]]
    s2$states = sd[[1]]

    return(s2)
}


#' Imposes a deadtime to a segment by removing any dwell that is shorter than the deadtime.
#'
#' The user specifies a deadtime in microseconds. The function effectively undoes
#' the work of the event detection algorithm by reverting the conductance level 
#' (of the brief dwell) back to the previous conductance level in the time sequence.
#' The function then returns a collapsed segment containing alternating dwells.
#' (See segment.consecutives_to_dwells for details about the collapsed segment.)
#' 
#'
#' @param  segment, the segment containing dwells and states.
#' @param  dead_time, the briefest possible event in microseconds.
#' @return A modified copy of the original segment 
#' @export
segment.impose_deadtime <- function(segment,dead_time){

    s2 = copy(segment)
    d  = s2$dwells
    s  = s2$states
    dt = dead_time / 1e6

    #If first dwell is < dead_time, then set the conductance level to 0.
    if (d[[1]] < dt){s[[1]] = 0}
       
    #For all other dwells in the segment, if the dwell is < dead_time, change the 
    #conductance level to the previous (in the time sequence) dwell's 
    #conductance level. 
    for (i in 2:length(d)){if (d[[i]] < dt){s[[i]] = s[[i-1]]}}

    #returns a collapsed segment of dwells with alternating conductance levels.
    return(segment.consecutives_to_dwells(s2))
}

#' Imposes a fixed conductance level (0 or 1) to all dwells with subconductance levels.
#'
#' The user specifies the desired level ('open' or 'closed'). The function will modify
#' any subconductance level (that is not 0 or 1) to be the desired level 1 for 'open'
#' or 0 for 'closed'. The function then reutrns a collapsed segment containing 
#' alternating dwells.
#' (See segment.consecutives_to_dwells for details about the collapsed segment.)
#'
#' @param  segment, the segment containing dwells and states.
#' @param  level, either 'open' or 'closed'
#' @return A modified copy of the original segment 
#' @export
segment.subconductance_as <- function(segment,level){

    s2 = copy(segment)
    d  = s2$dwells
    s  = s2$states


    #Sets desired conductance level to 1 or zero depending on users' choice.
    #Returns a warning message if level is not 'open' or 'closed'
    if      (level == 'open')   {l = 1}
    else if (level == 'closed') {l = 0}
    else    {return('Conductance level must be either \'open\' or \'closed\'.')}

    #For all the dwells in the segment, if the conductance level is not a 1 or a 0,
    #then set the conductance level as the desired level l.
    for (i in 1:length(s)){if ((s[[i]] != 0) & (s[[i]] != 1)){s[[i]] = l}}

    #Returns a collapsed segment of dwells with alternating conductance levels.
    return(segment.consecutives_to_dwells(s2))
}
