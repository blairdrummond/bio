## Copyright (c) 2019 Corrie daCosta, Blair Drummond
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



#' Read a Q-Matrix from a .bin file
#'
#' Reads one of the .bin model files from Mat's simulator
#'
#' @param filename The filename
#' @return The Q-Matrix
#' @examples
#'
#' #' # import some of the data included with the package
#' infile <- system.file("extdata", "example5.bin", package = "scbursts")
#'
#' qm <- qmatrix.read(infile)
#' qm
#' @export
qmatrix.read <- function (filename) {

    ## Should look like this
    ##
	## <%s30uM_model%>
	## {States:0,0,0,0,0,1,1,0}
	## {TRM:-4050,4050,0,0,0,0,0,0}
	## {TRM:1410,-5755,4020,325,0,0,0,0}
	## {TRM:0,25000,-43400,0,18400,0,0,0}
	## {TRM:0,34100,0,-107800,64500,9200,0,0}
	## {TRM:0,0,86400,17200,-228600,0,125000,0}
	## {TRM:0,0,0,27300,0,-27300,0,0}
	## {TRM:0,0,0,0,2100,0,-6600,4500}
	## {TRM:0,0,0,0,0,0,93000,-93000}
	## <%%>

    # load lines
    FileInput <- readLines(filename) 

    title <- sub("<%([^\n]*)%>","\\1", FileInput[[1]])

    tag     <- function (s) { sub("\\{([A-Za-z]*):([^\n]*)\\}.*","\\1", s) }
    content <- function (s) {
        stripped <- sub("\\{([A-Za-z]*):([^\n]*)\\}.*","\\2", s)
        split <- strsplit(stripped, ",")
        # Integer vector
        return(as.integer(unlist(split)))
    }

    mat <- c()
    for (i in 2:length(FileInput)) {

        ### end of file?
        if (FileInput[[i]] == "<%%>") {
            next
        }
        
        t <- tag(FileInput[[i]])
        c <- content(FileInput[[i]])

        if (tolower(t) == "states") {
            states <- as.integer(unlist(c))     # list of states
        } else {
            mat <- c(mat, as.integer(unlist(c))) # append new data
        }
    }

    mat <- qmatrix.create(mat, states)
    attr(mat, "title") <- title
    return(mat)
}


#' Reorder the states of a matrix so that
#' the top left corner contains all the open states.
#'
#' Necessary for HJCFIT.
#'
#' @param mat The matrix.
#' @return The Q-Matrix
#' @examples
#'
#' # import some of the data included with the package
#' infile <- system.file("extdata", "example5.bin", package = "scbursts")
#'
#' qm <- qmatrix.read(infile)
#' qm
#' 
#' qm_r <- qmatrix.reorder(qm)
#' qm_r
#' @export
qmatrix.reorder <- function (mat) {

    states <- qmatrix.states(mat)
    
    for (i in 1:(length(states)-1)) {
        if (states[i] != 1) {
            for (j in (i+1):length(states)) {
                if (states[j] == 1) {
                    warning(sprintf("Switching states %d <-> %d", i, j))
                    
                    states[i] <- 1
                    states[j] <- 0

                    ## row flip
                    temp    <- mat[i,]
                    mat[i,] <- mat[j,]
                    mat[j,] <- temp

                    ## col flip
                    temp    <- mat[,i]
                    mat[,i] <- mat[,j]
                    mat[,j] <- temp

                    # exit the j loop
                    break
                }
            }
        }
    }

    attr(mat, "states") <- states
    return(mat)

}


#' Return the vector of states
#'
#' @param mat The Q-matrix.
#' @return the states vector
#' @examples
#'
#' # import some of the data included with the package
#' infile <- system.file("extdata", "example5.bin", package = "scbursts")
#'
#' qm <- qmatrix.read(infile)
#' qm
#' 
#' qmatrix.states(qm)
#' @export
qmatrix.states <- function (mat) {
    attr(mat, "states")
}



#' Return the number of states (size of the matrix)
#'
#' @param mat The Q-matrix.
#' @return the number of states
#' @examples
#'
#' # import some of the data included with the package
#' infile <- system.file("extdata", "example5.bin", package = "scbursts")
#'
#' qm <- qmatrix.read(infile)
#' qm
#' 
#' qmatrix.states(qm)
#' qmatrix.num_states(qm)
#' @export
qmatrix.num_states <- function (mat) {
    length(attr(mat, "states"))
}


#' Return the number of open states.
#'
#' @param mat The Q-matrix.
#' @return the number of open states
#' @examples
#'
#' # import some of the data included with the package
#' infile <- system.file("extdata", "example5.bin", package = "scbursts")
#'
#' qm <- qmatrix.read(infile)
#' qm
#' 
#' qmatrix.states(qm)
#' qmatrix.num_open(qm)
#' @export
qmatrix.num_open   <- function (mat) {
    length(which(attr(mat, "states") != 0))
}



#' Create a Q-Matrix object from a kinetic rate matrix and vector
#' of associated conductance states.
#'
#' Note that the order of the states will be re-arranged, so that the
#' open states are the in the top-left of the matrix.
#'
#' @param mat The kinetic rate matrix.
#' @param states The associated conductance states.
#' @return the Q-Matrix object (ready for HJCFIT)
#' @examples
#'
#' m = matrix(c(
#'  -3050,  3050,   0,       0,       0,       0,      0,      0,
#'  910,    -4755,  3500,    345,     0,       0,      0,      0,
#'  0,      15000,  -23400,  0,       08400,   0,      0,      0,
#'  0,      24100,  0,       -87800,  54500,   9200,   0,      0,
#'  0,      0,      86400,   17200,   -28600,  0,      25000,  0,
#'  0,      0,      0,       7300,    0,       -7300,  0,      0,
#'  0,      0,      0,       0,       100,     0,      -4600,  4500,
#'  0,      0,      0,       0,       0,       0,      3000,   -3000
#' ), ncol=8, byrow=TRUE)
#' 
#' qm <- qmatrix.create(m, c(0,0,0,0,0,1,1,0))
#' 
#' qm
#' @export
qmatrix.create <- function (mat, states) {

    if (!is.matrix(mat))
        mat <- matrix(mat, ncol=length(states), byrow=TRUE)

    attr(mat, "states") <- states

    qmatrix.sanity_check(mat, no_order=TRUE)

    mat <- qmatrix.reorder(mat)

    return(mat)
}


#' Verify that this is a valid Q-Matrix.
#'
#' @param mat The Q-matrix.
#' @return TRUE if the matrix is ok, FALSE and a warning message otherwise.
#' @examples
#'
#' # import some of the data included with the package
#' infile <- system.file("extdata", "example5.bin", package = "scbursts")
#'
#' qm <- qmatrix.read(infile)
#' qm
#'
#' qmatrix.sanity_check(qm)
#'
#' # entry (3,3) is positive! It should be negative.
#' m = matrix(c(
#'  -3050,  3050,   0,       0,       0,       0,      0,      0,
#'  910,    -4755,  3500,    345,     0,       0,      0,      0,
#'  0,      15000,  23400,   0,       08400,   0,      0,      0,
#'  0,      24100,  0,       -87800,  54500,   9200,   0,      0,
#'  0,      0,      86400,   17200,   -28600,  0,      25000,  0,
#'  0,      0,      0,       7300,    0,       -7300,  0,      0,
#'  0,      0,      0,       0,       100,     0,      -4600,  4500,
#'  0,      0,      0,       0,       0,       0,      3000,   -3000
#' ), ncol=8, byrow=TRUE)
#'
#' # qmatrix.create calls sanity check for you, too.
#' qm <- qmatrix.create(m, c(0,0,0,0,0,1,1,0))
#' qm
#' 
#' qmatrix.sanity_check(qm)
#' @export
qmatrix.sanity_check <- function (mat, no_order=FALSE) {

    states <- qmatrix.states(mat)

    GOOD=TRUE

    ## Sanity check: number of states
    if (length(mat) != length(states)^2) {
        warning("Mismatch of matrix and the number of states.")
        return(FALSE)
    }
    
    ## All rows sum to zero
    ONCE=FALSE
    for(row in 1:nrow(mat)) {
        if (sum(mat[row,]) != 0 && !ONCE) {
            warning("Rows do not sum to zero!")
            ONCE=TRUE
        }
        GOOD=FALSE
    }

    ## Negatives only if on diagonal
    ONCE=FALSE
    for(row in 1:nrow(mat)) {
        for(col in 1:ncol(mat)) {
            if (mat[row,col] < 0 && row != col && !ONCE) {
                warning("Negative entries are not on the diagonal!")
                ONCE=TRUE
                GOOD=FALSE
            }
            if (mat[row,col] > 0 && row == col && !ONCE) {
                warning("Negative entries are not on the diagonal!")
                ONCE=TRUE
                GOOD=FALSE
            }
        }
    }
    
    ## Sanity check: No subconductive states.
    if (any(states != 0 & states != 1)) {
        warning("Subconductive states are not currently supported by hjcfit.")
        warning("We won't stop you, but use of this Q-Matrix comes without warranty.")
        GOOD=FALSE
    }

    if (no_order == FALSE) {
        if (states != sort(states, decreasing=TRUE)) {
            warning("Open states must be in the top left corser.")
            warning("Use qmatrix.reorder to fix this.")
            GOOD=FALSE
        }
    }

    return(GOOD)
}
