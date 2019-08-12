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



#' Something to do with qmatrix
#'
#' @param qmatrix your qmatrix.
#' @return A fitted model.
#' @importFrom Rcpp sourceCpp
#' @export
hjcfit.fit <- function (qmatrix) {

    ##Sys.setenv("PKG_CXXFLAGS"="-I/usr/local/include/dcprogs/ -I/usr/local/include/eigen3") 
    ##sourceCpp("src/hjcfit_test.cpp")
    ##message(swag_out())
    swag_out()
    return(1)
}
