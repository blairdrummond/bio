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



#' Computes the likelihood of a kinetic model given observed bursts,
#' with missed event correction, using the HJCFIT library. This
#' function is an interface for Log10Likelihood from HJCFIT.
#' 
#' https://github.com/DCPROGS/HJCFIT
#' 
#' See the link above and the references therein for details.
#'
#' @param qmatrix your qmatrix.
#' @param nopen number of open states.
#' @param bursts a list of bursts.
#' @param resolution Maximum length of the missed events / time resolution.
#' @param nmax The exact missed-event likelihood will be computed for 't < n_max * tau'
#' @param xtol Tolerance criteria for brentq().
#' @param rtol Tolerance criteria for brentq().
#' @param itermax Maximum number of iteration when calling brentq().
#' @return The log10-likelihood of the model given the burst sequence.
#' @export
hjcfit.likelihood <- function (qmatrix, nopen, bursts, resolution, nmax=2, xtol=1e-10, rtol=1e-10, itermax=100) {

    small_intervals <- function (seg) {
        length(subset(seg, dwells < resolution)$dwells) > 0
    }

    if (length(which(unlist(lapply(bursts, small_intervals)) == TRUE)) > 0) {
        warning("Bursts contain intervals small than the time resolution!")
    }
   
    return(cpp_hjcfit_likelihood(qmatrix, nopen, lapply(bursts, function (s) {s$dwells}), resolution, nmax, xtol, rtol, itermax))
    
}





#' A maximum likelihood estimate for a Q-Matrix with a fixed number
#' of states. Computes the likelihood of a kinetic model given
#' observed bursts, with missed event correction, using the HJCFIT
#' library. This function is an interface for Log10Likelihood from
#' HJCFIT.
#' 
#' https://github.com/DCPROGS/HJCFIT
#' 
#' See the link above and the references therein for details.
#'
#' @param qmatrix your initial-guess qmatrix.
#' @param nopen number of open states.
#' @param bursts a list of bursts.
#' @param resolution Maximum length of the missed events / time resolution.
#' @param nmax The exact missed-event likelihood will be computed for 't < n_max * tau'
#' @param xtol Tolerance criteria for brentq().
#' @param rtol Tolerance criteria for brentq().
#' @param itermax Maximum number of iteration when calling brentq().
#' @return A qmatrix kinetic model which numerically maximizes likelihood.
#' @export
hjcfit.fit_model <- function (qmatrix, nopen, bursts, resolution, nmax=2, xtol=10^(-10), rtol=10^(-10), itermax=100) {

    small_intervals <- function (seg) {
        length(subset(seg, dwells < resolution)$dwells) > 0
    }

    if (length(which(unlist(lapply(bursts, small_intervals)) == TRUE)) > 0) {
        warning("Bursts contain intervals small than the time resolution!")
    }
   
    message("Not implemented yet!!!")
    message("TODO: Implement/Find MLE algorithm in C++ and stuff.")
    return(cpp_hjcfit_likelihood_maximize(qmatrix, nopen, lapply(bursts, function (s) {s$dwells}, resolution, nmax, xtol, rtol, itermax)))
    # functional <- hjcfit_likelihood_function(nopen, lapply(bursts, function (s) {s$dwells}), resolution, nmax, xtol, rtol, itermax)
    
}
