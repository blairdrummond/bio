/*
 * Copyright (c) 2019 Corrie daCosta, Blair Drummond
 * 
 * This file is part of scbursts, a library for analysis and * sorting of single-channel bursts.
 * 
 * scbursts is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * scbursts is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with scbursts; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */


#include <iostream>
#include <typeinfo>

#include <exception>

#include "scbursts_types.h"
#include <DCProgsConfig.h>
#include <likelihood/qmatrix.h>
#include <likelihood/missed_eventsG.h>
#include <likelihood/idealG.h>
#include <likelihood/occupancies.h>
#include <likelihood/likelihood.h>

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace DCProgs;

// [[Rcpp::depends(RcppEigen)]]

//' A maximum likelihood estimate for a Q-Matrix with a fixed number
//' of states. Computes the likelihood of a kinetic model given
//' observed bursts, with missed event correction, using the HJCFIT
//' library. This function is an interface for Log10Likelihood from
//' HJCFIT.
//' 
//' https://github.com/DCPROGS/HJCFIT
//' 
//' See the link above and the references therein for details.
//'
//' @param qmatrix Initial guess qmatrix.
//' @param nopen Number of open states in the matrix
//' @param bursts A list of bursts
//' @param tau Maximum length of the missed events
//' @param nmax The exact missed-event likelihood will be computed for 't < n_max * tau'
//' @param xtol Tolerance criteria for brentq().
//' @param rtol Tolerance criteria for brentq().
//' @param itermax Maximum number of iteration when calling brentq().
//' @return A qmatrix kinetic model which numerically maximizes likelihood.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_hjcfit_likelihood_maximize (
		Eigen::MatrixXd qmatrix,
		int nopen,
		DCProgs::t_Bursts bursts,
		DCProgs::t_real tau,
		DCProgs::t_uint nmax=2,
		DCProgs::t_real xtol=0.0000000001,
		DCProgs::t_real rtol=0.0000000001,
		DCProgs::t_uint itermax=100) {
		
		DCProgs::Log10Likelihood likelihood( bursts, nopen, tau, /*tcrit*/ -1, nmax, xtol, rtol, itermax );
		
		Rcout << "Not implemented yet!!!" << std::endl;

		return qmatrix;
		
		//  // When giving names to elements
		//  List L = List::create(
		//  		Named("matrix") = idealG.get_matrix(),
		//  		Named("nopen")  = idealG.get_nopen(),
		//  		Named("nclosed")  = idealG.get_nshut());

}


// [[Rcpp::depends(RcppEigen)]]

//' Computes the likelihood of a kinetic model given observed bursts,
//' with missed event correction, using the HJCFIT library. This
//' function is an interface for Log10Likelihood from HJCFIT.
//' 
//' https://github.com/DCPROGS/HJCFIT
//' 
//' See the link above and the references therein for details.
//'
//' NOTE: It is slow to call this over and over again.
//'
//' @param qmatrix Initial guess qmatrix.
//' @param nopen Number of open states in the matrix
//' @param bursts A list of bursts
//' @param tau Maximum length of the missed events
//' @param nmax The exact missed-event likelihood will be computed for 't < n_max * tau'
//' @param xtol Tolerance criteria for brentq().
//' @param rtol Tolerance criteria for brentq().
//' @param itermax Maximum number of iteration when calling brentq().
//' @return The log10-likelihood of the model given the burst sequence.
//' @export
// [[Rcpp::export]]
DCProgs::t_real cpp_hjcfit_likelihood (
  		Eigen::MatrixXd qmatrix,
		int nopen,
		DCProgs::t_Bursts bursts,
		DCProgs::t_real tau,
		DCProgs::t_uint nmax=2,
		DCProgs::t_real xtol=0.0000000001,
		DCProgs::t_real rtol=0.0000000001,
		DCProgs::t_uint itermax=100) {

		DCProgs::Log10Likelihood likelihood( bursts, nopen, tau, /*tcrit*/ -1, nmax, xtol, rtol, itermax );
		DCProgs::QMatrix _qmatrix(qmatrix, nopen);
		return likelihood(_qmatrix);
		
		//  // When giving names to elements
		//  List L = List::create(
		//  		Named("matrix") = idealG.get_matrix(),
		//  		Named("nopen")  = idealG.get_nopen(),
		//  		Named("nclosed")  = idealG.get_nshut());

}
