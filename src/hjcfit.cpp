/*
 * Copyright (c) 2019 Corrie daCosta, Blair Drummond
 * 
 * This file is part of scbursts, a library for analysis and
 * sorting of single-channel bursts.
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

#include <nlopt.h>
#include <nlopt.hpp>
// #include <nloptrAPI.h>

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace DCProgs;
using namespace nlopt;

// [[Rcpp::depends(RcppEigen)]]


// used in likelihood wrapper
typedef struct {
	Eigen::MatrixXd qmatrix;
	int nopen;
	int size;
	int nrows;
	double **lookup;
	DCProgs::Log10Likelihood likelihood;
} qmatrix_struct;

// maximize this
double likelihood_wrapper(
	const std::vector<double> &vec,
	std::vector<double> &grad,
	void *data
	) {

	qmatrix_struct *fd = reinterpret_cast<qmatrix_struct*>(data);

	// [&lookup,size,&qmatrix,nrows,nopen,likelihood
	// update the qmatrix
	for (int k=0; k<fd->size; k++) {
		fd->qmatrix(fd->lookup[k][0], fd->lookup[k][1]) = vec.at(k);
	}
	
	// Add the constraint that rows sum to zero
	double row_sum = 0;
	for (int i=0; i<fd->nrows; i++) {
		row_sum = 0;
		for (int j=0; j<fd->nrows; j++) {
			if (i != j) // skip diagonal
				row_sum += fd->qmatrix(i,j);
		}
		fd->qmatrix(i,i) = - row_sum;
	}

	try {
		DCProgs::QMatrix new_qmatrix(fd->qmatrix, fd->nopen);
	} catch (const std::exception& e) { 
		Rcpp::Rcout << "hjcfit: error making matrix." << std::endl; 
		Rcpp::Rcout << "hjcfit: " << e.what() << std::endl; 
		Rcpp::Rcout << fd->qmatrix << std::endl;
		throw;
	}

	try {
		DCProgs::QMatrix new_qmatrix(fd->qmatrix, fd->nopen);
		double result = (double) fd->likelihood(new_qmatrix);
		return result;
	} catch (const std::exception& e) { 
		Rcpp::Rcout << "hjcfit: error computing likelihood." << std::endl; 
		Rcpp::Rcout << "hjcfit: " << e.what() << std::endl; 
		Rcpp::Rcout << fd->qmatrix << std::endl;
		throw;
	}
}


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
//' @param nopen Number of open states in the matrix.
//' @param bursts A list of bursts.
//' @param tau Maximum length of the missed events.
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
	DCProgs::t_uint itermax=100
	) {
		
	DCProgs::Log10Likelihood likelihood( bursts, nopen, tau, /*tcrit*/ -1, nmax, xtol, rtol, itermax );

	int ncols = qmatrix.cols();
	int nrows = qmatrix.rows();
	assert(nrows == ncols);
	
	// How many variable entries?
	int size = 0;
	double row_sum = 0;
	double matrix_max = -1;
	int bad_qmatrix = 0;
	for (int i=0; i<nrows; i++) {
		row_sum = 0;
		for (int j=0; j<ncols; j++) {
			row_sum += qmatrix(i,j);

			if (abs(qmatrix(i,j)) > matrix_max)
				matrix_max = abs(qmatrix(i,j));
			
			if ( qmatrix(i,j) != 0 && i !=j ) {
				if (qmatrix(i,j) < 0) {
					Rcpp::Rcerr << "Off diagonal rates should be positive!" << std::endl;
					bad_qmatrix = 1;
				}
				// Check reversibility
				if ( qmatrix(j,i) == 0 ) {
					Rcpp::Rcerr << "Matrix must be reversible, but " << i << " -> " << j << " is not reversible." << std::endl;
					bad_qmatrix = 1;
				}
				size++;
			}
		}
		if (row_sum != 0) {
			Rcpp::Rcerr << "Rows must sum to zero, but row " << i << " sums to " << row_sum << std::endl;
			bad_qmatrix = 1;
		}
	}
	if (bad_qmatrix) {
		Rcpp::Rcerr << std::endl << "QMatrix needs to be fixed. Returning he original matrix." << std::endl;
		return qmatrix;
	}
	

	/*
	  Now fill all the variable coordinates into a vector,
	  with another vector serving as a lookup table to
	  remember where the coordinates belong in the matrix.
	*/
	double **lookup;
	lookup = new double*[size];

	std::vector<double> maximizer;
	int k = 0;
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<nrows; j++) {
			if (qmatrix(i,j) != 0 && i != j) {
				maximizer.push_back(qmatrix(i,j));
				lookup[k] = new double[2];
				lookup[k][0] = i; 
				lookup[k][1] = j;
				k++;
			}
		}
	}

	qmatrix_struct qmatrix_data = {
		qmatrix, nopen, size, nrows, lookup, likelihood
	};
	
	Rcpp::Rcout << "nlopt: Starting with this matrix:" << std::endl;
	Rcpp::Rcout << (&qmatrix_data)->qmatrix << std::endl;
	Rcpp::Rcout << (&qmatrix_data)->lookup << std::endl;
		
	Rcpp::Rcout << "nlopt: Using COBYLA to maxmize hjcfit (-log_10) likelihood." << std::endl;
	nlopt::opt opt = nlopt::opt(LN_COBYLA, size); 
	
	nlopt::vfunc funky = likelihood_wrapper;

	Rcpp::Rcout << "nlopt: Defining the search space." << std::endl;
	Rcpp::Rcout << "nlopt: 0 <= a_ij < " << 10*matrix_max << std::endl;
	opt.set_lower_bounds(0);
	opt.set_upper_bounds(10*matrix_max);
	opt.set_xtol_rel(1e-4);

	
	Rcpp::Rcout << "nlopt: Setting up objective function." << std::endl;
	opt.set_max_objective(*funky, &qmatrix_data);
	
	double max_likelihood;
	nlopt::result result = opt.optimize(maximizer, max_likelihood);
	Rcpp::Rcout << "nlopt: Optimization finished." << std::endl;

	switch(result) {
	case NLOPT_SUCCESS:
        Rcpp::Rcout << "nlopt result: Success." << std::endl;
		break;

	case NLOPT_STOPVAL_REACHED:
        Rcpp::Rcout << "nlopt result: Optimization stopped because stopval (above) was reached." << std::endl;
		break;

	case NLOPT_FTOL_REACHED:
        Rcpp::Rcout << "nlopt result: Optimization stopped because ftol_rel or ftol_abs (above) was reached." << std::endl;
		break;

	case NLOPT_XTOL_REACHED:
        Rcpp::Rcout << "nlopt result: Optimization stopped because xtol_rel or xtol_abs (above) was reached." << std::endl;
		break;

	case NLOPT_MAXEVAL_REACHED:
        Rcpp::Rcout << "nlopt result: Optimization stopped because maxeval (above) was reached." << std::endl;
		break;

	case NLOPT_MAXTIME_REACHED:
        Rcpp::Rcout << "nlopt result: Optimization stopped because maxtime (above) was reached." << std::endl;
		break;
	}
	
	// Load final result into qmatrix
	for (int k=0; k<size; k++) {
		qmatrix(lookup[k][0], lookup[k][1]) = maximizer.at(k);
	}

	return qmatrix;
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
//' @return The negative log10-likelihood of the model given the burst sequence.
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


