/*
 * Copyright (c) 2019 Corrie daCosta, Blair Drummond
 * 
 * This file is part of scbursts, a library for analysis and * sorting of single-channel bursts.
 * This file also makes use of modified version of asa047.cpp, which is licensed under the LGPLv3. See the code and it's own copyright information below.
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

// #include "asa047.hpp"

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

/* FOR asa047.cpp */
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;
/* END of asa047.cpp includes. */





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
//' @param nopen Number of open states in the matrix.
//' @param bursts A list of bursts.
//' @param tau Maximum length of the missed events.
//' @param nmax The exact missed-event likelihood will be computed for 't < n_max * tau'
//' @param xtol Tolerance criteria for brentq().
//' @param rtol Tolerance criteria for brentq().
//' @param itermax Maximum number of iteration when calling brentq().
//' @param REQMIN Nelder-Mead parameter: stop if within required distance.
//' @param KONVGE Nelder-Mead parameter: check convergence every KONVGE steps.
//' @param KMAX Nelder-Mead parameter: maximum number of iterations.
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
	DCProgs::t_uint itermax=100,
	double REQMIN=0.000001,
	int KONVGE=10,
	int KMAX=1000
	) {
		
	DCProgs::Log10Likelihood likelihood( bursts, nopen, tau, /*tcrit*/ -1, nmax, xtol, rtol, itermax );

	int ncols = qmatrix.cols();
	int nrows = qmatrix.rows();
	const int size = ncols * nrows;
		
	double minimizer[size];
	double min_likelihood;

	// All values are <=1, of course.
	double simplex[size];
	for (int fill = fill; fill<size; fill++)
		simplex[fill] = 1.0;

	int how_many_iter = 0;
	int restarts = 0;
	int error_val = -1;


	
		
    // auto fn = [](double *a) { return a[0]*a[0] + a[1]*a[1]; };

	auto likelihood_wrapper = [=](double vec[]) {
								  DCProgs::QMatrix new_qmatrix(
									  Eigen::Map<Eigen::MatrixXd>(vec, nrows, ncols), nopen);
								  return (double) likelihood(new_qmatrix);
							  };
	
	auto     fn      = likelihood_wrapper;    /* function */
	int      n       = size;             /* size */
	double * start   = qmatrix.data();   /* start */
	double * xmin    = minimizer;        /* RESULT matrix */
	double * ynewlo  = &min_likelihood;  /* RESULT likelihood */
	double   reqmin  = REQMIN;           /* Stopping condition */
	double * step    = simplex;          /* parameters <=1 */
	int      konvge  = KONVGE;           /* Check every KONVGE iterations */
	int      kcount  = KMAX;             /* Max number of iterations */                
	int    * icount  = &how_many_iter;   /* Count the number of iterations */
	int    * numres  = &restarts;        /* How many algorithm restarts */
	int    * ifault  = &error_val;       /* 0=OK; 1=Bad Parameter, 2=No Convrgence */

	// // The Nelder-Mead call
	// nelmin(
	// 	X.likelihood_wrapper, /* function */
	// 	size,                  /* size */
	// 	qmatrix.data(),        /* start */
	// 	minimizer,             /* RESULT matrix */
	// 	&min_likelihood,       /* RESULT likelihood */
	// 	REQMIN,                /* Stopping condition */
	// 	simplex,               /* parameters <=1 */
	// 	KONVGE,                /* Check every KONVGE iterations */
	// 	KMAX,                  /* Max number of iterations */                
	// 	&how_many_iter,        /* Count the number of iterations */
	// 	&restarts,             /* How many algorithm restarts */
	// 	&error_val             /* exit code: 0=OK, 1=Illegal Parameter, 2=No Convergence */
	// 	);



	
	/* The following is copied and ever-so-slightly modeifed from asa047.cpp */
	
    // REMOVED BY BLAIR // # include <cstdlib>
    // REMOVED BY BLAIR // # include <iostream>
    // REMOVED BY BLAIR // # include <iomanip>
    // REMOVED BY BLAIR // # include <ctime>
    // REMOVED BY BLAIR // # include <cmath>
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR // using namespace std;
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR // # include "asa047.hpp"
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR // //****************************************************************************80
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR // void nelmin ( auto fn (double * x), int n, double start[], double xmin[], 
    // REMOVED BY BLAIR //   double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
    // REMOVED BY BLAIR //   int *icount, int *numres, int *ifault )
    
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    NELMIN minimizes a function using the Nelder-Mead algorithm.
    //
    //  Discussion:
    //
    //    This routine seeks the minimum value of a user-specified function.
    //
    //    Simplex function minimisation procedure due to Nelder+Mead(1965),
    //    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
    //    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
    //    25, 97) and Hill(1978, 27, 380-2)
    //
    //    The function to be minimized must be defined by a function of
    //    the form
    //
    //      function fn ( x, f )
    //      double fn
    //      double x(*)
    //
    //    and the name of this subroutine must be declared EXTERNAL in the
    //    calling routine and passed as the argument FN.
    //
    //    This routine does not include a termination test using the
    //    fitting of a quadratic surface.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    27 February 2008
    //
    //  Author:
    //
    //    Original FORTRAN77 version by R ONeill.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //
    //    John Nelder, Roger Mead,
    //    A simplex method for function minimization,
    //    Computer Journal,
    //    Volume 7, 1965, pages 308-313.
    //
    //    R ONeill,
    //    Algorithm AS 47:
    //    Function Minimization Using a Simplex Procedure,
    //    Applied Statistics,
    //    Volume 20, Number 3, 1971, pages 338-345.
    //
    //  Parameters:
    //
    //    Input, double FN ( double x[] ), the name of the routine which evaluates
    //    the function to be minimized.
    //
    //    Input, int N, the number of variables.
    //
    //    Input/output, double START[N].  On input, a starting point
    //    for the iteration.  On output, this data may have been overwritten.
    //
    //    Output, double XMIN[N], the coordinates of the point which
    //    is estimated to minimize the function.
    //
    //    Output, double YNEWLO, the minimum value of the function.
    //
    //    Input, double REQMIN, the terminating limit for the variance
    //    of function values.
    //
    //    Input, double STEP[N], determines the size and shape of the
    //    initial simplex.  The relative magnitudes of its elements should reflect
    //    the units of the variables.
    //
    //    Input, int KONVGE, the convergence check is carried out 
    //    every KONVGE iterations.
    //
    //    Input, int KCOUNT, the maximum number of function 
    //    evaluations.
    //
    //    Output, int *ICOUNT, the number of function evaluations 
    //    used.
    //
    //    Output, int *NUMRES, the number of restarts.
    //
    //    Output, int *IFAULT, error indicator.
    //    0, no errors detected.
    //    1, REQMIN, N, or KONVGE has an illegal value.
    //    2, iteration terminated because KCOUNT was exceeded without convergence.
    //
    // REMOVED BY BLAIR // {
      double ccoeff = 0.5;
      double del;
      double dn;
      double dnn;
      double ecoeff = 2.0;
      double eps = 0.001;
      int i;
      int ihi;
      int ilo;
      int j;
      int jcount;
      int l;
      int nn;
      double *p;
      double *p2star;
      double *pbar;
      double *pstar;
      double rcoeff = 1.0;
      double rq;
      double x;
      double *y;
      double y2star;
      double ylo;
      double ystar;
      double z;
    //
    //  Check the input parameters.
    //
      if ( reqmin <= 0.0 )
      {
        *ifault = 1;
		// REMOVED BY BLAIR // return;
		goto end_of_nelder_mead;
      }
    
      if ( n < 1 )
      {
        *ifault = 1;
		// REMOVED BY BLAIR // return;
		goto end_of_nelder_mead;
      }
    
      if ( konvge < 1 )
      {
        *ifault = 1;
		// REMOVED BY BLAIR // return;
		goto end_of_nelder_mead;
      }
    
      p = new double[n*(n+1)];
      pstar = new double[n];
      p2star = new double[n];
      pbar = new double[n];
      y = new double[n+1];
    
      *icount = 0;
      *numres = 0;
    
      jcount = konvge; 
      dn = ( double ) ( n );
      nn = n + 1;
      dnn = ( double ) ( nn );
      del = 1.0;
      rq = reqmin * dn;
    //
    //  Initial or restarted loop.
    //
      for ( ; ; )
      {
        for ( i = 0; i < n; i++ )
        { 
          p[i+n*n] = start[i];
        }
        y[n] = fn ( start );
        *icount = *icount + 1;
    
        for ( j = 0; j < n; j++ )
        {
          x = start[j];
          start[j] = start[j] + step[j] * del;
          for ( i = 0; i < n; i++ )
          {
            p[i+j*n] = start[i];
          }
          y[j] = fn ( start );
          *icount = *icount + 1;
          start[j] = x;
        }
    //                    
    //  The simplex construction is complete.
    //                    
    //  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
    //  the vertex of the simplex to be replaced.
    //                
        ylo = y[0];
        ilo = 0;
    
        for ( i = 1; i < nn; i++ )
        {
          if ( y[i] < ylo )
          {
            ylo = y[i];
            ilo = i;
          }
        }
    //
    //  Inner loop.
    //
        for ( ; ; )
        {
          if ( kcount <= *icount )
          {
            break;
          }
          *ynewlo = y[0];
          ihi = 0;
    
          for ( i = 1; i < nn; i++ )
          {
            if ( *ynewlo < y[i] )
            {
              *ynewlo = y[i];
              ihi = i;
            }
          }
    //
    //  Calculate PBAR, the centroid of the simplex vertices
    //  excepting the vertex with Y value YNEWLO.
    //
          for ( i = 0; i < n; i++ )
          {
            z = 0.0;
            for ( j = 0; j < nn; j++ )
            { 
              z = z + p[i+j*n];
            }
            z = z - p[i+ihi*n];  
            pbar[i] = z / dn;
          }
    //
    //  Reflection through the centroid.
    //
          for ( i = 0; i < n; i++ )
          {
            pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
          }
          ystar = fn ( pstar );
          *icount = *icount + 1;
    //
    //  Successful reflection, so extension.
    //
          if ( ystar < ylo )
          {
            for ( i = 0; i < n; i++ )
            {
              p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
            }
            y2star = fn ( p2star );
            *icount = *icount + 1;
    //
    //  Check extension.
    //
            if ( ystar < y2star )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+ihi*n] = pstar[i];
              }
              y[ihi] = ystar;
            }
    //
    //  Retain extension or contraction.
    //
            else
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+ihi*n] = p2star[i];
              }
              y[ihi] = y2star;
            }
          }
    //
    //  No extension.
    //
          else
          {
            l = 0;
            for ( i = 0; i < nn; i++ )
            {
              if ( ystar < y[i] )
              {
                l = l + 1;
              }
            }
    
            if ( 1 < l )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+ihi*n] = pstar[i];
              }
              y[ihi] = ystar;
            }
    //
    //  Contraction on the Y(IHI) side of the centroid.
    //
            else if ( l == 0 )
            {
              for ( i = 0; i < n; i++ )
              {
                p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
              }
              y2star = fn ( p2star );
              *icount = *icount + 1;
    //
    //  Contract the whole simplex.
    //
              if ( y[ihi] < y2star )
              {
                for ( j = 0; j < nn; j++ )
                {
                  for ( i = 0; i < n; i++ )
                  {
                    p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                    xmin[i] = p[i+j*n];
                  }
                  y[j] = fn ( xmin );
                  *icount = *icount + 1;
                }
                ylo = y[0];
                ilo = 0;
    
                for ( i = 1; i < nn; i++ )
                {
                  if ( y[i] < ylo )
                  {
                    ylo = y[i];
                    ilo = i;
                  }
                }
                continue;
              }
    //
    //  Retain contraction.
    //
              else
              {
                for ( i = 0; i < n; i++ )
                {
                  p[i+ihi*n] = p2star[i];
                }
                y[ihi] = y2star;
              }
            }
    //
    //  Contraction on the reflection side of the centroid.
    //
            else if ( l == 1 )
            {
              for ( i = 0; i < n; i++ )
              {
                p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
              }
              y2star = fn ( p2star );
              *icount = *icount + 1;
    //
    //  Retain reflection?
    //
              if ( y2star <= ystar )
              {
                for ( i = 0; i < n; i++ )
                {
                  p[i+ihi*n] = p2star[i];
                }
                y[ihi] = y2star;
              }
              else
              {
                for ( i = 0; i < n; i++ )
                {
                  p[i+ihi*n] = pstar[i];
                }
                y[ihi] = ystar;
              }
            }
          }
    //
    //  Check if YLO improved.
    //
          if ( y[ihi] < ylo )
          {
            ylo = y[ihi];
            ilo = ihi;
          }
          jcount = jcount - 1;
    
          if ( 0 < jcount )
          {
            continue;
          }
    //
    //  Check to see if minimum reached.
    //
          if ( *icount <= kcount )
          {
            jcount = konvge;
    
            z = 0.0;
            for ( i = 0; i < nn; i++ )
            {
              z = z + y[i];
            }
            x = z / dnn;
    
            z = 0.0;
            for ( i = 0; i < nn; i++ )
            {
              z = z + pow ( y[i] - x, 2 );
            }
    
            if ( z <= rq )
            {
              break;
            }
          }
        }
    //
    //  Factorial tests to check that YNEWLO is a local minimum.
    //
        for ( i = 0; i < n; i++ )
        {
          xmin[i] = p[i+ilo*n];
        }
        *ynewlo = y[ilo];
    
        if ( kcount < *icount )
        {
          *ifault = 2;
          break;
        }
    
        *ifault = 0;
    
        for ( i = 0; i < n; i++ )
        {
          del = step[i] * eps;
          xmin[i] = xmin[i] + del;
          z = fn ( xmin );
          *icount = *icount + 1;
          if ( z < *ynewlo )
          {
            *ifault = 2;
            break;
          }
          xmin[i] = xmin[i] - del - del;
          z = fn ( xmin );
          *icount = *icount + 1;
          if ( z < *ynewlo )
          {
            *ifault = 2;
            break;
          }
          xmin[i] = xmin[i] + del;
        }
    
        if ( *ifault == 0 )
        {
          break;
        }
    //
    //  Restart the procedure.
    //
        for ( i = 0; i < n; i++ )
        {
          start[i] = xmin[i];
        }
        del = eps;
        *numres = *numres + 1;
      }
      delete [] p;
      delete [] pstar;
      delete [] p2star;
      delete [] pbar;
      delete [] y;
    
    // REMOVED BY BLAIR //   return;
    // REMOVED BY BLAIR // }
    // REMOVED BY BLAIR // //****************************************************************************80
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR // void timestamp ( void )
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR // //****************************************************************************80
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //  Purpose:
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //    TIMESTAMP prints the current YMDHMS date as a time stamp.
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //  Example:
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //    31 May 2001 09:45:54 AM
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //  Licensing:
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //    This code is distributed under the GNU LGPL license. 
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //  Modified:
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //    24 September 2003
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //  Author:
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //    John Burkardt
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //  Parameters:
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // //    None
    // REMOVED BY BLAIR // //
    // REMOVED BY BLAIR // {
    // REMOVED BY BLAIR // # define TIME_SIZE 40
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR //   static char time_buffer[TIME_SIZE];
    // REMOVED BY BLAIR //   const struct tm *tm;
    // REMOVED BY BLAIR //   size_t len;
    // REMOVED BY BLAIR //   time_t now;
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR //   now = time ( NULL );
    // REMOVED BY BLAIR //   tm = localtime ( &now );
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR //   len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR //   cout << time_buffer << "\n";
    // REMOVED BY BLAIR // 
    // REMOVED BY BLAIR //   return;
    // REMOVED BY BLAIR // # undef TIME_SIZE
    // REMOVED BY BLAIR // }
    /* END of asa047.cpp */
end_of_nelder_mead:
	
	if (error_val == 1) 
		Rcout << "Illegal input values!!!" << std::endl;
	else if (error_val == 2)
		Rcout << "Failed to Converge!!!" << std::endl;

	return Eigen::Map<Eigen::MatrixXd>(minimizer, nrows, ncols);

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
