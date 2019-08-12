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

#include <likelihood/qmatrix.h>
#include <likelihood/missed_eventsG.h>
#include <likelihood/idealG.h>
#include <likelihood/occupancies.h>

#include <Rcpp.h>
using namespace Rcpp;
 
int main() {

  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;

  DCProgs::QMatrix qmatrix(matrix, /*nopen=*/2);
 
  Rcout << qmatrix.matrix << std::endl;
 
  return 0;
}

//' Testing documentation with Rcpp
//'
//' @return           A string
//' @export
// [[Rcpp::export]]
std::string hjcfit_test () {

  DCProgs::t_rmatrix matrix(5 ,5);
  matrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;

  DCProgs::QMatrix qmatrix(matrix, /*nopen=*/2);
 
  Rcout << qmatrix.matrix << std::endl;

  std::ostringstream ss;
  ss << qmatrix.matrix << std::endl;
  return ss.str();
}
