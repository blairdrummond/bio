#include <iostream>
#include <Eigen/Dense>


#include <Rcpp.h>
using namespace Rcpp;
using Eigen::MatrixXd;

//' Testing documentation with Rcpp
//'
//' @return           A number
//' @export
// [[Rcpp::export]]
int eigen_test()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  
  Rcout << m << std::endl;
  return 0;

}
