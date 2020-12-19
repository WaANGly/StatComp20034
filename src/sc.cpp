#include <Rcpp.h>
using namespace Rcpp;

//' @title  A Metropolis sampler using Rcpp
//' @description  A Metropolis sampler using Rcpp
//' @param N the number of samples
//' @param x the number of between-sample random number
//' @param sigma variance
//' @return a numericmatrix with a random sample of size \code{n} and k
//' @examples
//' \dontrun{
//' sigma <- 2
//' x <- 25
//' N <- 2000
//' res <- rwc(sigma, x, N)
//' }
//' @useDynLib StatComp20034, .registration = TRUE
//' @exportPattern "^[[:alpha:]]+"
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
NumericMatrix rwc(double sigma,int x,int N ) {
  NumericMatrix mat(N,2);
  mat(0,0) = x;
  double k =0;
  for (int i=1;i<N;i++){
    double u = runif(1)[0];
    double y = rnorm(1,mat(i-1,0),sigma)[0];
    double r1=fabs(mat(i-1,0));
    double r2=abs(y);
    double r=r1-r2;
    if (u <= exp(r)){
      mat(i,0) = y;
    }
    else {
      mat(i,0) = mat(i-1,0);
      k = k+1;
    }
    mat(i,1) = k;
  } 
  return mat;
}


