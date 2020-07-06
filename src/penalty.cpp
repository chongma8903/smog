//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins("cpp11")]]

#include <RcppArmadillo.h>
#include "penalty.h"

using namespace Rcpp;
using namespace arma;

//' proximal operator on L1 penalty
//' 
//' @param x numeric value.
//' @param lambda numeric value for the L1 penalty parameter.
//' @keywords internal
//' 
//[[Rcpp::export]]
double proxL1(const double &x, const double &lambda){
  double res = std::fabs(x) > lambda ? ( x > lambda ? x-lambda : x+lambda ) : 0;
  return res;
}

//' proximal operator on L2 penalty
//' 
//' @param x A numeric vector.
//' @param lambda numeric value for the L2 penalty parameter.
//' @keywords internal
//' 
//[[Rcpp::export]]
arma::vec proxL2(const arma::vec &x, const double &lambda){
  double thr = 1 - lambda*std::sqrt(x.n_elem)/arma::norm(x,2);
  if(thr > 0){
    return thr*x;
  }else{
    return arma::zeros(x.n_elem);
  }
}

//' proximal operator on the composite L2, L2-Square, and L1 penalties
//' 
//' @param x A numeric vector of two.
//' @param lambda a vector of three penalty parameters. \eqn{\lambda_1} and 
//'        \eqn{\lambda_2} are L2 and L2-Square (ridge) penalties for \eqn{x} in 
//'        a group level, and \eqn{\lambda_3} is the L1 penalty for \eqn{x_2}, respectively.
//' @param hierarchy a factor value in levels 0, 1, 2, which represent different
//'                  hierarchical structure in x, respectively. When \code{hierarchy=0},
//'                  \eqn{\lambda_2} and \eqn{\lambda_3} are forced to be zeroes; when
//'                  \code{hierarchy=1}, \eqn{\lambda_2} is forced to be zero; when 
//'                  \code{hierarchy=2}, there is no constraint on \eqn{\lambda}'s. 
//'                  See \code{\link{smog}}.
//' @param d indices for overlapped variables in x.   
//' 
//' @seealso \code{\link{cv.smog}}, \code{\link{smog.default}}, \code{\link{smog.formula}}, 
//'          \code{\link{predict.smog}}, \code{\link{plot.smog}}.
//' 
//' @author Chong Ma, \email{chong.ma@@yale.edu}.
//' @references \insertRef{ma2019structural}{smog}
//' 
//' 
//[[Rcpp::export]]
arma::vec prox(const arma::vec &x, const arma::vec &lambda, 
               const int &hierarchy, const arma::uvec &d){
  arma::vec res = x;
  switch(hierarchy){
  case 0:{
    res = proxL2(res,lambda[0]);
    break; 
  }
  case 1:{
    res.elem(d) = proxL2(res.elem(d),lambda[2]);
    res = proxL2(res,lambda[0]);
    break;
  }
  case 2:{
    res.elem(d) = proxL2(res.elem(d),lambda[2]);
    double thr = 1 - lambda[0]*std::sqrt(res.n_elem)/arma::norm(res,2);
    if(thr > 0){
      res = 1.0/(1+2*lambda[1])*thr*res;
    }else{
      res = arma::zeros(res.n_elem);
    }
    break;
  }
  default:{
    Rcpp::stop("hierarchy must be a value from 0, 1, 2");
  }
  }
  
  return res;
}


//' Penalty function on the composite L2, L2-Square, and L1 penalties
//' 
//' @param x A vector of two numeric values, in which \eqn{x_1} represents
//'          the prognostic effect, and \eqn{x_2} for the predictive effect, 
//'          respectively. 
//' @param lambda a vector of three penalty parameters. \eqn{\lambda_1} and 
//'        \eqn{\lambda_2} are L2 and L2-Square (ridge) penalties for \eqn{x} in 
//'        a group level, and \eqn{\lambda_3} is the L1 penalty for \eqn{x_2}, respectively.
//' @param hierarchy a factor value in levels 0, 1, 2, which represent different
//'                  hierarchical structure in x, respectively. When \code{hierarchy=0},
//'                  \eqn{\lambda_2} and \eqn{\lambda_3} are forced to be zeroes; when
//'                  \code{hierarchy=1}, \eqn{\lambda_2} is forced to be zero; when 
//'                  \code{hierarchy=2}, there is no constraint on \eqn{\lambda}'s. 
//'                  See \code{\link{smog}}.
//' @param d indices for overlapped variables in x. 
//' 
//' @seealso \code{\link{cv.smog}}, \code{\link{smog.default}}, \code{\link{smog.formula}}, 
//'          \code{\link{predict.smog}}, \code{\link{plot.smog}}.
//' 
//' @author Chong Ma, \email{chong.ma@@yale.edu}.
//' @references \insertRef{ma2019structural}{smog}
//' 
//[[Rcpp::export]]
double penalty(const arma::vec &x, const arma::vec &lambda, 
               const int &hierarchy, const arma::uvec &d){
  double res;
  switch(hierarchy){
  case 0:{
    res = arma::norm(x,2)*lambda[0];
    break;
  }
  case 1:{
    res = arma::norm(x,2)*lambda[0] + arma::norm(x.elem(d),2)*lambda[2];
    break;
  }
  case 2:{
    res = arma::norm(x,2)*lambda[0] + std::pow(arma::norm(x,2),2)*lambda[1] + arma::norm(x.elem(d),2)*lambda[2];
    break;
  }
  default:{
    Rcpp::stop("hierarchy must be a value from 0, 1, 2");
  }
  }
  
  return res;
}

