#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>



using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


struct DC : public Worker {
  
  // input matrix to read from
  const RMatrix<double>  rmat;
  
  // output matrix to write to
  RMatrix<double> D;
  RMatrix<double> C;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  DC(const NumericMatrix rmat, NumericMatrix D, NumericMatrix C)
    : rmat(rmat), D(D), C(C) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      arma::vec d = zeros<vec>(rmat.nrow() );
      for (std::size_t j = 0; j < rmat.nrow(); j++) {
        
        // rows we will operate on
        RMatrix<double>::Row row1 = rmat.row(i);
        RMatrix<double>::Row row2 = rmat.row(j);
        
        // compute the average using std::tranform from the STL
        
        for (size_t k = 0; k < row1.length(); k ++) {
          d[j] += ( (row1[k] - row2[k])*(row1[k] - row2[k]) );
        }
        
      }
      //d = sqrt(d)
      arma::uvec idx = sort_index(d);
      // And return x in that order
      
      for (size_t j = 0; j < C.ncol(); j ++)
      {
        C(i, j) = idx( j + 1 ) + 1;
        D(i, j) = d( idx( j + 1 ) );
      }
    }
  }
};

// [[Rcpp::export]]
List DC_para(NumericMatrix& mat, int k = 10) {
  
  // allocate the matrix we will return
  NumericMatrix D(mat.nrow(), k);
  NumericMatrix C(mat.nrow(), k);
  
  // create the worker
  DC DC(mat, D, C);
  
  // call it with parallelFor
  parallelFor(0, mat.nrow(), DC);
  
  List ret;
  ret["D"] = D;
  ret["C"] = C;
  return ret;
}

// Compute jaccard coefficient between nearest-neighbor sets
//
// Weights of both i->j and j->i are recorded if they have intersection. In this case
// w(i->j) should be equal to w(j->i). In some case i->j has weights while j<-i has no
// intersections, only w(i->j) is recorded. This is determinded in code `if(u>0)`. 
// In this way, the undirected graph is symmetrized by halfing the weight 
// in code `weights(r, 2) = u/(2.0*ncol - u)/2`.
//
// Author: Chen Hao, Date: 25/09/2015


// [[Rcpp::export]]
NumericMatrix jaccard_coeff(NumericMatrix idx) {
  int nrow = idx.nrow(), ncol = idx.ncol();
  NumericMatrix weights(nrow*ncol, 3);
  int r = 0;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int k = idx(i,j)-1;
      NumericVector nodei = idx(i,_);
      NumericVector nodej = idx(k,_);
      int u = intersect(nodei, nodej).size();  // count intersection number
      if(u>0){ 
        weights(r, 0) = i+1;
        weights(r, 1) = k+1;
        weights(r, 2) = u/(2.0*ncol - u)/2;  // symmetrize the graph
        r++;
      }
    }
  }
  
  return weights;
}

// [[Rcpp::export]]
arma::uvec non_zero_index(arma::mat& data) {
  arma::rowvec col_sum = zeros<rowvec>(data.n_cols);
  
  for (size_t i = 0; i < data.n_cols; i++)
  {
    for(size_t j = 0; j < data.n_rows; j++)
    {
      if(data(j,i) != 0)
      {
        col_sum(i) = 1;
        break;
      }
    }
  }
  
  arma::uvec idx = find(col_sum > 0);
  
  return idx;
}

// [[Rcpp::export]]
arma::sp_mat normalize_data_dense(arma::mat& data) {
  arma::rowvec col_sum = zeros<rowvec>(data.n_cols);
  
  for (size_t i = 0; i < data.n_cols; i++)
  {
    for(size_t j = 0; j < data.n_rows; j++)
    {
      if(data(j,i) != 0)
      {
        col_sum(i) = 1;
        break;
      }
    }
  }
  
  arma::uvec idx = find(col_sum > 0);
  data = data.cols(idx);
  
  arma::colvec tmp_max = arma::max(data,1);
  arma::colvec tmp_min = arma::min(data,1);
  arma::colvec tmp_sub = tmp_max-tmp_min;
  
  
  
  data.each_col() -= tmp_min;
  data.each_col() /= tmp_sub;
  
  arma::sp_mat result(data);
  
  return result;
}

// [[Rcpp::export]]
arma::sp_mat normalize_data_sparse(arma::sp_mat& data) {
  arma::rowvec col_sum = zeros<rowvec>(data.n_cols);
  
  for (size_t i = 0; i < data.n_cols; i++)
  {
    for(size_t j = 0; j < data.n_rows; j++)
    {
      if(data(j,i) != 0)
      {
        col_sum(i) = 1;
        break;
      }
    }
  }
  
  arma::uvec idx = find(col_sum == 0);
  for (size_t i = 0; i < idx.n_elem; i++)
  {
    data.shed_col(idx(i) - i);
  }
  
  
  arma::colvec tmp_max = arma::colvec(arma::max(data,1));
  arma::colvec tmp_min = arma::colvec(arma::min(data,1));
  arma::colvec tmp_sub = tmp_max-tmp_min;
  
  for (arma::sp_mat::const_iterator i = data.begin(); i != data.end(); ++i) {
    data(i.row(), i.col()) = (*i - tmp_min(i.row())) / tmp_sub(i.row());
  }
  
  return data;
}