#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double auc_cpp(NumericVector times, NumericVector outputs, NumericVector range = NumericVector::create()) {
  double start_time, end_time;
  
  if(times.size() != outputs.size()) stop("times and output must have same length");
  
  if(range.size() == 0)
  {
    start_time = min(times);
    end_time = max(times);
  }
  else
  {
    start_time = range[0];
    end_time = range[1];
  }
  
  int n = times.size();
  double auc = 0;
  
  for(int i = 0; i < n-1; ++i)
  {
    if(times[i+1] > end_time) return(auc);
    if(times[i] < start_time) continue;
    if(R_IsNA(outputs[i])) return(NA_REAL);
    
    auc += (outputs[i] + outputs[i+1])*(times[i+1] - times[i])/2;
  }
  
  return auc;
}