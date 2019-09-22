#include <iostream>
#include <cmath>
#include <range/v3/all.hpp>

template <typename Range>
auto INTG(Range X,Range Q){

  int i, N = X.size();
  double A = 0.0;
  for (size_t i = 1; i < X.size(); ++i){
    A += (Q[i]+Q[i-1])*(X[i]-X[i-1])*0.5;
  }
  return A;
}


