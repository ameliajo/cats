//#ifndef CATS_TRAPZ
//#define CATS_TRAPZ

template <typename Range,typename Function>
auto trapz(Range x, Function func){
  double integral = 0.0;
  for (size_t i = 0; i < x.size()-1; ++i){
    integral += 0.5*(func(i)+func(i+1))*(x[i+1]-x[i]);
  }
  return integral;
}

//#endif
