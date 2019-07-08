#ifndef CATS_INTERPOLATE
#define CATS_INTERPOLATE
template <typename Range, typename Float>
auto interpolate(Range betas, Range rho, Float t){
  for (size_t i = 0; i < betas.size()-1; ++i){
    if (betas[i] <= t and t <= betas[i+1]){
      return (rho[i+1]-rho[i])/(betas[i+1]-betas[i])*(t-betas[i])+rho[i];
    }
  }
  return 0.0;
}
#endif 
