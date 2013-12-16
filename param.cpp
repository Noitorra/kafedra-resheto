#include "param.h"

double mod_vec(std::vector<double>& vec) {
  double result = 0.0;
  std::vector<double>::iterator iter;
  for(iter=vec.begin();iter!=vec.end();iter++) {
    result += std::pow((*iter), 2);
  }
  return std::sqrt(result);
}
