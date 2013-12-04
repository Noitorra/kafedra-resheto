#ifndef SPEED_H
#define SPEED_H

#include "header.h"
#include "cmath"

template<typename type>
inline type mod_vec(std::vector<type>& vec) {
  type result = type(0);
  std::vector<type>::iterator iter;
  for(iter=vec.begin();iter!=vec.end();iter++) {
    result += std::pow((*iter), 2);
  }
  return std::sqrt(result);
}

class Gas {
public:
    Gas(double getMass)
        :mass(getMass)
    {}
    double mass;

};

class Impulse {
public:
  Impulse(unsigned int size, double cutSpeed)
    :n(size),
     cut(cutSpeed) {
    double* line_impulse = new double[n];
    for(unsigned int i=0;i<n;i++) {
      line_impulse[i] = cut*(2.0*i/(n-1)-1);
      cout << "Inpulse[" << i << "] = " << line_impulse[i] << endl;
    }
    for(unsigned int x=0;x<n;x++)
      for(unsigned int y=0;y<n;y++)
        for(unsigned int z=0;z<n;z++) {
          std::vector<double> impulse3d;
          impulse3d.push_back(line_impulse[x]);
          impulse3d.push_back(line_impulse[y]);
          impulse3d.push_back(line_impulse[z]);
          if(mod_vec(impulse3d) < cut) {
            value.push_back(impulse3d);
          }
        }


  }
  ~Impulse() {}
  // Impulse access
  double cut;
  unsigned int n;
  std::vector<std::vector<double>> value;
  double dP;
  double d3P;
protected:

};

class Param
{
public:
    Param() {}
    static Param* instance() {
        static Param *s_Param = new Param();
        return s_Param;
    }


    // Impulse
    Impulse* impulse;
    // Global Program Params
    std::vector<Gas*> gas;
    double timestep;
};

#define P Param::instance()

#endif // SPEED_H
