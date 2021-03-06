#ifndef SPEED_H
#define SPEED_H

#include "header.h"
#include "cmath"

double mod_vec(std::vector<double>& vec);

class Gas {
public:
    Gas(double getMass)
        :mass(getMass)
    {}
    double mass;
};

class Impulse {
public:
  // Impulse access
  double cut;
  unsigned int n;
  std::vector< std::vector<double> > value;
  double dP;
  double d3P;
  int*** xyz2i;
public:
  Impulse(unsigned int size, double cutSpeed)
    :n(size),
     cut(cutSpeed)
  {
    dP = 2*cut/(n-1);
    d3P = std::pow(dP, 3);
    xyz2i = new int**[n];
    for(int x=0;x<n;x++) {
      xyz2i[x] = new int*[n];
      for(int y=0;y<n;y++) {
        xyz2i[x][y] = new int[n];
        for(int z=0;z<n;z++) {
          xyz2i[x][y][z] = -1;
        }
      }
    }
    double* line_impulse = new double[n];
    for(unsigned int i=0;i<n;i++) {
      //cut*(2.0*i/(n-1)-1);
      line_impulse[i] = dP*i - cut;
      //cout << "Inpulse[" << i << "] = " << line_impulse[i] << endl;
    }
    for(unsigned int x=0;x<n;x++)
      for(unsigned int y=0;y<n;y++)
        for(unsigned int z=0;z<n;z++) {
          std::vector<double> impulse3d;
          impulse3d.push_back(line_impulse[x]);
          impulse3d.push_back(line_impulse[y]);
          impulse3d.push_back(line_impulse[z]);
          if(mod_vec(impulse3d) < cut) {
            xyz2i[x][y][z] = value.size();
            value.push_back(impulse3d);
          }
        }


  }
  ~Impulse() {}
protected:

};

class Param {
public:
  enum MPI_Type { MPI_LEFT, MPI_NORMAL, MPI_RIGHT };
    Param() {}
    static Param* instance() {
        static Param *s_Param = new Param();
        return s_Param;
    }
    // MPI Info
    int mpi_rank;
    int mpi_size;
    MPI_Type mpi_type;

    // Impulse
    Impulse* impulse;
    // Global Program Params
    std::vector<Gas*> gas;
    double timestep;
};

#define P Param::instance()

#endif // SPEED_H
