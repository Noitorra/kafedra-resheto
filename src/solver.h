#ifndef SOLVER_H
#define SOLVER_H

#include <mpi.h>

#include "header.h"
#include "cell.h"

namespace Solver {

template<typename T>
class Vector2<T> {
public:
  Vector2() {
    set( T(0), T(0) );
  }
  Vector2(const T& x, const T& y) {
    set( x, y );
  }

  void set(const T& x, const T& y) {
    m_x = x;
    m_y = y;
  }
  T x() { return m_x; }
  T y() { return m_y; }
protected:
  T m_x;
  T m_y;
};

typedef Vector2<double> Vector2d;
typedef Vector2<float> Vector2f;
typedef Vector2<int> Vector2i;
typedef Vector2<unsigned int> Vector2u;

struct SolverData {
  Vector2d H;
  Vector2u N;
  unsigned int VesselLength;
};

struct MPIData {

};


class Solver {
public:
  virtual void ReadConfig() {
    m_solverData.H.set(1.0, 1.0);
    m_solverData.N.set(50, 32);
    m_solverData.VesselLength = 10;
  }

  virtual void Init();
  virtual void Run();
  virtual void SaveData();
protected:
  Vector2u m_size;
  SolverData m_solverData;
  std::vector<Cell*> m_cells;
};

class SolverMPI {
public:

};

}

class Solver
{
public:
  unsigned int originNx;
  unsigned int originNy;
  unsigned int Nx;
  unsigned int Ny;
  unsigned int iter;
  unsigned int max_iter;
public:
    Solver();

    void Initialize();
    void Run();
    void Load(const std::string& filename);
    void Save(const std::string& filename);
    // var
    Cell* GetCell(unsigned int x, unsigned int y);
private:
    std::vector<Cell*> m_cell;

    int chBlock(unsigned int x, unsigned int y,
        unsigned int bx, unsigned int by,
        unsigned int bw, unsigned int bh);

    void saveMacroData();

    void syncCreateLeft();
    void syncCreateRight();

    void syncLeft();
    void syncRight();

    void syncCreate();
    void sync();

    void syncSaveMacro();
    void writeMacroData(std::vector< std::vector<double> >& data, int gas, int data_type);

    void makeIntegral(int gas1, int gas2);
};

#endif // SOLVER_H
