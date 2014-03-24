#ifndef SOLVER_H
#define SOLVER_H

#include <mpi.h>

#include "header.h"
#include "cell.h"
#include "types.h"

struct SolverData {
  Vector2d H;
  Vector2i N;
  unsigned int VesselLength;
};

struct MPIData {
  int size;
  int rank;
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
  Vector2i m_size;
  SolverData m_solverData;
  std::vector<Cell*> m_cells;
};

class SolverMPI {
public:

};

//class Solver
//{
//public:
//  unsigned int originNx;
//  unsigned int originNy;
//  unsigned int Nx;
//  unsigned int Ny;
//  unsigned int iter;
//  unsigned int max_iter;
//public:
//    Solver();

//    void Initialize();
//    void Run();
//    void Load(const std::string& filename);
//    void Save(const std::string& filename);
//    // var
//    Cell* GetCell(unsigned int x, unsigned int y);
//private:
//    std::vector<Cell*> m_cell;

//    int chBlock(unsigned int x, unsigned int y,
//        unsigned int bx, unsigned int by,
//        unsigned int bw, unsigned int bh);

//    void saveMacroData();

//    void syncCreateLeft();
//    void syncCreateRight();

//    void syncLeft();
//    void syncRight();

//    void syncCreate();
//    void sync();

//    void syncSaveMacro();
//    void writeMacroData(std::vector< std::vector<double> >& data, int gas, int data_type);

//    void makeIntegral(int gas1, int gas2);
//};

#endif // SOLVER_H
