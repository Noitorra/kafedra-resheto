#ifndef SOLVER_H
#define SOLVER_H

#include <mpi.h>

#include "header.h"
#include "cell.h"

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
    void writeMacroData(std::vector< std::vector<double> >& data);
};

#endif // SOLVER_H
