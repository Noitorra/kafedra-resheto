#ifndef SOLVER_H
#define SOLVER_H

#include "header.h"
#include "cell.h"

class Solver
{
public:
    unsigned int Nx;
    unsigned int Ny;
public:
    Solver();

    void Initialize();
    void Load(const std::string& filename);
    void Save(const std::string& filename);
    // var
private:
    std::vector<Cell*> m_cell;
};

#endif // SOLVER_H
