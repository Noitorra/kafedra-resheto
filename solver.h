#ifndef SOLVER_H
#define SOLVER_H

#include "header.h"
#include "cell.h"

class Solver
{
public:
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
    void saveMacroData();
};

#endif // SOLVER_H
