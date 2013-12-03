#include "solver.h"

Solver::Solver()
{
    Nx = 100;
    Ny = 100;
}

void Solver::Initialize() {
    for(unsigned int i=0;i<Nx*Ny;i++) {
        m_cell.push_back(new Cell());
    }
}

void Solver::Load(const std::string& filename) {

}

void Solver::Save(const std::string &filename) {

}
