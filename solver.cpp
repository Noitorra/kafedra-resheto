#include "solver.h"
#include <fstream>
#include <sstream>

template <typename T>
std::string ToString(T x) {
  std::ostringstream os;
  os << x;
  return os.str();
}

Solver::Solver()
{
    Nx = 50;
    Ny = 50;
    max_iter = 10;
    iter = 0;
}

void Solver::Initialize() {
  for(unsigned int i=0;i<Nx*Ny;i++) {
    std::cout << "\r";
    Cell* cell = new Cell();
    cell->T = 1.0;
    cell->m_h.push_back(1.0);
    cell->m_h.push_back(1.0);
    cell->gasIndex = 0;


    m_cell.push_back(cell);

    std::cout << "Initialization [" << i+1 << "/" << Nx*Ny << "]";
  }
  std::cout << std::endl;

  for(unsigned int x=1;x<Nx-1;x++) {
    for(unsigned int y=1;y<Ny-1;y++) {
      GetCell(x, y)->m_next[Cell::X] = GetCell(x+1, y);
      GetCell(x, y)->m_prev[Cell::X] = GetCell(x-1, y);
      GetCell(x, y)->m_next[Cell::Y] = GetCell(x, y+1);
      GetCell(x, y)->m_prev[Cell::Y] = GetCell(x, y-1);

      GetCell(x, y)->T = 0.5;
      GetCell(x, y)->Init();
    }
  }

  for(unsigned int x=1;x<Nx-1;x++) {
    GetCell(x, 0)->m_next[Cell::Y] = GetCell(x, 1);
    GetCell(x, Ny-1)->m_prev[Cell::Y] = GetCell(x, Ny-2);

    GetCell(x, 0)->Init();
    GetCell(x, Ny-1)->Init();
  }

  for(unsigned int y=1;y<Ny-1;y++) {
    GetCell(0, y)->m_next[Cell::X] = GetCell(1, y);
    GetCell(Nx-1, y)->m_prev[Cell::X] = GetCell(Nx-2, y);

    GetCell(0, y)->Init();
    GetCell(Nx-1, y)->Init();
  }

  //saveMacroData();
}

void Solver::Load(const std::string& filename) {

}

void Solver::Save(const std::string &filename) {

}

Cell* Solver::GetCell(unsigned int x, unsigned int y) {
  return m_cell[x + y*Nx];
}

void Solver::Run() {
  for(iter=0;iter<max_iter;iter++) {
    std::cout << "\r";
    for(unsigned int i=0;i<Nx*Ny;i++) {
      m_cell[i]->ComputeHalf(Cell::X);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      m_cell[i]->ComputeValue(Cell::X);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      m_cell[i]->ComputeHalf(Cell::Y);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      m_cell[i]->ComputeValue(Cell::Y);
    }
          
    saveMacroData();

    std::cout << "Run [" << iter+1 << "/" << max_iter << "]";
  }
  std::cout << std::endl;
}

void Solver::saveMacroData() {
  std::string filename;
  filename = "data/Den" + ToString(iter) + ".bin";

  std::ofstream fs(filename, std::ios::out | std::ios::binary);
  double* density = new double[Nx*Ny];
  for(unsigned int i=0;i<Nx*Ny;i++) {
    density[i] = m_cell[i]->getDensity();
    //std::cout << i << " on "<< density[i] << std::endl;
    fs.write(reinterpret_cast<const char*>(&density[i]), sizeof(double));
  }
  
  //for(unsigned int i=0;i<Nx*Ny;i++) {
  //  filestream << ToString(i) << " " << ToString(m_cell[i]->getDensity()) << std::endl;
  //}
  fs.close();
}

