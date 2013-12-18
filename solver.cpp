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
  m_cell.resize(Nx*Ny, NULL);

  for(unsigned int x=0;x<Nx;x++) {
    for(unsigned int y=0;y<Ny;y++) {
      // Block #1
      if( (Nx*2/6 < x && x < Nx*5/6) && (Ny*1/8 < y && y < Ny*3/8) ) {

      } else
      // Block #2
      if( (Nx*1/6 < x && x < Nx*4/6) && (Ny*5/8 < y && y < Ny*7/8) ) {

      } else
      // InnerSpace
      {
        Cell* cell = new Cell();
        cell->m_h.push_back(1.0);
        cell->m_h.push_back(1.0);
        cell->gasIndex = 0;
        // Upper and lower lines
        if ( (x==0 || x==Nx-1) && (0 < y && y < Ny-1) ) {
          cell->T = 1.0;
        } else 
        // Left and Right lines
        if( (y==0 || y==Ny-1) && (0 < x && x < Nx-1) ) {
          cell->T = 1.0;
        } else
        // Normal Cells
        {
          cell->T = 0.5;
        }
        m_cell[x+y*Nx] = cell;
      }
    }
  }


  for(unsigned int i=0;i<Nx*Ny;i++) {
    std::cout << "\r";
    if(m_cell[i]) m_cell[i]->Init();
    std::cout << "Initialization [" << i+1 << "/" << Nx*Ny << "]";
  }
  std::cout << std::endl;

  for(unsigned int x=1;x<Nx-1;x++) {
    for(unsigned int y=1;y<Ny-1;y++) {
      if(GetCell(x, y)) {
        GetCell(x, y)->m_next[Cell::X] = GetCell(x+1, y);
        GetCell(x, y)->m_prev[Cell::X] = GetCell(x-1, y);
        GetCell(x, y)->m_next[Cell::Y] = GetCell(x, y+1);
        GetCell(x, y)->m_prev[Cell::Y] = GetCell(x, y-1);
      }
    }
  }

  for(unsigned int x=1;x<Nx-1;x++) {
    GetCell(x, 0)->m_next[Cell::Y] = GetCell(x, 1);
    GetCell(x, Ny-1)->m_prev[Cell::Y] = GetCell(x, Ny-2);
  }

  for(unsigned int y=1;y<Ny-1;y++) {
    GetCell(0, y)->m_next[Cell::X] = GetCell(1, y);
    GetCell(Nx-1, y)->m_prev[Cell::X] = GetCell(Nx-2, y);
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
      if(m_cell[i]) m_cell[i]->ComputeHalf(Cell::X);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      if(m_cell[i]) m_cell[i]->ComputeValue(Cell::X);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      if(m_cell[i]) m_cell[i]->ComputeHalf(Cell::Y);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      if(m_cell[i]) m_cell[i]->ComputeValue(Cell::Y);
    }
          
    saveMacroData();

    std::cout << "Run [" << iter+1 << "/" << max_iter << "]";
  }
  std::cout << std::endl;
}

void Solver::saveMacroData() {
  std::string filename;
  filename = "data/Den/" + ToString(iter) + ".bin";

  std::ofstream fs(filename, std::ios::out | std::ios::binary);
  double density;
  for(unsigned int i=0;i<Nx*Ny;i++) {
    if(m_cell[i]) {
      density = m_cell[i]->getDensity();
    } else {
      density = 0.0;
    }
    fs.write(reinterpret_cast<const char*>(&density), sizeof(double));
  }
  fs.close();
}

