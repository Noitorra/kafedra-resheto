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
    originNx = 80;
    originNy = 40;
    max_iter = 10;
    iter = 0;
}

int Solver::chBlock(unsigned int x, unsigned int y,
        unsigned int bx, unsigned int by,
        unsigned int bw, unsigned int bh) {
  if( (x>=bx && x<=(bx+bw)) && (y>=by && y<=(by+bh)) ) {
    if( x==bx || x==(bx+bw) || y==by || y==(by+bh) ) {
      // Block Edge
      return 2;
    } else {
      // Block Inside
      return 1;
    }
  } else {
    // Not in block
    return 0;
  }
}

void Solver::Initialize() {
  Nx = originNx/P->mpi_size;
  Ny = originNy;
  for(unsigned int y=0;y<Ny;y++) {
    for(unsigned int x=P->mpi_rank*Nx;x<(P->mpi_rank+1)*Nx;x++) {
      int block1 = chBlock(x, y, originNx*3/8, originNy*1/8,
                            originNx*4/8, originNy*2/8);
      int block2 = chBlock(x, y, originNx*1/8, originNy*5/8,
                            originNx*4/8, originNy*2/8);
      
      Cell* cell;
      if(block1 == 1 || block2 == 1) {
        cell = NULL;
      } else if(block1 == 2) {
        cell = new Cell();
        cell->m_h.push_back(1.0);
        cell->m_h.push_back(1.0);
        cell->gasIndex = 0;
        cell->T = 0.5;
      } else if(block2 == 2) {
         cell = new Cell();
        cell->m_h.push_back(1.0);
        cell->m_h.push_back(1.0);
        cell->gasIndex = 0;
        cell->T = 1.0;
      } else {
        cell = new Cell();
        cell->m_h.push_back(1.0);
        cell->m_h.push_back(1.0);
        cell->gasIndex = 0;
        cell->T = 0.5;
      }
      m_cell.push_back(cell);
    }
  }

  cout << "Rank=" << P->mpi_rank << " [Nx:Ny]=" << Nx << ":" << Ny << endl;
  for(unsigned int y=0;y<Ny;y++) {
    for(unsigned int x=0;x<Nx;x++) {
      if(GetCell(x,y))
        cout << 0;
      else
        cout << 1;
    }
    cout << endl;
  }

  //if(P->mpi_type == Param::MPI_LEFT)
  //  for(unsigned int y=0;y<Ny;y++) {
  //    GetCell(0, y)->T = 1.0;
  //  }

  for(unsigned int i=0;i<m_cell.size();i++) {
    if(m_cell[i]) m_cell[i]->Init();
  }

  //syncSaveMacro();
  
  for(unsigned int x=0;x<Nx;x++) {
    for(unsigned int y=0;y<Ny;y++) {
      if(GetCell(x, y)) {
        GetCell(x, y)->m_next[Cell::X] = GetCell(x+1, y);
        GetCell(x, y)->m_prev[Cell::X] = GetCell(x-1, y);
        GetCell(x, y)->m_next[Cell::Y] = GetCell(x, y+1);
        GetCell(x, y)->m_prev[Cell::Y] = GetCell(x, y-1);
      }
    }
  }

  for(unsigned int x=0;x<Nx;x++) {
    GetCell(x, 0)->m_prev[Cell::Y] = GetCell(x, Ny-1);
    GetCell(x, Ny-1)->m_next[Cell::Y] = GetCell(x, 0);
  }

  switch (P->mpi_type) {
    case Param::MPI_LEFT:
      for(unsigned int y=0;y<Ny;y++) {
        GetCell(0, y)->m_next[Cell::Y] = NULL;
        GetCell(0, y)->m_prev[Cell::Y] = NULL;
      }
      break;
    case Param::MPI_NORMAL:
      break;
    case Param::MPI_RIGHT:
      for(unsigned int y=0;y<Ny;y++) {
        GetCell(Nx-1, y)->m_next[Cell::Y] = NULL;
        GetCell(Nx-1, y)->m_prev[Cell::Y] = NULL;
      }
      break;
  }

  if(P->mpi_size > 1) {
    syncCreate();

    std::cout << "Init core [" << P->mpi_rank
    << "/" << P->mpi_size << "]"
    << " size=" << m_cell.size()
    << " type=" << P->mpi_type
    << std::endl;

    MPI::COMM_WORLD.Barrier();


  }
}

void Solver::syncCreateLeft() {
  int* recv = new int[Ny];
  MPI::COMM_WORLD.Recv(recv, Ny, MPI::INT, P->mpi_rank - 1, 10);
  for(int y=0;y<Ny;y++) {
    if(recv[y] == 1) {
      Cell* cell = new Cell();
      cell->InitBase();
      GetCell(0, y)->m_prev[Cell::X] = cell;
    }
  }

  int* send1 = new int[Ny];
  int* send2 = new int[Ny];
  for(int y=0;y<Ny;y++) {
    if(GetCell(0, y)) {
      if(GetCell(1, y)) {
        send2[y] = 1;
      } else {
        send2[y] = 0;
      }
      send1[y] = 1;
    } else {
      send1[y] = 0;
      send2[y] = 0;
    }
  }
  MPI::COMM_WORLD.Ssend(send1, Ny, MPI::INT, P->mpi_rank - 1, 10);
  MPI::COMM_WORLD.Ssend(send2, Ny, MPI::INT, P->mpi_rank - 1, 10);
}

void Solver::syncCreateRight() {
  int* send = new int[Ny];
  for(int y=0;y<Ny;y++) {
    if(GetCell(Nx-1, y)) {
      send[y] = 1;
    } else {
      send[y] = 0;
    }
  }
  MPI::COMM_WORLD.Ssend(send, Ny, MPI::INT, P->mpi_rank + 1, 10);

  int* recv1 = new int[Ny];
  int* recv2 = new int[Ny];
  MPI::COMM_WORLD.Recv(recv1, Ny, MPI::INT, P->mpi_rank + 1, 10);
  MPI::COMM_WORLD.Recv(recv2, Ny, MPI::INT, P->mpi_rank + 1, 10);
  for(int y=0;y<Ny;y++) {
    if(recv1[y] == 1) {
      Cell* cell1 = new Cell();
      cell1->InitBase();
      GetCell(Nx-1, y)->m_next[Cell::X] = cell1;
      if(recv2[y] == 1) {
        Cell* cell2 = new Cell();
        cell2->InitBase();
        cell1->m_next[Cell::X] = cell2;
      }
    }
  }
}

void Solver::syncLeft() {
  for(int y=0;y<Ny;y++) {
    Cell* cell = GetCell(0, y);
    if(cell) {
      cell = cell->m_prev[Cell::X];
      if(cell) {
        MPI::COMM_WORLD.Recv(cell->m_half, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank - 1, 20);
        MPI::COMM_WORLD.Recv(cell->m_value, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank - 1, 20);
      }
    }
  }

  for(int y=0;y<Ny;y++) {
    Cell* cell1 = GetCell(0, y);
    if(cell1) {
      MPI::COMM_WORLD.Ssend(cell1->m_half, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank - 1, 20);
      MPI::COMM_WORLD.Ssend(cell1->m_value, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank - 1, 20);
      Cell* cell2 = GetCell(1, y);
      if(cell2) {
        MPI::COMM_WORLD.Ssend(cell2->m_half, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank - 1, 20);
        MPI::COMM_WORLD.Ssend(cell2->m_value, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank - 1, 20);
      }
    }
  }
}

void Solver::syncRight() {
  for(int y=0;y<Ny;y++) {
    Cell* cell = GetCell(Nx-1, y);
    if(cell) {
      MPI::COMM_WORLD.Ssend(cell->m_half, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank + 1, 20);
      MPI::COMM_WORLD.Ssend(cell->m_value, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank + 1, 20);
    }
  }

  for(int y=0;y<Ny;y++) {
    Cell* cell1 = GetCell(Nx-1, y);
    if(cell1) {
      cell1 = cell1->m_next[Cell::X];
      if(cell1) {
        MPI::COMM_WORLD.Recv(cell1->m_half, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank + 1, 20);
        MPI::COMM_WORLD.Recv(cell1->m_value, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank + 1, 20);
        Cell* cell2 = cell1->m_next[Cell::X];
        if(cell2) {
          MPI::COMM_WORLD.Recv(cell2->m_half, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank + 1, 20);
          MPI::COMM_WORLD.Recv(cell2->m_value, P->impulse->value.size(), MPI::DOUBLE, P->mpi_rank + 1, 20);
        }
      }
    }
  }
}

void Solver::syncCreate() {
  MPI::COMM_WORLD.Barrier();
  switch(P->mpi_type) {
    case Param::MPI_LEFT:
      syncCreateRight();
      break;
    case Param::MPI_NORMAL:
      syncCreateLeft();
      syncCreateRight();
      break;
    case Param::MPI_RIGHT:
      syncCreateLeft();
      break;
  }
  MPI::COMM_WORLD.Barrier();
}

void Solver::sync() {
  MPI::COMM_WORLD.Barrier();
  switch(P->mpi_type) {
    case Param::MPI_LEFT:
      syncRight();
      break;
    case Param::MPI_NORMAL:
      syncLeft();
      syncRight();
      break;
    case Param::MPI_RIGHT:
      syncLeft();
      break;
  }
  MPI::COMM_WORLD.Barrier();
}

void Solver::Load(const std::string& filename) {

}

void Solver::Save(const std::string &filename) {

}

Cell* Solver::GetCell(unsigned int x, unsigned int y) {
  if(0<=x && x<Nx && 0<=y && y<Ny)
    return m_cell[x + y*Nx];
  else
    return NULL;
}

void Solver::Run() {
  for(iter=0;iter<max_iter;iter++) {
    MPI::COMM_WORLD.Barrier();

    sync();
    for(int i=0;i<m_cell.size();i++) {
      if(m_cell[i]) m_cell[i]->ComputeHalf(Cell::X);
    }

    sync();
    for(unsigned int i=0;i<Nx*Ny;i++) {
      if(m_cell[i]) m_cell[i]->ComputeValue(Cell::X);
    }

    for(unsigned int i=0;i<Nx*Ny;i++) {
      if(m_cell[i]) m_cell[i]->ComputeHalf(Cell::Y);
    }
    for(unsigned int i=0;i<Nx*Ny;i++) {
      if(m_cell[i]) m_cell[i]->ComputeValue(Cell::Y);
    }

    //saveMacroData();
    syncSaveMacro();
    MPI::COMM_WORLD.Barrier();
    if(P->mpi_rank==0) {
      cout << "Run [" << 1.0*(iter+1)/max_iter*100 << "%]" << std::endl;;
    }
  }
  if(P->mpi_rank==0) {
    cout << "Done" << std::endl;
  }
}

void Solver::syncSaveMacro() {
  double* density = new double[m_cell.size()];
  for(unsigned int i=0;i<m_cell.size();i++) {
    if(m_cell[i]) {
      density[i] = m_cell[i]->getDensity();
    } else {
      //cout << "OPA:" << i <<endl;
      density[i] = 0.0;
    }
    //cout << density[i];
  }
  if(P->mpi_rank == 0) {
    std::vector< std::vector<double> > vvden;
    vvden.resize(Ny);
    for(int rank=0;rank<P->mpi_size;rank++) {
      std::vector<double> loc_vec;
      if(rank==0) {
        // interesting !!
        for(unsigned int y=0;y<Ny;y++) {
          for(unsigned int x=0;x<Nx;x++) {
            vvden[y].push_back(density[y*Nx+x]);
          }
        }
      } else {
        int rNx;
        MPI::COMM_WORLD.Recv(&rNx, 1, MPI::INT, rank, 30);
        double* denr = new double[rNx*Ny];
        MPI::COMM_WORLD.Recv(denr, rNx*Ny, MPI::DOUBLE, rank, 31);
        // same interest !!
         for(unsigned int y=0;y<Ny;y++) {
          for(unsigned int x=0;x<rNx;x++) {
            vvden[y].push_back(denr[y*rNx+x]);
          }
        }
        delete []denr;
      }
    }
    // writing whole thing to file
    writeMacroData(vvden);

  } else {
    int sNx = Nx;
    MPI::COMM_WORLD.Send(&sNx, 1, MPI::INT, 0, 30);
    MPI::COMM_WORLD.Send(density, m_cell.size(), MPI::DOUBLE, 0, 31);
  }
  delete []density;
}

void Solver::writeMacroData(std::vector< std::vector<double> >& data) {
  std::string filename;
  filename = "data/Den/" + ToString(iter) + ".bin";

  std::ofstream fs(filename.c_str(), std::ios::out | std::ios::binary);
  double density;
  for(unsigned int y=0;y<originNy;y++) {
    for(unsigned int x=0;x<originNx;x++) {
      //cout << data[y][x];
      fs.write(reinterpret_cast<const char*>(&data[y][x]), sizeof(double));
    }
  }
  fs.close();
}

void Solver::saveMacroData() {
  std::string filename;
  filename = "data/Den/" + ToString(iter) + ".bin";

  std::ofstream fs(filename.c_str(), std::ios::out | std::ios::binary);
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

