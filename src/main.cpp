#include <iostream>
#include "solver.h"
#include "param.h"

using namespace std;

int main(int argc, char* argv[])
{
//  double start_time, end_time;
//  start_time = MPI::Wtime();

  MPI::Init(argc, argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  std::cout << "Processor: " << rank << " of " << size << std::endl;

  if( size > 1 ) {
    std::cout << "MPI avaliable, using SolverMPI." << std::endl;
  } else {
    std::cout << "No MPI avaliable, using Solver." << std::endl;
  }



  MPI::Finalize();

//  P->mpi_size = MPI::COMM_WORLD.Get_size();
//  P->mpi_rank = MPI::COMM_WORLD.Get_rank();

//  if(P->mpi_rank == 0)
//    P->mpi_type = Param::MPI_LEFT;
//  else if(P->mpi_rank == P->mpi_size - 1)
//    P->mpi_type = Param::MPI_RIGHT;
//  else
//    P->mpi_type = Param::MPI_NORMAL;

//  P->impulse = new Impulse(20, 4.8);
//  P->timestep = 0.01;
//  P->gas.push_back(new Gas(1.0));
//  P->gas.push_back(new Gas(0.5));

//  Solver solver;
//  solver.originNx = 120;
//  solver.originNy = 60;
//  solver.max_iter = 1000;
//  solver.Initialize();
//  solver.Run();

//  MPI::Finalize();
//  end_time = MPI::Wtime();
//  std::cout << "ExTime: " << end_time - start_time << "s" << std::endl;
  return 0;
}

