#include <iostream>
#include "solver.h"
#include "param.h"

using namespace std;

int main()
{
  P->impulse = new Impulse(20, 4.8);
  P->timestep = 0.1;
  P->gas.push_back(new Gas(1.0));
  P->gas.push_back(new Gas(0.5));

  Solver solver;
  solver.max_iter = 500;
  solver.Initialize();
  solver.Run();
  return 0;
}

