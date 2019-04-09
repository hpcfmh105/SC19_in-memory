#pragma once
#include <mpi.h>
#include <string>

namespace lammpsBridge
{
  void Initialize(MPI_Comm world,int step, int nlocal, double **x, 
               double xsublo, double xsubhi, 
               double ysublo, double ysubhi, 
               double zsublo, double zsubhi,
               const std::string& config_file);
  void update_bridge(int step, double *buffer);
  void Finalize();
}

