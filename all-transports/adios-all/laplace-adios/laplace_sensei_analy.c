/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.  *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
 * This code is using the generic read API, which can read in
 * arbitrary slices of an array and thus we can read in an array
 * on arbitrary number of processes (provided our code is smart 
 * enough to do the domain decomposition).
 *
 * Run this example after adios_global, which generates 
 * adios_global.bp. Run this example on equal or less 
 * number of processes since we decompose only on one 
 * dimension of the global array here. 
*/

#define CLOG_MAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "ds_adaptor.h"
#include <assert.h>
#include "transports.h"

#include "utility.h"
#include "nmoments-anal/run_analysis.h"

static char var_name[STRING_LENGTH];
static char var_size[STRING_LENGTH];
static size_t elem_size=sizeof(double);

void analyze(int nsteps/*,int lp*/, double *data, uint64_t nodeedge, int rank, MPI_Comm comm, int timestep)
{
  int lp = N_LP;
  uint64_t nn = 0;
  double sum_vx[NMOMENT], sum_vy[NMOMENT];
  /*char nodename[256];
  int nodename_length;
  MPI_Get_processor_name(nodename, &nodename_length );*/

  if(rank == 0)
  {
    PINF("stat: Consumer start at %lf \n", MPI_Wtime());
  }
  double t_start = MPI_Wtime();
  /******************** configuration stop ***********/
#ifdef ENABLE_TIMING
  double t0, t1, t2, t3, t4;
  double t_prepare, t_get, t_close, t_analy;
  t_prepare=0;
  t_get = 0;
  t_close = 0;
  t_analy = 0;
#endif

  /*
   * set bounds and dspaces variables
   */
  sprintf(var_name, "atom");
  sprintf(var_size, "atom_size");
  if (data == NULL)
  {
    size_t allc_size= sizeof (double) ;
    PERR("malloc failed with %ld bytes", allc_size);
    TRACE();
    MPI_Abort(comm, -1);
  }
  run_analysis(data, nodeedge/SIZE_ONE, lp, sum_vx,sum_vy);

  t4 = MPI_Wtime();
  t_analy += t4-t3;
  PDBG("previous step released");
    // advance to (1)the next availibale step (2)blocked if not unavailble
  PDBG("successfully step into next available step");
  if(rank == 0)
    PINF("rank %d: Step %d moments calculated, t_read %lf, t_advance %lf, t_analy %lf\n", rank, timestep, t2-t1, t3-t2, t4-t3);

#ifdef ENABLE_TIMING
  printf("[rank %d]:analysis_time %.3lf \n", rank, t_analy);
  MPI_Barrier(comm);
  double t_end = MPI_Wtime();
  double global_t_analy=0;
  MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  if(rank == 0)
  {
    PINF("stat:Consumer end  at %lf \n", t_end);
    PINF("stat:max time for analyst %f s\n",global_t_analy);
  }
#endif 
  printf("rank %d: exit\n", rank);

  //return 0;
}

