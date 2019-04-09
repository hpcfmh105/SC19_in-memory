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
static transport_method_t transport;

#ifdef V_T
#include <VT.h>
int class_id;
int analysis_id;
#endif


static char var_name[STRING_LENGTH];
static char var_size[STRING_LENGTH];
static size_t elem_size=sizeof(double);

//#define DEBUG_Feng


int main (int argc, char ** argv) 
{
     /*
     * @input
     * @param NSTOP
     * @param FILESIZE2PRODUCE
     */
  if(argc !=3)
  {
    printf("heat3d nsteps\n");
    exit(-1);
  }

  int nsteps = atoi(argv[1]);
  int slice_size = atoi(argv[2]);

  int lp = N_LP;

  int iconf[7];
  double conf[2];
  uint64_t nx, ny, nodeedge, nblock, nproc;
  uint64_t nn = 0;
  double sum_vx[NMOMENT], sum_vy[NMOMENT];

  /******************** configuration stop ***********/

  int rank, nprocs;
  int has_keep=0;
  double * data = NULL;
  uint64_t start[3], count[3];
  MPI_Init (&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &nprocs);

  if(rank==0)
    readParam(iconf, conf);

  /* Broadcast input parameters */
  MPI_Bcast(iconf,7,MPI_INT,0,comm);
  MPI_Bcast(conf,2,MPI_DOUBLE,0,comm);

  /* Assign input parameters to variables */
  nodeedge = iconf[0];
  nodeedge = iconf[1];
  nblock = nprocs;
  int maxStep   = iconf[6];
  nx = nodeedge * nblock;
  ny = nodeedge;
  nn = nx * ny;
  char nodename[256];
  int nodename_length;
  MPI_Get_processor_name(nodename, &nodename_length );

  /*
   * define the trace
   */

  int r;
  char *filepath = getenv("BP_DIR");
  if(filepath == NULL)
  {
    fprintf(stderr, "scratch dir is not set!\n");
  }

    /*
     * get transport method
     */
  if(rank == 0)
  {
    PINF("stat: Consumer start at %lf \n", MPI_Wtime());
  }
  char filename[256];
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
  // how many lines in global
  //uint64_t gdims[2] = {2, global_size};
  //dspaces_define_gdim(var_name, 2,gdims);
  uint64_t bounds_size[6] = {0};
  uint64_t bounds[6] = {0};
  double time_comm;
  uint64_t* global_size = (uint64_t *) malloc(sizeof(long));

    /* Using less readers to read the global array back, i.e., non-uniform */
   // uint64_t slice_size = *global_size/nprocs;
   // if (rank == nprocs-1) /* last rank may read more lines */
    //    slice_size = slice_size + *global_size%nprocs;
  slice_size = slice_size * nodeedge;
  bounds[0] = slice_size * rank;
  bounds[3] = bounds[0] + slice_size - 1;
  bounds[4] = nodeedge - 1;
  data = (double *)malloc (slice_size * nodeedge * sizeof(double));
  if (data == NULL)
  {
    size_t allc_size=slice_size * sizeof (double) ;
    PERR("malloc failed with %ld bytes", allc_size);
    TRACE();
    MPI_Abort(comm, -1);
  }
  int timestep = 0;

  while(timestep < nsteps)
  {
    PINF("rank %d: Step %d start\n", rank, timestep);
           /* Read a subset of the temperature array */
        // 0:not used for strea; 1: must be set in stream
        // analysis
    run_analysis(data, slice_size * (nodeedge)/SIZE_ONE, lp, sum_vx,sum_vy);

    t4 = MPI_Wtime();
    t_analy += t4-t3;

    PDBG("previous step released");
    // advance to (1)the next availibale step (2)blocked if not unavailble
    PDBG("successfully step into next available step");
    if(rank == 0)
      PINF("rank %d: Step %d moments calculated, t_read %lf, t_advance %lf, t_analy %lf\n", rank, timestep, t2-t1, t3-t2, t4-t3);
    timestep ++;
  }

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
  free (data);
  free (global_size);

  /*
   * close logger
   */
  MPI_Finalize ();
  printf("rank %d: exit\n", rank);

  return 0;
}
