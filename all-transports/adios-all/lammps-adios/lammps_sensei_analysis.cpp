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

/*#define CLOG_MAIN*/
//#include "utility.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "adios_adaptor.h"
//#include "adios_helper.h"
#include "adios.h"
#include <assert.h>
#include "transports.h"
#include "utility.h"
#include "msd-anal/run_msd.h"
#include "ds_adaptor.h"
#include "adios_adaptor.h"

#include "lammps_sensei_analysis.h"

#define SIZE_ONE (5)

static char var_name[STRING_LENGTH];
static char var_size[STRING_LENGTH];
static size_t elem_size=sizeof(double);
#include "lammps_sensei_analysis.h"

#define SIZE_ONE (5)

/*double **init_analyze(int nsteps, int size_one)
{
  double **msd;
  size_one = SIZE_ONE;
  msd =  init_msd(nsteps, size_one);
  return msd;
}*/

/*void analyze(int nsteps, double *data,  int size_one)
{
  double **msd;
  //size_one = SIZE_ONE;
  msd =  init_msd(nsteps, size_one);

  int ii = 0;
  sizes[0] = nprocs;
  sizes[1] = size_one;
  if(tmplocal == 0)
    sizes[2] = 512000;
  else
    sizes[2] = tmplocal;
  data = (double *)malloc (slice_size * sizes[2] * SIZE_ONE* sizeof (double));

  while(ii < nsteps){
    size_t nelem;
    uint64_t bounds[6] = {0};
                                        // xmin
    bounds[2] = 0;
    bounds[0]=slice_size * rank;
    bounds[1]= 0;
                                        // ymin
    bounds[5] = sizes[2] - 1;
                                        // xmax
    bounds[3]=slice_size * (rank + 1) -1;
                                        // ymax
    bounds[4]= EVAL(size_one)-1;


    if(S_OK != get_common_buffer(transport_minor, ii,3, bounds,rank, var_name, (void **)&data, elem_size, &time_comm, &comm)){
        TRACE();
        printf("ERROR, after get_common_buffer in main rank %d \n", rank);
        MPI_Abort(comm, -1);
        } 
    printf("after get common buffer step: %d \n", ii);
    t1 =MPI_Wtime();
    printf("before calc_msd %ld, %ld, %d \n",  sizes[2], size_one, ii);
    calc_msd(msd, data, slice_size * sizes[2], size_one, ii);
    t2 =MPI_Wtime();
    t_analy += t2-t1;
    ii++;
    }
    perform_msd_reduce(msd, nsteps, comm);
    free_msd(msd, size_one);

    double t_end = MPI_Wtime();
    double global_t_analy;
    
    //output

    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);
    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
        printf("[msd]: max t_analysis %.3fs, now existing\n", global_t_analy);
    }

    MPI_Finalize ();
    free(sizes);

}*/

void analyze(int nsteps, int ii, double *data,  int size_one, /*uint64_t bounds[6],*/ int nlocal, double **msd)
{
  /**In this function, the nlocal means the sizes[2]**/
  //double **msd;
  // msd =  init_msd(nsteps, size_one);
  //int ii = 0;
  //sizes[0] = nprocs;
  //sizes[1] = size_one;
  //if(tmplocal == 0)
    //sizes[2] = 512000;
  //else
    //sizes[2] = tmplocal;
  //data = (double *)malloc (nlocal * size_one * sizeof (double));

  //while(ii < nsteps){
    //size_t nelem;
    //uint64_t bounds[6] = {0};
                                        // xmin
    //bounds[2] = 0;
    //bounds[0]=slice_size * rank;
    //bounds[1]= 0;
                                        // ymin
    //bounds[5] = sizes[2] - 1;
                                        // xmax
    //bounds[3]=slice_size * (rank + 1) -1;
                                        // ymax
    //bounds[4]= EVAL(size_one)-1;


    /*if(S_OK != get_common_buffer(transport_minor, ii,3, bounds,rank, var_name, (void **)&data, elem_size, &time_comm, &comm)){
        TRACE();
        printf("ERROR, after get_common_buffer in main rank %d \n", rank);
        MPI_Abort(comm, -1);
        }*/ 
    double t1, t2, t3, t4;
    double t_read_1, t_read_2, t_analy;
    t_read_1 = 0;
    t_read_2 = 0;
    t_analy = 0;

    printf("after get common buffer step: %d \n", ii);
    t1 =MPI_Wtime();
    printf("before calc_msd %ld, %ld, %d \n",  nlocal, size_one, ii);
    calc_msd(msd, data, nlocal, size_one, ii);
    t2 =MPI_Wtime();
    t_analy += t2-t1;
    //ii++;
    //}
    //perform_msd_reduce(msd, nsteps, comm);
    //free_msd(msd, size_one);

    /*double t_end = MPI_Wtime();
    double global_t_analy;
    //output

    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);
    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
        printf("[msd]: max t_analysis %.3fs, now existing\n", global_t_analy);
    }*/

    //MPI_Finalize ();
   // free(sizes);

}

/*void finalize_analyze(int nsteps, int size_one, double **msd)
{ 
    perform_msd_reduce(msd, nsteps, comm);
    free_msd(msd, size_one);

    double t_end = MPI_Wtime();
    double global_t_analy;
 
    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);  
    printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);
    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
        printf("[msd]: max t_analysis %.3fs, now existing\n", global_t_analy);
    }

    MPI_Finalize ();
    free(sizes);
}*/

