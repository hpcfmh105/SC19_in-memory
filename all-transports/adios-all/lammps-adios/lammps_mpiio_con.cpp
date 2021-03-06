/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "adios_read.h"
#include "adios_error.h"
//#include "adios_read_global.h"
#include "utility.h"
#include "msd-anal/run_msd.h"

#include "adios_adaptor.h"

#ifdef V_T
#include <VT.h>
int class_id;
int analysis_id;
#endif


        

#define DEBUG_Feng

int main (int argc, char ** argv) 
{
    /*
     * @input
     * @param NSTOP
     */
    if(argc !=2){
        printf("need to specify timstop\n");
        exit(-1);
    }
	int nsteps = atoi(argv[1]);
    double **msd;

    /******************** configuration stop ***********/
#ifdef ENABLE_TIMING
    double t0, t1, t2, t3, t4;
    double t_prepare, t_get, t_close, t_analy;
    t_prepare=0;
    t_get = 0;
    t_close = 0;
    t_analy = 0;
#endif




    int         rank, nprocs;
    MPI_Comm    comm = MPI_COMM_WORLD;
    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_DATASPACES;
    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_DIMES;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel = NULL;
    void * data = NULL;
    uint64_t start[2], count[2];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nprocs);

	double t_start = MPI_Wtime();
	// init msd
    uint64_t size_one = 5;

    msd =  init_msd(nsteps, size_one);

    /*
     * define the trace
     */

#ifdef V_T
     VT_classdef( "Analysis", &class_id );
     VT_funcdef("ANL", class_id, &analysis_id);
#endif
    
    int r;

    char *filepath = getenv("BP_DIR");
    if(filepath == NULL){
        fprintf(stderr, "scratch dir is not set!\n");
    }
    int timestep = -1;
    
    /**** use index file to keep track of current step *****/
    char step_index_file [256];
    int time_stamp = -1;// flag read from producer , this tells which is the largest version we can read
    int time_stamp_old=-1;
    sprintf(step_index_file, "%s/stamp.file", filepath);
	
    MPI_Barrier(comm);

    adios_read_init_method (method, comm, "verbose=3");
    PINF("rank %d:reader init\n", rank);

    char filename_atom[256];
    ADIOS_FILE * f;

    ADIOS_VARINFO * v;

    /* Using less readers to read the global array back, i.e., non-uniform */
    uint64_t slice_size;
    //for(timestep = 0; timestep < 10;){
    //
    // write will send -2 to mark end

    // flag whether producer finishes all the steps
    int has_more = 1;
    while(has_more){
        time_stamp_old = time_stamp;

        adios_adaptor_get_avail_version(comm, step_index_file, &time_stamp, nsteps);

        // has nothing yet
        if(time_stamp==-1){
            sleep(1);
            continue;
        }

        if(time_stamp == nsteps -1){
           has_more = 0;
        }

        // new step avaible from producer
        while(timestep < time_stamp){
            timestep++;
            if(rank ==0)
                PINF("----reader opened step %d\n", timestep);
            sprintf(filename_atom, "%s/atom_%d.bp", filepath, timestep);

            t0 = MPI_Wtime();
            f = adios_read_open_file(filename_atom, method, comm);
			size_t nelem = 1;

            // the first timestep will get variable metadata
            if(timestep == 0){
                 if (f == NULL)
                {
                    printf ("%s\n", adios_errmsg());
                    return -1;
                }
                if(S_OK != query_select_lammps( f, rank, nprocs, &sel, &nelem ) || sel == NULL){
                    PERR("query not succeed"); 
                    TRACE();
                    MPI_Abort(comm, -1);
                }
                

                data = malloc( nelem* sizeof (double));
                if (data == NULL)
                {
                    PERR("malloc failed.\n");
                    TRACE();
                    MPI_Abort(comm, -1);
                }

            }

            /* Read a subset of the temperature array */
            adios_schedule_read (f, sel, "var_2d_array", 0, 1, data);

            // timer for open and schedule
            t1 = MPI_Wtime();
            t_prepare+= t1-t0;
            
            // timer for actual read
            adios_perform_reads (f, 1);
            t2 = MPI_Wtime();
            t_get += t2-t1;

            // timer for closing file
            adios_read_close (f);
            t3 = MPI_Wtime();
            t_close += t3-t2;

            if(rank ==0)
                PINF("Step %d read\n", timestep);
            // analysis
#ifdef V_T
      VT_begin(analysis_id);
#endif
	calc_msd(msd, (double *)data, nelem/size_one, size_one, timestep);

            //run_analysis(data, slice_size, lp, sum_vx,sum_vy);
#ifdef V_T
      VT_end(analysis_id);
#endif

            t4 = MPI_Wtime();
            t_analy += t4-t3;
            PINF("rank %d: Step %d moments calculated,t_prepare %lf, t_read %lf, t_close %lf, t_analy %lf, time%lf\n", rank, timestep, t1-t0, t2-t1, t3-t2, t4-t3, t4);
        }
    }
    perform_msd_reduce(msd, nsteps, comm);
    free_msd(msd, size_one);

    free (data);
    adios_selection_delete (sel);
    //adios_free_varinfo (v);
	double t_end = MPI_Wtime();

	printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

#ifdef ENABLE_TIMING

    printf("[rank %d]:analysis_time %.3lf \n", rank, t_analy);
    MPI_Barrier(comm);
    t_end = MPI_Wtime();


        double global_t_prepare=0;
        double global_t_get=0;
        double global_t_close=0;
        double global_t_analy=0;
        MPI_Reduce(&t_prepare, &global_t_prepare, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&t_get, &global_t_get, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&t_close, &global_t_close, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
 
    if(rank == 0){
      PINF("stat:Consumer end  at %lf \n", t_end);
      PINF("stat:time for prepare %fs, read %f s; time for close %f s; max time for analy %f s\n",global_t_prepare/nprocs, global_t_get/nprocs, global_t_close/nprocs, global_t_analy);
    }
#endif
    

    MPI_Barrier (comm);
    adios_read_finalize_method (method);

    PINF("rank %d: exit\n", rank);

  /*
   * close logger
   */
    MPI_Finalize ();
    return 0;
}
