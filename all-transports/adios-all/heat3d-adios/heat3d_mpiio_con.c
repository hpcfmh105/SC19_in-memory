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

#include "adios_read.h"
#include "adios_error.h"
#include "adios_adaptor.h"
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



//#define DEBUG_Feng


int main (int argc, char ** argv) 
{
     /*
     * @input
     * @param NSTOP
     * @param FILESIZE2PRODUCE
     */
    if(argc !=2){
        printf("need to specify timstop\n");
        exit(-1);
    }
    int nstop = atoi(argv[1]);

    int lp = N_LP;

    double sum_vx[NMOMENT], sum_vy[NMOMENT];

    /******************** configuration stop ***********/
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

    int has_keep=0;

    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_DIMES;
    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;

    ADIOS_SELECTION * sel;
    void * data = NULL;
    uint64_t start[3], count[3];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nprocs);

    char nodename[256];
    int nodename_length;
    MPI_Get_processor_name(nodename, &nodename_length );

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

   int has_more = 1;
    while(has_more){
        time_stamp_old = time_stamp;

        adios_adaptor_get_avail_version(comm, step_index_file, &time_stamp, nstop);

        // has nothing yet
        if(time_stamp==-1){
            sleep(1);
            continue;
        }

        if(time_stamp == nstop -1){
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
			uint64_t size_tot = 0;

            // the first timestep will get variable metadata
            if(timestep == 0){
                 if (f == NULL)
                {
                    printf ("%s\n", adios_errmsg());
                    return -1;
                }


                ADIOS_VARINFO * v = adios_inq_var (f, "var_3d_array");

                if(rank ==0)
                PINF("reader opened the stream, dims = %ld, %ld, %ld\n", v->dims[0], v->dims[1], v->dims[2]);


                /* Using less readers to read the global array back, i.e., non-uniform */
                uint64_t slice_size = v->dims[0]/nprocs;
                start[0] = slice_size * rank;
                if (rank == nprocs-1) /* last rank may read more lines */
                    slice_size = slice_size + v->dims[0]%nprocs;
                count[0] = slice_size;

                start[1] = 0;
                count[1] = v->dims[1];
                start[2] = 0;
                count[2] = v->dims[2];
                size_tot = slice_size * v->dims[1] * v->dims[2];
                PINF("rank %d: start: (%ld, %ld, %ld), count:( %ld, %ld, %ld)\n", rank, start[0], start[1], start[2], count[0], count[1], count[2]);


                data = malloc (slice_size * v->dims[1] * v->dims[2] * sizeof (double));
                if (data == NULL)
                {

                    size_t allc_size=slice_size * v->dims[1] * v->dims[2] * sizeof (double);
                    PERR("malloc failed with %ld bytes", allc_size);
                    MPI_Abort(comm, -1);
                }

                sel = adios_selection_boundingbox (v->ndim, start, count);

            }

            /* Read a subset of the temperature array */
            adios_schedule_read (f, sel, "var_3d_array", 0, 1, data);

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
            run_analysis(data, size_tot, lp, sum_vx,sum_vy);
#ifdef V_T
      VT_end(analysis_id);
#endif

            t4 = MPI_Wtime();
            t_analy += t4-t3;
            PINF("rank %d: Step %d moments calculated,t_prepare %lf, t_read %lf, t_close %lf, t_analy %lf, time%lf\n", rank, timestep, t1-t0, t2-t1, t3-t2, t4-t3, t4);
        }
    }

    free (data);
    adios_selection_delete (sel);
    //adios_free_varinfo (v);

#ifdef ENABLE_TIMING

    printf("[rank %d]:analysis_time %.3lf \n", rank, t_analy);
    MPI_Barrier(comm);
    double t_end = MPI_Wtime();


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
