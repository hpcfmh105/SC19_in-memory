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
static transport_method_t transport;

#define SIZE_ONE (5)

#ifdef V_T
#include <VT.h>
int class_id;
int analysis_id;
#endif

static char var_name[STRING_LENGTH];
static char var_size[STRING_LENGTH];
static size_t elem_size=sizeof(double);

//#define DEBUG_Feng
int main (int argc, char ** argv){

     /*
     * @input
     * @param NSTOP
	 * @nprod/ncon
	 * @nlocal
     */
    if(argc > 4){
        printf("lammps_adios_con nsteps nprod/ncon nlocal\n");
        exit(-1);
    }

    int nsteps = atoi(argv[1]);
    double **msd;
 	
	uint64_t slice_size = 1;
	if(argc > 2){
		slice_size = atoi(argv[2]);

	}
	uint64_t tmplocal = 0;
	if(argc == 4)
		tmplocal = atoi(argv[3]);
	uint64_t* sizes = (uint64_t *) malloc(3 * sizeof(uint64_t));
    /* init adios */
    int         rank, nprocs;
    MPI_Comm    comm = MPI_COMM_WORLD;


    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_DIMES;
    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    double * data = NULL;
    uint64_t start[3], count[3];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nprocs);

    char nodename[256];
    int nodename_length;
    MPI_Get_processor_name(nodename, &nodename_length );


#ifdef V_T
     VT_classdef( "Analysis", &class_id );
     VT_funcdef("ANL", class_id, &analysis_id);
#endif



    /******************** configuration stop ***********/

#ifdef ENABLE_TIMING
    double t1, t2, t3, t4;
    double t_read_1, t_read_2, t_analy;
    t_read_1 = 0;
    t_read_2 = 0;
    t_analy = 0;
#endif
    int step =0;

    /*
     * get transport method
     */
    transport = get_current_transport();
    uint8_t transport_major = get_major(transport);
    uint8_t transport_minor = get_minor(transport);
    printf("%s:I am rank %d of %d, tranport code %x-%x\n",
            nodename, rank, nprocs,
            get_major(transport), get_minor(transport) );

    if(rank == 0){
      printf("stat: Consumer start at %lf \n", MPI_Wtime());
    }
    assert(transport_major == NATIVE_STAGING);

    char filename[256];

    if(S_OK != ds_adaptor_init_client(nprocs, 2, &comm, NULL)){
        TRACE();
        MPI_Abort(comm, -1);
    }

    double t_start = MPI_Wtime();
    /*
    * set bounds and dspaces variables
    */
    sprintf(var_name, "atom_0");
	sprintf(var_size, "atom_size");
    // how many lines in global
    //long * y_size = (long *) malloc(sizeof (long));
    //uint64_t gdims[2] = {2, global_size};
    //dspaces_define_gdim(var_name, 2,gdims);
    
	uint64_t bounds_size[6] = {0};
	bounds_size[3] =2; 
	double time_comm;

	int size_one = SIZE_ONE;

    msd =  init_msd(nsteps, size_one);
    int ii = 0;
    double      *t = NULL;
	sizes[0] = nprocs;
	sizes[1] = size_one;
	if(tmplocal == 0)
		sizes[2] = 512000;
	else
		sizes[2] = tmplocal;
   	data = (double *)malloc (slice_size * sizes[2] * SIZE_ONE* sizeof (double));
    while(ii < nsteps){       

        size_t nelem;
   		         

		//if(S_OK != get_common_buffer(transport_minor, ii,3, bounds_size,rank, var_size, (void **)& sizes, sizeof(uint64_t), &time_comm, &comm)){
		//	  TRACE();
		//	  MPI_Abort(comm, -1);
		//}
		// Dan if prod procs / con procs = 2.

		//uint64_t slice_size = (uint64_t) (sizes[0])/(nprocs);
		
		//printf("total sizes: %ld, %ld, %ld \n", sizes[0], sizes[1], sizes[2]);
		uint64_t bounds[6] = {0};
					// xmin
		bounds[2] = 0;
		bounds[1]=slice_size * rank;
		bounds[0]= 0;
					// ymin
		bounds[5] = sizes[2] - 1;
					// xmax
		bounds[4]=slice_size * (rank + 1) -1;
					// ymax
		bounds[3]= EVAL(size_one)-1;


         if(S_OK != get_common_buffer(transport_minor, ii,3, bounds,rank, var_name, (void **)&data, elem_size, &time_comm, &comm)){
          TRACE();
			printf("ERROR, after get_common_buffer in main rank %d \n", rank);
          MPI_Abort(comm, -1);
        }
		/*
		bounds[1] = slice_size * (rank + 1) - 1;
    	bounds[4]=slice_size * (rank + 1) - 1;
		double * data2 = data + sizes[0] * SIZE_ONE* sizeof (double); 
         if(S_OK != get_common_buffer(transport_minor, ii,3, bounds,rank, var_name, (void **)&data2, elem_size, &time_comm)){
          TRACE();
          MPI_Abort(comm, -1);
        }
		
		*/
		printf("after get common buffer step: %d \n", ii);
        //sleep(20);
    
#if 0
        printf("Rank=%d: test_scalar: %d step: %d, t[0,5+x] = [", rank, test_scalar, ii);
        for(j=0; j<nelem; j++) {
            printf(", %6.2f", t[j]);
        }
        printf("]\n\n");
#endif
        t1 =MPI_Wtime(); 

#ifdef V_T
      VT_begin(analysis_id);
#endif
		printf("before calc_msd %ld, %ld, %d \n",  sizes[2], size_one, ii);
        calc_msd(msd, data, slice_size * sizes[2], size_one, ii);
#ifdef V_T
      VT_end(analysis_id);
#endif
        t2 =MPI_Wtime(); 

        t_analy += t2-t1;

        ii++;
        //MPI_Barrier (comm);
        //sleep(1);
    }
    //
  /*
     * reduce
     */
    perform_msd_reduce(msd, nsteps, comm);
    free_msd(msd, size_one);

    double t_end = MPI_Wtime();
    double global_t_analy;


    /*
     * output
     */
    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    
    printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
        printf("[msd]: max t_analysis %.3fs, now existing\n", global_t_analy);
    }

    if(t){
        free(t);
    }
    MPI_Finalize ();
    return 0;


//#ifdef RAW_DSPACES
    dspaces_finalize();
//#endif


    /* 
     * get the dimension info
     */
#if 0
    uint64_t slice_size = v->dims[1]/nprocs;
    start[1] = slice_size * rank;
    if (rank == nprocs-1) /* last rank may read more lines */
        slice_size = slice_size + v->dims[1]%nprocs;
    count[1] = slice_size;


    start[0] = 0;
    count[0] = v->dims[0];

    start[2] = 0;
    count[2] = v->dims[2];

    printf("rank %d: start: (%ld, %ld, %ld), count:( %ld, %ld, %ld)\n", rank, start[0], start[1], start[2], count[0], count[1], count[2]);

    int size_one = v->dims[0];
    int nlines = slice_size*(v->dims[2]);
       
    msd =  init_msd(nsteps, size_one);

    data = (double *)malloc (nlines*size_one* sizeof (double));
    memset(data, 0, nlines*size_one*sizeof(double));
    if (data == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    sel = adios_selection_boundingbox (v->ndim, start, count);


    printf("rank %d: adios init complete\n", rank);

    double t_start = MPI_Wtime();
    while(adios_errno != err_end_of_stream){
        if(rank == 0){
            printf("rank %d: Step %d start\n", rank, step);
        }

        // read
        adios_schedule_read (f, sel, "array", 0, 1, data);
        adios_perform_reads (f, 1);
        adios_release_step(f);
        adios_advance_step(f, 0, -1);

        t1 =MPI_Wtime(); 

#ifdef V_T
      VT_begin(analysis_id);
#endif
        calc_msd(msd, data, nlines, size_one, step);
#ifdef V_T
      VT_end(analysis_id);
#endif
        t2 =MPI_Wtime(); 

        t_analy += t2-t1;

        step ++;
    }

    /*
     * reduce
     */
    perform_msd_reduce(msd, nsteps, comm);
    free_msd(msd, size_one);

    double t_end = MPI_Wtime();
    double global_t_analy;


    /*
     * output
     */
    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    
    printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
		printf("[msd]: max t_analysis %.3fs, now existing\n", global_t_analy);
	}


    /*
     * finalize adios
     */

    free (data);
    adios_selection_delete (sel);
    adios_free_varinfo (v);
    adios_read_close (f);
    MPI_Barrier (comm);
#endif
	free(sizes);
    MPI_Finalize ();
    return 0;
}
