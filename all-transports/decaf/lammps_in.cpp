//---------------------------------------------------------------------------
//
// lammps example
//
// 4-node workflow
//
//          print (1 proc)
//        /
//    lammps (4 procs)
//        \
//          print2 (1 proc) - print (1 proc)
//
//  entire workflow takes 10 procs (1 dataflow proc between each producer consumer pair)
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
//--------------------------------------------------------------------------
#include <decaf/decaf.hpp>
#include <bredala/data_model/pconstructtype.h>
#include <bredala/data_model/vectorfield.hpp>
#include <bredala/data_model/boost_macros.h>

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <utility>
#include <map>

// lammps includes
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "msd-anal/run_msd.h"
#define SIZE_ONE (5)
//#define NSTEPS (100)

#ifdef V_T
#include <VT.h>
int class_id, class_id2;
int advance_step_id, get_buffer_id, put_buffer_id;
int analysis_id;
#endif



using namespace decaf;
using namespace LAMMPS_NS;
using namespace std;

struct lammps_args_t                         // custom args for running lammps
{
    LAMMPS* lammps;
    string infile;
};

struct pos_args_t                            // custom args for atom positions
{
    int natoms;                              // number of atoms
    double* pos;                             // atom positions
};

// runs lammps and puts the atom positions to the dataflow at the consumer intervals
//void lammps(Decaf* decaf, int nsteps, int analysis_interval, string infile)

// gets the atom positions and prints them
void msd_noedge(Decaf* decaf, int nsteps)
{
    double global_t_analy = 0, t_analy = 0;
    double t1, t2;
    int slice_size = 0;
    double *buffer;
    MPI_Comm comm;
    int rank;
    int step;

#ifdef V_T
     VT_classdef( "Analysis", &class_id2 );
     VT_funcdef("ANL", class_id2, &analysis_id);
#endif

    comm = decaf->con_comm_handle();
    rank = decaf->con_comm()->rank();

    vector< pConstructData > in_data;

    step = 0;


    if(rank == 0){
        printf("[msd]: application starts\n");
    }


    /* msd required*/
    double **msd;
    //int nsteps = NSTEPS;
    //int timestep = 0;
    int size_one = SIZE_ONE;
    msd =  init_msd(nsteps, size_one);

    //MPI_Barrier(comm);
    double t_start = MPI_Wtime();
      if(rank == 0){
        printf("[msd]: after init msd starts rank %d\n", rank);
    }


	  while (decaf->get(in_data))
    {
        // get the values

		// dan debug
		printf("decaf get data size %d at rank %d \n, in_data.size, rank");	
        for (size_t i = 0; i < in_data.size(); i++)
        {
            VectorFliedd pos = in_data[i]->getFieldData<VectorFliedd>("pos");
            if (pos)
            {
                // debug
                slice_size = pos.getNbItems();

                if(rank == 0){
                fprintf(stderr, "[msd]: consumer processing %d atoms at step %d\n",
                        slice_size,
                        step);
                }


                buffer = &pos.getVector()[0];

                t1 =MPI_Wtime(); 

#ifdef V_T
      VT_begin(analysis_id);
#endif
                calc_msd(msd, buffer, slice_size, size_one, step);

#ifdef V_T
      VT_end(analysis_id);
#endif

                //run_analysis(buffer, slice_size, lp, sum_vx,sum_vy);

                t2 =MPI_Wtime(); 
                t_analy += t2-t1;

            }
            else
                fprintf(stderr, "[msd]: Error: null pointer in node2\n");
        }
        step+=1;
    }

    perform_msd_reduce(msd, nsteps, comm);
    free_msd(msd, size_one);

   printf("[msd]: t_analy %lf\n", t_analy);
    //MPI_Barrier(comm);
    double t_end = MPI_Wtime();

    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    
    printf("[msd]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
		printf("[msd]: max t_analysis %.3fs, now existing\n", global_t_analy);
	}

    decaf->terminate();
}


void con(Workflow& workflow,              // workflow
        int nsteps,
         string infile)                      // lammps input config file*/
{
    MPI_Init(NULL, NULL);
	int rank;
	 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);
	printf("enter run workflow rank %d \n", rank);
    // run workflow node tasks
    // decaf simply tells the user whether this rank belongs to a workflow node
    // how the tasks are called is entirely up to the user
    // e.g., if they overlap in rank, it is up to the user to call them in an order that makes
    // sense (threaded, alternting, etc.)
    // also, the user can define any function signature she wants
        msd_noedge(decaf, nsteps);

    // MPI_Barrier(MPI_COMM_WORLD);

    // cleanup
	printf("before delete workflow rank %d \n", rank);

    delete decaf;
    MPI_Finalize();
}

// test driver for debugging purposes
// normal entry point is run(), called by python
int main(int argc,
         char** argv)
{
    printf("main function launched\n");
    Workflow workflow;
    Workflow::make_wflow_from_json(workflow, "lammps_decaf_noedge.json");

    // run decaf
char * prefix         = getenv("DECAF_PREFIX");
    if (prefix == NULL)
    {
        fprintf(stderr, "ERROR: environment variable DECAF_PREFIX not defined. Please export "
                "DECAF_PREFIX to point to the root of your decaf install directory.\n");
        exit(1);
    }

    if(argc !=3){
        fprintf(stderr, "[lammps]: cmd nsteps infile\n");
    }
    int nsteps = atoi(argv[1]);
    string infile = argv[2];

    con(workflow, nsteps, infile);

           
    return 0;
}
