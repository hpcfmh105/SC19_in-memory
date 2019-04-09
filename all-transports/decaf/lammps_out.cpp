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
void lammps_noedge(Decaf* decaf, int nsteps, string infile)
{

    //int nsteps = NSTEPS;
    int rank;
    int line;

	/*those are all the information io libaray need to know about*/
    MPI_Comm comm;
	int nlocal; //nlines processed by each process
	int size_one = SIZE_ONE; // each line stores 2 doubles
	double *buffer; // buffer address
    double **x;// all the atom values

    double t_start = MPI_Wtime();
#ifdef V_T
      
      //VT_initialize(NULL, NULL);
      printf("[decaf]: trace enabled and initialized\n");
      VT_classdef( "Computation", &class_id );
      VT_funcdef("ADVSTEP", class_id, &advance_step_id);
      VT_funcdef("GETBUF", class_id, &get_buffer_id);
      VT_funcdef("PUT", class_id, &put_buffer_id);
#endif



    LAMMPS* lps = new LAMMPS(0, NULL, decaf->prod_comm_handle());
    lps->input->file(infile.c_str());
    printf("prod lammps_decaf started with input %s\n", infile.c_str() );

    rank = decaf->prod_comm()->rank();


    for (int timestep = 0; timestep < nsteps; timestep++)
    {

#ifdef V_T
      VT_begin(advance_step_id);
#endif
        lps->input->one("run 1 pre no post no");

#ifdef V_T
      VT_end(advance_step_id);
#endif
        //int natoms = static_cast<int>(lps->atom->natoms);
        //lammps_gather_atoms(lps, (char*)"x", 1, 3, x);

        //extract "value"
#ifdef V_T
      VT_begin(get_buffer_id);
#endif
        x = (double **)(lammps_extract_atom(lps,(char *)"x"));
        nlocal = static_cast<int>(lps->atom->nlocal); // get the num of lines this rank have
        if(x == NULL){
            fprintf(stderr, "extract failed\n");
            break;
        }

        buffer = new double[size_one * nlocal];

#if DEBUG
        printf("step %d i have %d lines\n",timestep, nlocal);
#endif
        for(line = 0; line < nlocal; line++){
            buffer[line*size_one] = line;
            buffer[line*size_one+1] = 1;
            buffer[line*size_one+2] = x[line][0];
            buffer[line*size_one+3] = x[line][1];
            buffer[line*size_one+4] = x[line][2];
        }

#ifdef V_T
      VT_end(get_buffer_id);
#endif

#ifdef V_T
      VT_begin(put_buffer_id);
#endif
        /* decaf put */
        if (1)
        {
            pConstructData container;

            // lammps gathered all positions to rank 0
            //if (decaf->prod_comm()->rank() == 0)
            if (rank == 0)
            {
                fprintf(stderr, "lammps producing time step %d with %d lines\n",
                        timestep, nlocal);

            }
                // debug
                //         for (int i = 0; i < 10; i++)         // print first few atoms
                //           fprintf(stderr, "%.3lf %.3lf %.3lf\n",
                // x[3 * i], x[3 * i + 1], x[3 * i + 2]);

            VectorFliedd data(buffer, size_one * nlocal, size_one);

            container->appendData("pos", data,
                                      DECAF_NOFLAG, DECAF_PRIVATE,
                                      DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);
			if (rank == 0)
			{
				printf("after lammps append data into decaf\n");
			}
            /*else*/
            //{
                //vector<double> pos;
                //VectorFliedd data(pos, 3);
                //container->appendData("pos", data,
                                      //DECAF_NOFLAG, DECAF_PRIVATE,
                                      //DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);
            /*}*/

            decaf->put(container);
			if (rank == 0)
			{
				printf("after lammps put data into decaf\n");
			}
        }

#ifdef V_T
      VT_end(put_buffer_id);
#endif

       free(buffer);
    }

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    //
    double t_end = MPI_Wtime();
    printf("[lammps]:total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

#ifdef V_T
    /*VT_finalize();*/
    /*printf("[itac]: trace finalized");*/
#endif

    if(rank == 0){
        fprintf(stderr, "[lammps] now terminating\n");
    }
    decaf->terminate();

    delete lps;
}


void prod(Workflow& workflow,              // workflow
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
        lammps_noedge(decaf, nsteps, infile);

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

    prod(workflow, nsteps, infile);

           
    return 0;
}
