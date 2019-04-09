
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

# include <stdlib.h>
# include <stdio.h>



#include "readParam.hpp"
#include "transports.h"
#include "utility.h"
static transport_method_t transport;
#include "nmoments-anal/run_analysis.h"

#ifdef V_T
#include <VT.h>
#endif
//# define n 48            /* matrix is nxn, excluding boundary values     */
//# define nodeedge 24     /* a task works on a nodeedge x nodeedge matrix */
//# define nblock n/nodeedge   /* number of tasks per row of matrix            */
//# define nproc nblock*nblock /* total number of tasks (processors)           */


using namespace decaf;
using namespace std;


  uint64_t nx, ny, nodeedge, nblock, nproc;
    static char var_name[STRING_LENGTH];                                                                                                           
    static char var_size[STRING_LENGTH];                                                                                                           
    static size_t elem_size=sizeof(double);  
	int maxStep = 10;
int main ( int argc, char **argv );
void con (Decaf* decaf, int nsteps); 

void local_con (Decaf* decaf, int nsteps) ; 
/******************************************************************************/


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
        local_con(decaf, nsteps);

    // MPI_Barrier(MPI_COMM_WORLD);

    // cleanup
	printf("before delete workflow rank %d \n", rank);

    delete decaf;
    MPI_Finalize();
}

void local_con (Decaf* decaf, int nsteps)  
{
     /*
     * @input
     * @param NSTOP
     * @param FILESIZE2PRODUCE
     */

    int lp = N_LP;
	uint64_t slice_size = 0;
	double * buffer;
  	int iconf[7];
	int ntasks_prod;
  	double conf[2];
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
	MPI_Comm comm_con = decaf->con_comm_handle();
	int rank = decaf->con_comm()->rank();
	int ntasks_con;
    MPI_Comm_size ( decaf->con_comm_handle(), &ntasks_con);



   /* Get input parameters */
   if(rank==0)
    readParam(iconf, conf);

   /* Broadcast input parameters */
   MPI_Bcast(iconf,7,MPI_INT,0,comm_con);
   MPI_Bcast(conf,2,MPI_DOUBLE,0,comm_con);

   /* Assign input parameters to variables */
   nodeedge    = iconf[0];
   nodeedge    = iconf[1];
   nblock = ntasks_prod;
   int maxStep   = iconf[6];
   ny = nodeedge;	
    char nodename[256];
    int nodename_length;
	double start_t = MPI_Wtime();
    /*
     * define the trace
     */

#ifdef V_T
     VT_classdef( "Analysis", &class_id );
     VT_funcdef("ANL", class_id, &analysis_id);
#endif
    

    vector< pConstructData > in_data;

	int step = 0;
  while (decaf->get(in_data))
    {
        // get the values

		// dan debug
		printf("[nmoments]: decaf get in_data size %d at rank %d \n",in_data.size(), rank);	
        for (size_t i = 0; i < in_data.size(); i++)
        {
            VectorFliedd pos = in_data[i]->getFieldData<VectorFliedd>("pos");
            if (pos)
            {
                // debug
                slice_size = pos.getNbItems();

                if(rank == 0){
                fprintf(stderr, "[nmoments]: consumer processing %d atoms at step %d\n",
                        slice_size,
                        step);
                }


                buffer = &pos.getVector()[0];

                t1 =MPI_Wtime(); 

#ifdef V_T
      VT_begin(analysis_id);
#endif

            run_analysis(buffer, slice_size/2, lp, sum_vx,sum_vy);
#ifdef V_T
      VT_end(analysis_id);
#endif

                //run_analysis(buffer, slice_size, lp, sum_vx,sum_vy);

                t2 =MPI_Wtime(); 
                t_analy += t2-t1;

            }
            else
                fprintf(stderr, "[nmoments]: Error: null pointer in node2\n");
        }
        step+=1;
    }




    //adios_free_varinfo (v);

#ifdef ENABLE_TIMING

    printf("[nmoments]: [rank %d]:analysis_time %.3lf \n", rank, t_analy);
    MPI_Barrier(comm_con);
    double t_end = MPI_Wtime();


        double global_t_prepare=0;
        double global_t_get=0;
        double global_t_close=0;
        double global_t_analy=0;
        MPI_Reduce(&t_prepare, &global_t_prepare, 1, MPI_DOUBLE, MPI_SUM, 0, comm_con);
        MPI_Reduce(&t_get, &global_t_get, 1, MPI_DOUBLE, MPI_SUM, 0, comm_con);
        MPI_Reduce(&t_close, &global_t_close, 1, MPI_DOUBLE, MPI_SUM, 0, comm_con);
        MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm_con);

    if(rank == 0){
      PINF("[nmoments]: stat:Consumer end  at %lf \n", t_end);
      PINF("stat:time for prepare %fs, read %f s; time for close %f s; max time for analy %f s\n",global_t_prepare/ntasks_con, global_t_get/ntasks_con, global_t_close/ntasks_con, global_t_analy);
    }
#endif

                                         
	printf("total analysis time : %f \n", MPI_Wtime() - start_t);
    //MPI_Barrier (comm_con);

    PINF("rank %d: exit\n", rank);

  /*
   * close logger
   */
	decaf->terminate();


}
// test driver for debugging purposes
// normal entry point is run(), called by python
int main(int argc,
         char** argv)
{
    printf("main function launched\n");
    Workflow workflow;
    Workflow::make_wflow_from_json(workflow, "laplace_decaf_noedge.json");

    // run decaf
char * prefix         = getenv("DECAF_PREFIX");
    if (prefix == NULL)
    {
        fprintf(stderr, "ERROR: environment variable DECAF_PREFIX not defined. Please export "
                "DECAF_PREFIX to point to the root of your decaf install directory.\n");
        exit(1);
    }

    if(argc !=3){
        fprintf(stderr, "[laplace]: cmd nsteps infile\n");
    }
    int nsteps = atoi(argv[1]);
    string infile = argv[2];

    con(workflow, nsteps, infile);

           
    return 0;
}


