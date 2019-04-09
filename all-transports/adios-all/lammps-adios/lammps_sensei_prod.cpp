#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <utility>
#include <map>


#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#define SIZE_ONE (5)

//#include "utility.h"

//#include "adios_adaptor.h"
//#include "adios_helper.h"
#include "adios.h"
#include "adios_error.h"
#include "adios_adaptor.h"
#include "ds_adaptor.h"

#include "transports.h"
#include "utility.h"
#include "msd-anal/run_msd.h"

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <mpi.h>
#include "library.h"        /* this is a LAMMPS include file */
#include "domain.h"         /* this is a LAMMPS include file */

#include "lammps.h"
#include "modify.h"
#include "fix.h"
#include "fix_external.h"


//#include "senseiConfig.h"   /* SENSEI bridge */
//#ifdef ENABLE_SENSEI
//#include "lammpsBridge.h"
//#endif
//#else
#include "lammps_sensei_analysis.h"
//#endif

//#include <Timer.h>
#include "msd-anal/run_msd.h"

static transport_method_t transport;

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


/** Native staging need to use these*/
static char var_name[STRING_LENGTH];
static char var_size[STRING_LENGTH];
static size_t elem_size=sizeof(double);

/*double xsublo = 0.0;
double ysublo = 0.0;
double zsublo = 0.0;

double xsubhi = 0.0;
double ysubhi = 0.0;
double zsubhi = 0.0;*/

int main(int argc, char * argv[])
{
    /*******************************************************/
	if(argc > 4)
	{
		printf("read argc = %d,lammps_adios_prod nstep inputfile nlocal\n",
                            argc);
        	exit(-1);
	}
       // #ifdef ENABLE_SENSEI
        //std::string sensei_xml;
        //sensei_xml="";
	//sensei_xml = argv[4];
       // #endif


	/*****************************************************/
	int nsteps; //run how many steps

	nsteps = atoi(argv[1]);
	string infile = argv[2];
	uint64_t tmplocal = 0;
	if(argc == 4)
	{
		tmplocal = atoi(argv[3]);
	}

	MPI_Init(&argc, &argv);

	MPI_Comm comm = MPI_COMM_WORLD;
	int rank, nprocs;
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &nprocs);

	char nodename[256];
	int nodename_length;
	MPI_Get_processor_name(nodename, &nodename_length);

	uint64_t * sizes = (uint64_t *) malloc(3 * sizeof(uint64_t));

        if(rank == 0)
	{
		printf("stat: Producer start at %lf \n", MPI_Wtime());
        }

    	/* * lammps*/
    	int step;
    	int line;
    	int nlocal=0; //nlines processed by each process
    	int size_one = SIZE_ONE; // each line stores 2 doubles
    	double *buffer; // buffer address
    	double **x /*= NULL*/;// all the atom values

    	LAMMPS* lps = new LAMMPS(0, NULL, comm);
    	lps->input->file(infile.c_str());
    	printf("prod lammps_decaf started with input %s\n", infile.c_str() );
        
        //timer::Initialize();

        //#ifdef ENABLE_SENSEI
	//lammpsBridge::Initialize(comm, nsteps, nlocal, x, xsublo, xsubhi,
                                 //ysublo, ysubhi, zsublo, zsubhi, sensei_xml);
        //#endif


    	double t_start = MPI_Wtime();
    	nlocal = static_cast<int>(lps->atom->nlocal);
    	buffer = (double *)malloc(size_one * 512000 *sizeof(double));
    	// output dir for adios
        double **msd;
        msd =  init_msd(nsteps, size_one);
    	for (step = 0; step < nsteps; step++)
    	{ 
                double t_simu_start, t_simu_end, /*t_anal_start,*/ t_anal_end;
                t_simu_start = MPI_Wtime();
		if(rank == 0)
		{
			printf("before lps input one\n");
		}
        	lps->input->one("run 1 pre no post no");
        	x = (double **)(lammps_extract_atom(lps,(char *)"x"));
        	nlocal = static_cast<int>(lps->atom->nlocal); // get the num of lines this rank have
        	if(x == NULL)
		{
            		fprintf(stderr, "extract failed\n");
            		break;
       		}	

        	nlocal = 512000;
                int line_buffer = nlocal;
        	if(rank == 0)
		{
            		printf("step %d i have %d lines\n",step, nlocal);
		}
        	for(line = 0; line < nlocal && line < line_buffer; line++)
		{
	            	buffer[line*size_one] = line;
        	    	buffer[line*size_one+1] = 1;
	            	buffer[line*size_one+2] = x[line][0];
	            	buffer[line*size_one+3] = x[line][1];
	            	buffer[line*size_one+4] = x[line][2];
	        }
   
		
        	
            	uint64_t bounds[6] = {0};
            	uint64_t bounds_size[6] = {0};
            	bounds_size[3] = 2;
            	double time_comm;
            	//Dan for small test;
            	if(tmplocal == 0)
			nlocal = 512000;
            	else
			nlocal = tmplocal;
            	sizes[2] = nlocal;
            	sizes[1] = nprocs;
            	sizes[0] = size_one;
            	// xmin
            	bounds[5]=(nlocal) - 1;
            	// ymin
            	bounds[2]=0;
            	bounds[3] = rank;
            	bounds[0] = rank;
            	// xmax
            	bounds[1]= 0;
            	// ymax
            	bounds[4]= (size_one)-1;
//            	put_common_buffer(transport_minor, step ,3, bounds,rank , var_name, (void **)&buffer, elem_size, &time_comm, &comm);
               
               // #ifdef ENABLE_SENSEI
               // lammpsBridge::update_bridge(step, buffer);
               // #else
                t_simu_end = MPI_Wtime();
                printf("---This is the %d step, %d rank, the simulation time is  %.3f.\n", step, rank, t_simu_end-t_simu_start);   
                analyze( nsteps, step, buffer, size_one, sizes[2], msd);
                t_anal_end = MPI_Wtime();
                printf("---This is the %d step, %d rank, the analysis time is  %.3f ---.\n", step, rank, t_anal_end-t_simu_end);
               // #endif
   	 }
         perform_msd_reduce(msd, nsteps, comm);
         free_msd(msd, size_one);
	 free(buffer);
	 free(sizes);
   	 double t_end = MPI_Wtime();
    	 printf("[lammps]:total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

    	 t_end = MPI_Wtime();
    	 MPI_Barrier(comm);
       	if (rank == 0)
	{
		printf("rank 0 prod total time: %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);
	}

    	/* finalize*/
       // #ifdef ENABLE_SENSEI
	//lammpsBridge::Finalize();
       // #endif
	//timer::Finalize();

	MPI_Finalize();
    	printf("rank %d: exit\n", rank);
    	return 0;
} 
