	// lammps includes
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
	static transport_method_t transport;

	#ifdef V_T
	#include <VT.h>
	int class_id;
	int advance_step_id, get_buffer_id, put_buffer_id;
	//int analysis_id;
	#endif




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

	/*
	 * clog identifier
	 */

	/*
	 * Native staging need to use these
	 */
	static char var_name[STRING_LENGTH];
	static char var_size[STRING_LENGTH];
	static size_t elem_size=sizeof(double);

	int main(int argc, char * argv[]){

		/*
		 * @input
		 * @param NSTOP
		 * @param FILESIZE2PRODUCE
		 * @nlocal
		 */
		if(argc > 4){
			printf("read argc = %d,lammps_adios_prod nstep inputfile nlocal\n",
					argc);
			exit(-1);
		}
		int nsteps; //run how many steps

		nsteps = atoi(argv[1]);
		string infile = argv[2];
		uint64_t tmplocal = 0;
		if(argc == 4){
			tmplocal = atoi(argv[3]);
		}

		MPI_Init(&argc, &argv);

		MPI_Comm comm = MPI_COMM_WORLD;
		int         rank, nprocs;
		MPI_Comm_rank (comm, &rank);
		MPI_Comm_size (comm, &nprocs);

		char nodename[256];
		int nodename_length;
		MPI_Get_processor_name(nodename, &nodename_length);

		uint64_t * sizes = (uint64_t *) malloc(3 * sizeof(uint64_t));


		/*
		 * get transport method from env variable
		 */
		transport = get_current_transport();

		uint8_t transport_major = get_major(transport);
		uint8_t transport_minor = get_minor(transport);
		printf("%s:I am rank %d of %d, tranport code %x-%x\n",
				nodename, rank, nprocs,
				get_major(transport), get_minor(transport) );
		if(rank == 0){
			printf("stat: Producer start at %lf \n", MPI_Wtime());
		}


		if(transport_major == ADIOS_DISK || transport_major == ADIOS_STAGING){

			char xmlfile[256], trans_method[256];

			if(transport_major == ADIOS_DISK){
				strcpy(trans_method, "mpiio");
			}
			else{
				if(transport_minor == DSPACES)
					strcpy(trans_method, "dataspaces");
				else if(transport_minor == DIMES) 
					strcpy(trans_method, "dimes");

				else if(transport_minor == FLEXPATH)
					strcpy(trans_method, "flexpath");
			}

			//sprintf(xmlfile,"xmls/dbroker_%s.xml", trans_method);
			sprintf(xmlfile,"xmls/lammps_%s.xml", trans_method);
			if(rank == 0)
				printf("[r%d] try to init with %s\n", rank, xmlfile);
			const char* pnid = getenv("P2TNID");
			printf("[r%d] show P2TNID: %s\n", rank,  pnid);
			if(adios_init (xmlfile, comm) != 0 ){
				fprintf( stderr, "[r%d] ERROR: adios init err with %s\n", rank, trans_method);
				fprintf(stderr, "[r%d] ERR: %s\n", rank, adios_get_last_errmsg());
				return -1;
			}
			else{
				//if(rank ==0)
				printf( "rank %d : adios init complete with %s\n", rank, trans_method);
			}
			MPI_Barrier(comm);
		} //use ADIOS_DISK or ADIOS_STAGING

		else  if(transport_major == NATIVE_STAGING){

	        if(S_OK != ds_adaptor_init_client(nprocs, 1, &comm, NULL)){
            	TRACE();
				printf( "dataspaces init error");
            	MPI_Abort(comm, -1);
        	}	

        	sprintf(var_name, "atom_0");
			sprintf(var_size, "atom_size");

			printf("trying init dspaces for %d process\n", nprocs);
			printf("dspaces init successfuly \n");
		}


		/* 
		 * lammps
		 */
		int step;
		int line;
		int nlocal; //nlines processed by each process
		int size_one = SIZE_ONE; // each line stores 2 doubles
		double *buffer; // buffer address
		double **x;// all the atom values

	#ifdef V_T

		//VT_initialize(NULL, NULL);
		printf("[decaf]: trace enabled and initialized\n");
		VT_classdef( "Computation", &class_id );
		VT_funcdef("ADVSTEP", class_id, &advance_step_id);
		VT_funcdef("GETBUF", class_id, &get_buffer_id);
		VT_funcdef("PUT", class_id, &put_buffer_id);
	#endif




		LAMMPS* lps = new LAMMPS(0, NULL, comm);
		lps->input->file(infile.c_str());
		printf("prod lammps_decaf started with input %s\n", infile.c_str() );

		double t_start = MPI_Wtime();
		nlocal = static_cast<int>(lps->atom->nlocal); 	
		buffer = (double *)malloc(size_one * 512000 *sizeof(double));
		// output dir for adios
		char *filepath = getenv("BP_DIR");
		//FILE * fdm = fopen ("prodsize","w"); 
		for (step = 0; step < nsteps; step++)
		{

	#ifdef V_T
			VT_begin(advance_step_id);
	#endif
			if(rank == 0)
				printf("before lps input one\n");
			lps->input->one("run 1 pre no post no");
	#ifdef V_T
			VT_end(advance_step_id);
	#endif

	#ifdef V_T
			VT_begin(get_buffer_id);
	#endif
			x = (double **)(lammps_extract_atom(lps,(char *)"x"));

			nlocal = static_cast<int>(lps->atom->nlocal); // get the num of lines this rank have
			if(x == NULL){
				fprintf(stderr, "extract failed\n");
				break;
			}
			int natoms = static_cast<int>(lps->atom->natoms);
			int navg = (int) natoms/nprocs; //avg atoms per proc
			// TODO: precise!
	#ifdef PRECISE
			int line_buffer = nlocal; // how many lines for buffer
	#else
			//Dan modify int line_buffer = navg;
			int line_buffer = nlocal;
			if(rank == 0)
				printf("[warning]: use estimate lines line buffer: %d\n", line_buffer);
	#endif

			// buffer = (double *)malloc(size_one * line_buffer*sizeof(double));

			nlocal = 512000;
			if(rank == 0)
				printf("step %d i have %d lines\n",step, nlocal);
			for(line = 0; line < nlocal && line < line_buffer; line++){
				buffer[line*size_one] = line;
				buffer[line*size_one+1] = 1;
				buffer[line*size_one+2] = x[line][0];
				buffer[line*size_one+3] = x[line][1];
				buffer[line*size_one+4] = x[line][2];
			}

	#ifdef V_T
			VT_end(get_buffer_id);
	#endif

			/*
			 * insert intto adios
			 */



	#ifdef V_T
			VT_begin(put_buffer_id);
	#endif
			char        filename [256];
			int         offset, size_y;
			int         NX = size_one; // this is the longest dimension
			int        NY = 1;
			int        NZ = nlocal;
			int64_t     adios_handle;

			offset = rank*NY;
			size_y = nprocs*NY;

			if (rank == 0)
				printf("before save data to adios\n");
			if(transport_major == ADIOS_DISK || transport_major == ADIOS_STAGING){
				if(transport_major == ADIOS_STAGING)
					sprintf(filename, "%s/atom.bp", filepath);
				else
					sprintf(filename, "%s/atom_%d.bp", filepath, step);

				if (rank == 0)
					printf("before adios_open\n");
				if(err_no_error != adios_open (&adios_handle, "temperature", filename, "w", comm)){
					PERR("cannot open");
					TRACE();
					MPI_Abort(comm, -1);

				}
				if (rank == 0)
					printf("before adios write\n");
				adios_write (adios_handle, "/scalar/dim/NX", &NX);
				adios_write (adios_handle, "/scalar/dim/NY", &NY);
				adios_write (adios_handle, "/scalar/dim/NZ", &NZ);
				//adios_write (adios_handle, "test_scalar", &test_scalar);
				adios_write (adios_handle, "size", &nprocs);
				adios_write (adios_handle, "rank", &rank);
				adios_write (adios_handle, "offset", &offset);
				adios_write (adios_handle, "size_y", &size_y);
				adios_write (adios_handle, "var_2d_array", buffer);

				if(rank == 0)
					printf("before adios close\n");
				if(err_no_error != adios_close (adios_handle)){
					PERR("cannot close");
					TRACE();
					MPI_Abort(comm, -1);
				}
				if(rank == 0)
					printf("after adios close\n");
				if(transport_major == ADIOS_DISK){
					
					/**** use index file to keep track of current step *****/
					char step_index_file[256];
					sprintf(step_index_file, "%s/stamp.file", filepath);

					if(S_OK != adios_adaptor_update_avail_version(comm, step_index_file, step, nsteps)){
						PERR("index file not found");
						TRACE();
						MPI_Abort(comm, -1);
					}
				}
			}
			else if(transport_major == NATIVE_STAGING){
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
				/*if(step == 0){
					fprintf (fdm, "%ld\n", sizes[0]);
					fprintf (fdm, "%ld\n", sizes[1]);
					fprintf (fdm, "%ld\n", sizes[2]);
					fflush(fdm);
				}	*/	
//			status_t put_common_buffer(uint8_t transport_minor, int timestep,int ndim, uint64_t bounds[6], int rank,char * var_name, void  **p_buffer,size_t elem_size, double *p_time_used)	
				//if(rank == 0 && step == 0)
			//	put_common_buffer(transport_minor, step,3, bounds_size,rank, var_size, (void **) &sizes, sizeof(uint64_t), &time_comm, &comm);
				put_common_buffer(transport_minor, step ,3, bounds,rank , var_name, (void **)&buffer, elem_size, &time_comm, &comm);
				//t_put+=time_comm;
			}


			else{
				PERR("transport %u:%u is not not supported", transport_major, transport_minor);
				TRACE();
				MPI_Abort(comm, -1);
			}



			/*// dimension sizexnlocalxsize_one*/
			//adios_open (&adios_handle, "atom", filename, "w", comm);

			//adios_write (adios_handle, "size", &nprocs);
			//adios_write (adios_handle, "rank", &rank);
			//adios_write (adios_handle, "nlocal", &const_nlocal);
			//adios_write (adios_handle, "size_one", &size_one);
			//adios_write (adios_handle, "array", buffer);

			//// end of gwrite_atom

			//adios_close (adios_handle);


			//insert_into_Adios(transport, var_name, step,nsteps, line_buffer, size_one, buffer,"w" , &comm);

	#ifdef V_T
			VT_end(put_buffer_id);
	#endif

		}
		//fclose(fdm);
		free(buffer);
		free(sizes);
		double t_end = MPI_Wtime();
		printf("[lammps]:total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

		t_end = MPI_Wtime();
		MPI_Barrier(comm);
		if (rank == 0)
			printf("rank 0 prod total time: %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

		/*
		 * finalize
		 */
		//Dan modify , not finalize dataspace due to dataspaces stucked.
//        MPI_Finalize();
  //      printf("rank %d: exit\n", rank);
    //    return 0;


		if(transport_major == ADIOS_DISK || transport_major == ADIOS_STAGING){
			adios_finalize (rank);
			fprintf(stderr, "rank %d: adios finalize complete\n", rank); 
		}

		else if(transport_major == NATIVE_STAGING){
			
			/* dimes needs to flush last step */
			if(transport_minor == DIMES){
				ds_adaptor_flush_dimes(var_name, comm, nsteps);
			}

			double t_end = MPI_Wtime();
			  if(rank == 0){
				  PINF("stat:Simulation stop at %lf \n", t_end);
			  }

			dspaces_put_sync();
			dspaces_finalize();
		}

		MPI_Finalize();
		printf("rank %d: exit\n", rank);
		return 0;
	}

