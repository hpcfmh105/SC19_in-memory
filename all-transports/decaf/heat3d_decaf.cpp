
#include <decaf/decaf.hpp>
#include <bredala/data_model/pconstructtype.h>
#include <bredala/data_model/vectorfield.hpp>
#include <bredala/data_model/boost_macros.h>

#include <assert.h>
#include <string.h>
#include <utility>
#include <map>
#include "nmoments-anal/run_analysis.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "transports.h"
#include "utility.h"
static transport_method_t transport;

#ifdef V_T
#include <VT.h>
#endif


#define min(a,b) a <= b ? a : b

extern "C" void readParam(int* iconf, double* conf);

extern "C" void processToMap(int me, int* xs, int *xe, int *ys, int *ye, int *zs, int *ze,
                  int xcell, int ycell, int zcell, int x_domains, int y_domains, int z_domains);

extern "C" void initValues(int nb, double*** tab, int a, int b, int c, double temp1, double temp2);

extern "C" void updateBound(double*** x, int sizex, int sizey, int sizez, int* NeighBor, MPI_Comm comm,
                 MPI_Datatype type1, MPI_Datatype type2, MPI_Datatype type3, int current,
                 int* x1, int* y1, int* z1, int* x2, int* y2, int* z2);

extern "C" void computeNext(double*** x0, double*** x, int sizex, int sizey, int sizez, double dt,
                 double hx, double hy, double hz, double* res, int me, int* x1, int* y1, int* z1,
                 int* x2, int* y2, int* z2, double k);

using namespace decaf;
using namespace std;

void prod(Decaf* decaf, int nsteps)
{

   /* Sizes for discretization */
   int size_x, size_y, size_z, me, x_domains, y_domains, z_domains;
   int size_x_glo, size_y_glo, size_z_glo;

   /* Arrays */
   double ***x;
   double ***x0;
   double *x_all;
   double *x0_all;
   double *x_alloc;
   double *x0_alloc;
   double *xfinal;
   double *xtemp;

   /* For reading parameters */
   int iconf[7];
   double conf[2];

   /* Spacing and time steps */
   double dt, dt1, dt2, hx, hy, hz, min1, min2;

   /* Current error and limit convergence */
   double resLoc, result, epsilon;

   /* Output file descriptor */
   FILE* file;

   /* Convergence pseudo-boolean */
   int convergence = 0;

   /* Index variables */
   int i, j, k, l, p, v, m;

   /* Time and step variables */
   double t;
   int step;

   /* Max step */
   int maxStep;

   /* Variables for clock */
   double time_init, time_final;
   double elapsed_time;

   /* Number of initial borders layers to */
   /* avoid artifacts problem on corners  */
   int nb_layers = 1;

   /* temp1_init: temperature init on edges */
   double temp1_init = 10.0;

   /* temp2_init: temperature init inside */
   double temp2_init = -10.0;

   /* Diffusivity coefficient */
   double k0 = 1;

   /* MPI variables */
   int sizes[3], subsizes1[3], subsizes2[3], subsizes3[3], starts[3];
   int nproc, ndims;
   MPI_Comm comm, comm3d;
   int dims[3];
   int periods[3];
   int reorganisation = 0;
   MPI_Datatype matrix_type_oxz, matrix_type_oxy, matrix_type_oyz;
   int S=0, E=1, N=2, W=3, Zd=4, Zu=5;
   int NeighBor[6];
   int xcell, ycell, zcell, size_tot_x, size_tot_y, size_tot_z;
   int *xs, *ys, *zs, *xe, *ye, *ze;

   /* MPI initialization */





   /* Get input parameters */
   if(me==0)
    readParam(iconf, conf);

    me = decaf->prod_comm()->rank();
	comm = decaf->prod_comm_handle();
   /* Broadcast input parameters */
   MPI_Bcast(iconf,7,MPI_INT,0,comm);
   MPI_Bcast(conf,2,MPI_DOUBLE,0,comm);

   /* Assign input parameters to variables */
   size_x    = iconf[0];
   size_y    = iconf[1];
   size_z    = iconf[2];
   x_domains = iconf[3];
   y_domains = iconf[4];
   z_domains = iconf[5];
   maxStep   = iconf[6];
   dt1       = conf[0];
   epsilon   = conf[1];

   /* Warning message if dimensions and number of processes don't match */
   if((me==0) && (nproc!=(x_domains*y_domains*z_domains)))
    printf("Number of processes not equal to Number of subdomains\n");

   /* Various other variables */
   size_x_glo = size_x+2;
   size_y_glo = size_y+2;
   size_z_glo = size_z+2;
   hx = 1.0/(double)(size_x_glo);
   hy = 1.0/(double)(size_y_glo);
   hz = 1.0/(double)(size_z_glo);
   min1 = min(hx,hy);
   min2 = min(min1,hz);
   dt2  = 0.125*min2*min2/k0;
   size_tot_x = size_x+2*x_domains+2;
   size_tot_y = size_y+2*y_domains+2;
   size_tot_z = size_z+2*z_domains+2;

   /* Take a right time step for convergence */
   if(dt1>=dt2)
   {
    if(me==0)
    {
     printf("\n");
     printf("  Time step too large in 'param' file - Taking convergence criterion\n");
    }
    dt = dt2;
   }
   else dt = dt1;

   /* Allocate 3D Contiguous arrays */
   xfinal = (double*)malloc(size_x*size_y*size_z*sizeof(*xfinal));
   x_all = (double *)malloc(size_tot_x*size_tot_y*size_tot_z*sizeof(*x_all)); 
   x0_all = (double *) malloc(size_tot_x*size_tot_y*size_tot_z*sizeof(*x0_all));
   x_alloc = x_all;
   x0_alloc = x0_all;

   /* Allocate size_y rows */
   x = (double ***)malloc(size_tot_x*sizeof(*x));
   x0 = (double ***)malloc(size_tot_x*sizeof(*x0));

   /* Loop on rows */
   for (i=0;i<size_tot_x;i++)
   {
    /* Allocate size_tot_y columns for each row */
    x[i] = (double **)malloc(size_tot_y*sizeof(**x));
    x0[i] = (double **)malloc(size_tot_y*sizeof(**x0));
    /* Loop on columns */
    for(j=0;j<size_tot_y;j++)
    {
     /* Increment size_z block on x0[i][j] address */
     x[i][j] = x_alloc;
     x0[i][j] = x0_alloc;
     x_alloc += size_tot_z;
     x0_alloc += size_tot_z;
    }
   }

   /* Allocate coordinates of processes */
   xs = (int *)malloc(nproc*sizeof(int));
   xe = (int *)malloc(nproc*sizeof(int));
   ys = (int *)malloc(nproc*sizeof(int));
   ye = (int *)malloc(nproc*sizeof(int));
   zs = (int *)malloc(nproc*sizeof(int));
   ze = (int *)malloc(nproc*sizeof(int));

   /* Create 3D cartesian grid */
   periods[0] = 0;
   periods[1] = 0;
   periods[2] = 0;

   ndims = 3;
   dims[0] = x_domains;
   dims[1] = y_domains;
   dims[2] = z_domains;

   MPI_Cart_create(comm, ndims, dims, periods, reorganisation, &comm3d);

   /* Identify neighbors */
   NeighBor[0] = MPI_PROC_NULL;
   NeighBor[1] = MPI_PROC_NULL;
   NeighBor[2] = MPI_PROC_NULL;
   NeighBor[3] = MPI_PROC_NULL;
   NeighBor[4] = MPI_PROC_NULL;
   NeighBor[5] = MPI_PROC_NULL;

   /* Left/West and right/Est neigbors */
   MPI_Cart_shift(comm3d, 0, 1, &NeighBor[W], &NeighBor[E]);

   /* Bottom/South and Upper/North neigbors */
   MPI_Cart_shift(comm3d, 1, 1, &NeighBor[S], &NeighBor[N]);

   /* Zdown/South and Zup/North neigbors */
   MPI_Cart_shift(comm3d, 2, 1, &NeighBor[Zd], &NeighBor[Zu]);

   /* Size of each cell */
   xcell = (size_x/x_domains);
   ycell = (size_y/y_domains);
   zcell = (size_z/z_domains);

   /* Allocate subdomain */
   xtemp = (double *)malloc(xcell*ycell*zcell*sizeof(*xtemp));

   /* Compute xs, xe, ys, ye, zs, ze for each cell on the grid */
   processToMap(me, xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains);

   /* Create matrix data types to communicate */
   sizes[0] = size_tot_x;
   sizes[1] = size_tot_y;
   sizes[2] = size_tot_z;

   starts[0] = 0;
   starts[1] = 0;
   starts[2] = 0;

   /* Create matrix data type to communicate on vertical Oxz plane */
   subsizes1[0] = xcell;
   subsizes1[1] = 1;
   subsizes1[2] = zcell;

   MPI_Type_create_subarray(3, sizes, subsizes1, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxz);
   MPI_Type_commit(&matrix_type_oxz);

   /* Create matrix data type to communicate on vertical Oyz plane */
   subsizes2[0] = 1;
   subsizes2[1] = ycell;
   subsizes2[2] = zcell;

   MPI_Type_create_subarray(3, sizes, subsizes2, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oyz);
   MPI_Type_commit(&matrix_type_oyz);

   /* Create matrix data type to communicate on vertical Oxy plane */
   subsizes3[0] = xcell;
   subsizes3[1] = ycell;
   subsizes3[2] = 1;

   MPI_Type_create_subarray(3, sizes, subsizes3, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxy);
   MPI_Type_commit(&matrix_type_oxy);

   /* Initialize values */
   initValues(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init);

   /* Update the boundaries */
   updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d,
               matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

   int size_x_me = xe[me] - xs[me] + 1;
   int size_y_me = ye[me] - ys[me] + 1;
   int size_z_me = ze[me] - zs[me] + 1;
    double* buffer = (double *) malloc(size_x_me * size_y_me * size_z_me * sizeof(double));
    int xx = 0;
   for (int i = xs[me]; i <= xe[me]; i++){
        for (int j = ys[me]; j <= ye[me]; j++){
            for(int k = zs[me]; k <= ze[me]; k++){
                buffer[xx] = x[i][j][j];
                xx++;
            }
        }
    }

   /* Initialize step and time */
   step = 0;
   t = 0.0;

   /* Starting time */
   time_init = MPI_Wtime();

   /* Main loop */
   while(!convergence)
   {  
      /* Increment step and time */
      step = step + 1;
      t = t + dt ;

      /* Perform one step of the explicit scheme */
      computeNext(x0, x, size_tot_x, size_tot_y, size_tot_z, dt, hx, hy, hz, &resLoc, me, xs, ys, zs, xe, ye, ze, k0);

      /* Update the partial solution along the interface */
      updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d, 
                  matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

      /* Sum reduction to get error */
      MPI_Allreduce(&resLoc, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

      /* Current error */
      result = sqrt(result);

        /* decaf put */
        if (1)
        {
            pConstructData container;

            // lammps gathered all positions to rank 0
            //if (decaf->prod_comm()->rank() == 0)
            if (me == 0)
            {
                fprintf(stderr, "heat3d producing time step %d with %d bytes\n",
                        step, size_x_me * size_y_me * size_z_me * 4);

            }
                // debug
                //         for (int i = 0; i < 10; i++)         // print first few atoms
                //           fprintf(stderr, "%.3lf %.3lf %.3lf\n",
                // x[3 * i], x[3 * i + 1], x[3 * i + 2]);
           xx = 0;
           for (int i = xs[me]; i <= xe[me]; i++){
                for (int j = ys[me]; j <= ye[me]; j++){
                    for(int k = zs[me]; k <= ze[me]; k++){
                        buffer[xx] = x[i][j][j];
                        xx++;
                    }
                }
            }

            VectorFliedd data(buffer, size_x_me * size_y_me * size_z_me, 1);

            container->appendData("pos", data,
                                      DECAF_NOFLAG, DECAF_PRIVATE,
                                      DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);
            /*else*/
            //{
                //vector<double> pos;
                //VectorFliedd data(pos, 3);
                                                                //container->appendData("pos", data,
                                      //DECAF_NOFLAG, DECAF_PRIVATE,
                                      //DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);
            /*}*/
            decaf->put(container);
        }

      /* Break conditions of main loop */
      if ((step>maxStep)) break;

   }

   /* Gather all subdomains */
   i = 1;
   for(k=zs[me];k<=ze[me];k++)
   {
    l = 1;
    for(j=ys[me];j<=ye[me];j++)
    {
     for(m=0;m<=xcell-1;m++)
       xtemp[(l-1)*xcell+(i-1)*xcell*ycell+m] = x0[xs[me]+m][j][k];
     l = l+1;
    }
    i = i+1;
   }



    /* Dan comment
   MPI_Gather(xtemp, xcell*ycell*zcell, MPI_DOUBLE, xfinal, xcell*ycell*zcell, MPI_DOUBLE, 0, comm);
    */
   /* Ending time */
   time_final = MPI_Wtime();
   /* Elapsed time */
   elapsed_time = time_final - time_init;

   /* Print results */
   if(me == 0)
   {
    printf("\n");
    printf("  Time step = %3.18f\n",dt);
    printf("\n");
    printf("  Convergence = %11.9f after %d steps\n",epsilon,step);
    printf("\n");
    printf("  Problem size = %d\n",size_x*size_y*size_z);
    printf("\n");
    printf("  Wall Clock = %15.6f\n",elapsed_time);
    printf("\n");
    printf("  Computed solution in outputPar.dat\n");
    printf("\n");

    /* Store solution into output file */
   /*
     file=fopen("outputPar.dat","w");

    for(j=1;j<=size_y+2;j++) {
      for(i=1;i<=size_x+1;i++) {
        fprintf(file,"%15.11f ",temp1_init);
      }
      fprintf(file,"%15.11f\n",temp1_init);
    }

    fprintf(file,"\n");

    for(p=1;p<=size_z;p++) {
      for(v=1;v<=size_x+1;v++)
        fprintf(file,"%15.11f ",temp1_init);
      fprintf(file,"%15.11f\n",temp1_init);
      for(i=1;i<=y_domains;i++) {
        for(j=1;j<=ycell;j++) {
          fprintf(file,"%15.11f ",temp1_init);
          for(k=0;k<=x_domains-1;k++) {
            for(l=0;l<=xcell-1;l++) {
              fprintf(file,"%15.11f ",xfinal[(j-1)*xcell+l+(y_domains-i)*(z_domains*xcell*ycell*zcell)+k*(y_domains*z_domains*xcell*ycell*zcell)+(p-1)*xcell*ycell]);
            }
          }
          fprintf(file,"%15.11f\n",temp1_init);
        }
      }
      for(m=1;m<=size_x+1;m++)
        fprintf(file,"%15.11f ",temp1_init);
      fprintf(file,"%15.11f\n\n",temp1_init);
    }

    for(j=1;j<=size_y+2;j++) {
      for(i=1;i<=size_x+1;i++)
        fprintf(file,"%15.11f ",temp1_init);
      fprintf(file,"%15.11f\n",temp1_init);
    }

    fclose(file);
    */
   }


  /*
   * finalize
   */


   /* Free arrays */
   for (i=0;i<=size_tot_x-1;i++)
   {
    free(x[i]);
    free(x0[i]);
   }
   
   free(x);
   free(x0);
   free(x_all);
   free(x0_all);
   free(xfinal);
   free(xtemp);
   free(xs);
   free(xe);
   free(ys);
   free(ye);
   free(zs);
   free(ze);
   free(buffer);
   /* Free matrices type */
   MPI_Type_free(&matrix_type_oxz);
   MPI_Type_free(&matrix_type_oxy);
   MPI_Type_free(&matrix_type_oyz);

   MPI_Finalize();

}

// gets the atom positions and prints them
void con(Decaf* decaf, int nsteps)
{

     /*
     * @input
     * @param NSTOP
     * @param FILESIZE2PRODUCE
     */

    int lp = N_LP;
    double global_t_analy = 0, t_analy = 0;

    double sum_vx[NMOMENT], sum_vy[NMOMENT];

    /******************** configuration stop ***********/

#ifdef ENABLE_TIMING
    double t1, t2, t3, t4;
    double t_read_1, t_read_2;
    t_read_1 = 0;
    t_read_2 = 0;
                              
    t_analy = 0;
#endif


    int         rank, nprocs;

    int has_keep=0;

    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_DIMES;
    //enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    void * data = NULL;
    uint64_t start[3], count[3];


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


	MPI_Comm comm = decaf->con_comm_handle();
    rank = decaf->con_comm()->rank();

    vector< pConstructData > in_data;

    int step = 0;


    if(rank == 0){
        printf("[nmoment]: application starts\n");
    }


    int timestep = 0;
    //MPI_Barrier(comm);
    double t_start = MPI_Wtime();
	int slice_size = 0;
	double * buffer;
    while (decaf->get(in_data))
    {
        // get the values
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
                run_analysis(buffer, slice_size, lp, sum_vx,sum_vy);

#ifdef V_T
      VT_end(analysis_id);
#endif

                t2 =MPI_Wtime();
                t_analy += t2-t1;

            }
            else
                fprintf(stderr, "[noments]: Error: null pointer in node2\n");
        }

        //printf("[nmoments]: Step %d,t_analy %lf\n", step, t2-t1);

        step +=1;
    }

   printf("[nmoments]: t_analy %lf\n", t_analy);
    //MPI_Barrier(comm);
    double t_end = MPI_Wtime();

    MPI_Reduce(&t_analy, &global_t_analy, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    printf("[noments]: total-start-end %.3f %.3f %.3f\n", t_end- t_start, t_start, t_end);

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    if(rank == 0){
        printf("[nmoments]: max t_analysis %.3f s,now terminates\n", global_t_analy);
    }


    decaf->terminate();
}

// forwards the atom positions in this example
// in a more realistic example, could filter them and only forward some subset of them
void print2(Decaf* decaf)
{
    vector< pConstructData > in_data;
while (decaf->get(in_data))
    {

        // get the values and add them
        for (size_t i = 0; i < in_data.size(); i++)
        {
            fprintf(stderr, "print2 forwarding positions\n");
            decaf->put(in_data[i]);
        }
    }

    // terminate the task (mandatory) by sending a quit message to the rest of the workflow
    fprintf(stderr, "print2 terminating\n");
    decaf->terminate();
}

void run(Workflow& workflow,         // workflow
        int nsteps
        )
     /*    int lammps_nsteps,                  // number of lammps timesteps to execute*/
         //int analysis_interval,              // number of lammps timesteps to skip analyzing
         /*string infile)        */              // lammps input config file
{
    MPI_Init(NULL, NULL);



    Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);

    // run workflow node tasks
    // decaf simply tells the user whether this rank belongs to a workflow node
    // how the tasks are called is entirely up to the user
    // e.g., if they overlap in rank, it is up to the user to call them in an order that makes
    // sense (threaded, alternting, etc.)
    // also, the user can define any function signature she wants
    if (decaf->my_node("prod"))
        prod(decaf, nsteps);
    if (decaf->my_node("con"))
        con(decaf, nsteps);
    if (decaf->my_node("print2"))
        print2(decaf);

    // MPI_Barrier(MPI_COMM_WORLD);

    // cleanup
    delete decaf;

    MPI_Finalize();
}

// normal entry point is run(), called by python
#if 1
int main(int argc,
         char** argv)
{
    Workflow workflow;
    Workflow::make_wflow_from_json(workflow, "heat3d_decaf.json");

    if(argc != 2){
        fprintf(stderr, "[heat3d]: need steps\n");
        return -1;
    }
    int nsteps = atoi(argv[1]);

    // run decaf
    run(workflow, nsteps);

    return 0;

    /*int lammps_nsteps     = 1;*/
    //int analysis_interval = 1;
    //char * prefix         = getenv("DECAF_PREFIX");
    //if (prefix == NULL)
    //{
        //fprintf(stderr, "ERROR: environment variable DECAF_PREFIX not defined. Please export "
                //"DECAF_PREFIX to point to the root of your decaf install directory.\n");
}
#endif
          
