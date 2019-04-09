
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
void doblack ( double w, double** M);
void dored ( double w, double** M );
void exchange ( double** M, int comm[], int rank, MPI_Comm comm_prod);
void iterate ( double w, double** M, int rank, int comm[], MPI_Comm comm_prod, Decaf * decaf);
void setcomm ( int rank, int comm[] );
void setex ( double ex[], double** M, int which );
void initialize_matrix ( double** M );
void unpack ( double** M, int where, double in[] );

/******************************************************************************/

void prod (Workflow& workflow, int nsteps, string infile)

/******************************************************************************/
/*
  Purpose:

    LAPLACE_MPI solves Laplace's equation on a rectangle, using MPI.

  Discussion:

    This program uses a finite difference scheme to solve
    Laplace's equation for a square matrix distributed over a square
    (logical) processor topology.  A complete description of the algorithm
    is found in Fox.

    This program works on the SPMD (single program, multiple data)
    paradigm.  It illustrates 2-d block decomposition, nodes exchanging
    edge values, and convergence checking.

    Each matrix element is updated based on the values of the four
    neighboring matrix elements.  This process is repeated until the data
    converges, that is, until the average change in any matrix element (compared
    to the value 20 iterations previous) is smaller than a specified value.

    To ensure reproducible results between runs, a red/black
    checkerboard algorithm is used.  Each process exchanges edge values
    with its four neighbors.  Then new values are calculated for the upper
    left and lower right corners (the "red" corners) of each node's
    matrix.  The processes exchange edge values again.  The upper right
    and lower left corners (the "black" corners) are then calculated.

    The program is currently configured for a 48x48 matrix
    distributed over four processors.  It can be edited to handle
    different matrix sizes or number of processors, as long as the matrix
    can be divided evenly between the processors.

  Modified:

    14 November 2011

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Reference:

    Geoffrey Fox, Mark Johnson, Gregory Lyzenga, Steve Otto, John Salmon, 
    David Walker,
    Solving Problems on Concurrent Processors,
    Volume 1: General Techniques and Regular Problems, 
    Prentice Hall, 1988,
    ISBN: 0-13-8230226,
    LC: QA76.5.F627.

  Local parameters:

    Local, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.

    Local, FILE *fp, a pointer to the output file.

    Local, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Local, double RESULT[n][n], the results for the complete problem,
    kept by process 0.

    Local, double W, the SOR factor, which must be strictly between 0 and 2.
*/ 
{
  int comm[4];
  FILE *fp;
  int i;
  int j;
  int ntasks;
  int rank;
  double w;
  double wtime;
  int iconf[7];
  double conf[2];
  MPI_Init(NULL, NULL);  

  MPI_Comm comm_world = MPI_COMM_WORLD;
    Decaf* decaf = new Decaf(MPI_COMM_WORLD, workflow);
 MPI_Comm comm_prod = decaf->prod_comm_handle();
  MPI_Comm_rank ( comm_prod, &rank );
  MPI_Comm_size ( comm_prod, &ntasks );

  wtime = MPI_Wtime ( );
   /* Get input parameters */
   if(rank==0)
    readParam(iconf, conf);

   /* Broadcast input parameters */
   MPI_Bcast(iconf,7,MPI_INT,0,comm_prod);
   MPI_Bcast(conf,2,MPI_DOUBLE,0,comm_prod);

   /* Assign input parameters to variables */
   nodeedge    = iconf[0];
   nodeedge    = iconf[1];
   nblock = ntasks;
   nproc = ntasks;
   maxStep   = iconf[6];
   nx = nodeedge * nblock;
   ny = nodeedge;
  double  **M = new double*[nodeedge + 2];
  for(int i = 0; i < nodeedge + 2; ++i) {
    M[i] = new double[nodeedge + 2];
  }
   /* Warning message if dimensions and number of processes don't match */
   if((rank==0) && (ntasks!=(1 * nblock)))
    printf("Number of processes not equal to Number of subdomains\n");

    /*
     * get transport method from env variable
     */


  if ( rank == 0 ) 
  {
    printf ( "\n" );
    printf ( "LAPLACE_MPI:\n" );
    printf ( "  C/MPI version\n" );
    printf ( "  Solve the Laplace equation using MPI.\n" );
  }

  if ( ntasks != nproc )
  {
    if ( rank == 0 ) 
    {
      printf ( "\n" );
      printf ( "Fatal error!\n" );
      printf ( "  MP_PROCS should be set to %i!\n", nproc );
    }
    MPI_Finalize ( );
    exit ( 1 );
  }

  if ( rank == 0 ) 
  {
    printf ( "\n" );
    printf ( "  MPI has been set up.\n" );
  }
/* 
  Initialize the matrix M.
*/
  if ( rank == 0 ) 
  {
    printf ( "  Initialize the matrix M.\n" );
  }
  initialize_matrix ( M );
/* 
  Figure out who I communicate with.
*/
  if ( rank == 0 ) 
  {
    printf ( "  Set the list of neighbors.\n" );
  }
  setcomm ( rank, comm );
/* 
  Update M, using SOR value W, until convergence.
*/
  if ( rank == 0 ) 
  {
    printf ( "  Begin the iteration.\n" );
  }
  w = 1.2;
  iterate ( w, M,  rank, comm, comm_prod, decaf );
/* 
  Report timing 
*/ 
  wtime = MPI_Wtime ( ) - wtime;

  printf ( "  Task %i took %6.3f seconds\n", rank, wtime );
/*
  Write the solution to a file.
*/
/*
  if ( rank == 0 )
  {
    fp = fopen ( "laplace_solution.txt", "w" );

    for ( i = 0; i < n; i++ ) 
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf ( fp, "%f \n", result[i][j] );
      }  
    }
    fclose ( fp );
    printf ( "  Solution written to \"laplace_solution.txt\".\n" );
  }
*/
/*
  Terminate MPI.
*/
  if ( rank == 0 )
  {
    printf ( "\n" );
    printf ( "LAPLACE_MPI:\n" );
    printf ( "  Normal end of execution.\n" );
  }

	decaf->terminate();
}
/******************************************************************************/

void doblack ( double w, double** M )

/******************************************************************************/
/*
  Purpose:

    DOBLACK iterates on the upper right and lower left corners of my matrix.

  Modified:

    16 February 2013

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Input, double W, the SOR factor, which must be strictly between 0 and 2.

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.
*/
{
  int i;
  int j;
/*
  Upper right corner.
*/
  for ( i = 1; i <= nodeedge / 2; i++ )
  {
    for ( j = nodeedge / 2 + 1; j <= nodeedge; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
/*
  Lower left corner.
*/
  for ( i = nodeedge / 2 + 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge / 2; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}
/******************************************************************************/

void dored ( double w, double** M )

/******************************************************************************/   
/*
  Purpose:

    DORED iterates on the upper left and lower right corners of my matrix.

  Modified:

    16 February 2013

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Input, double W, the SOR factor, which must be strictly between 0 and 2.

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.
*/  
{
  int i;
  int j;
/*
  Upper left corner.
*/
  for ( i = 1; i <= nodeedge / 2; i++ )
  {
    for ( j = 1; j <= nodeedge / 2; j++ ) 
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
/*
  Lower right corner.
*/
  for ( i = nodeedge / 2 + 1; i <= nodeedge; i++ )
  {
    for ( j = nodeedge / 2 + 1; j <= nodeedge; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}
/******************************************************************************/

void exchange ( double** M, int comm[], int rank, MPI_Comm comm_prod)

/******************************************************************************/
/*
  Purpose:

   EXCHANGE trades edge values with up to four neighbors.

  Discussion:

    Up to 4 MPI sends are carried out, and up to 4 MPI receives.

  Modified:

    14 November 2011

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Input, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.

    Input, int RANK, the rank of this process.
*/
{
  if(rank==0)
	printf("before new ex and in \n");
  double* ex0 = new double[nodeedge];
  double* ex1 = new double[nodeedge];
  double* ex2 = new double[nodeedge];
  double* ex3 = new double[nodeedge];
  int i;
  double* in0 = new double[nodeedge];
  double* in1 = new double[nodeedge];
  double* in2 = new double[nodeedge];
  double* in3 = new double[nodeedge];
  if(rank==0)
    printf("after new ex and in \n");

  int partner;
  MPI_Request requests[8];
  MPI_Status status[8];
  int tag;
/* 
  Initialize requests.
*/
  for ( i = 0; i < 8; i++ ) 
  {
    requests[i] = MPI_REQUEST_NULL; 
  }
/* 
  Receive from UP neighbor (0).
*/
  if ( comm[0] == 1 )
  {
    partner = rank - nblock;
    tag = 0;
    MPI_Irecv ( in0, nodeedge, MPI_DOUBLE, partner, tag, comm_prod, 
      &requests[0] );
  }
/*
  Receive from RIGHT neighbor (1).
*/
  if ( comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 1;
    MPI_Irecv ( in1, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[1] );
  }
/*
  Receive from DOWN neighbor (2).
*/
  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 2;
    MPI_Irecv ( in2, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[2] );
  }
/*
  Receive from LEFT neighbor (3).
*/
  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 3;
    MPI_Irecv ( in3, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[3] );
  }
/*
  Send up from DOWN (2) neighbor.
*/
  if ( comm[0] == 1 )
  {
    partner = rank - nblock;
    tag = 2;
    setex ( ex0, M, 0 );
    MPI_Isend ( ex0, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[4] );
  }
/*
  Send right form LEFT (3) neighbor.
*/
  if (comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 3;
    setex ( ex1, M, 1 );
    MPI_Isend ( ex1, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[5] );
  }
/*
  Send down from UP (0) neighbor.
*/
  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 0;
    setex ( ex2, M, 2 );
    MPI_Isend ( ex2, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[6] );
  }
/*
  Send left from RIGHT (1) neighbor.
*/
  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 1;
    setex ( ex3, M, 3 );
    MPI_Isend ( ex3, nodeedge, MPI_DOUBLE, partner, tag, comm_prod,
      &requests[7] );
  }
/* 
  Wait for all communication to complete.
*/ 
  MPI_Waitall ( 8, requests, status );
/*
  Copy boundary values, sent by neighbors, into M.
*/
  if ( comm[0] == 1 ) 
  {
    unpack ( M, 0, in0 );
  }
  if ( comm[1] == 1 ) 
  {
    unpack ( M, 1, in1 );
  }
  if ( comm[2] == 1 ) 
  {
    unpack ( M, 2, in2 );
  }
  if ( comm[3] == 1 ) 
  {
    unpack ( M, 3, in3 );
  }

  if(rank==0)
    printf("before delete ex and in \n");

  delete [] ex0;
  delete [] ex1;
  delete [] ex2;
  delete [] ex3;
  delete [] in0;
  delete [] in1;
  delete [] in2;
  delete [] in3;
  if(rank==0)
    printf("after delete ex and in \n");

  return;
}
/******************************************************************************/

void initialize_matrix ( double** M )

/******************************************************************************/
/*
  Purpose:

    INITIALIZE_MATRIX initializes the partial results array M.

  Modified:

    10 January 2012

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Output, double M[nodeedge+2][nodeedge+2], the initialized partial 
    results array.
*/
{
  double avg;
  double bv[4];
  int i;
  int j;

  bv[0] = 100.0;
  bv[1] = 0.0;
  bv[2] = 0.0;
  bv[3] = 0.0;
/* 
  Put the boundary values into M.
*/ 
  for ( i = 1; i <= nodeedge; i++ )
  { 
    M[0][i] =          bv[0];
    M[i][nodeedge+1] = bv[1];
    M[nodeedge+1][i] = bv[2];
    M[i][0] =          bv[3];
  }
/* 
  Set all interior values to be the average of the boundary values.
*/ 
  avg = ( bv[0] + bv[1] + bv[2] + bv[3] ) / 4.0;

  for ( i = 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge; j++ )
    {
      M[i][j] = avg;
    }
  }

  return;
}
/******************************************************************************/

void iterate ( double w, double** M, int rank, 
  int comm[], MPI_Comm comm_prod, Decaf * decaf)

/******************************************************************************/
/*
  Purpose:

    ITERATE controls the iteration, including convergence checking.

  Modified:

    16 February 2013

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Input, double W, the SOR factor, which must be strictly between 0 and 2.

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Output, double RESULT[n][n], the results for the complete problem,
    kept by process 0.

    Input, int RANK, the rank of the process.

    Input, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.

  Local parameters:

    Local, int COUNT, the length, in elements, of messages.

    Local, double DIFF, the average absolute difference in elements
    of M since the last comparison.

    Local, int IT, the iteration counter.

    Local, double MM[n*n], a vector, used to gather the data from
    all processes.  This data is then rearranged into a 2D array.
*/
{
  int count;
  double diff;
  int done;
  double ediff;
  int i;
  double in;
  int index;
  int step;
  int j;
  int k;
  int l;
  //double MM[n*n];
  double  **mold = new double*[nodeedge + 2];
  for(int i = 0; i < nodeedge + 2; ++i) {
    mold[i] = new double[nodeedge + 2];
  }
 
  double  **send = new double*[nodeedge];
  for(int i = 0; i < nodeedge; ++i) {
    send[i] = new double[nodeedge];
  }

  step = 0;
  done = 0;
  for ( i = 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge; j++ )
    {
      mold[i][j] = M[i][j];
    }
  }
  while ( step < maxStep )
  {
/*
  Exchange values with neighbors, update red squares, exchange values
  with neighbors, update black squares.
*/
    exchange ( M, comm, rank, comm_prod );
    dored ( w, M );
    exchange ( M, comm, rank , comm_prod);
    doblack ( w, M );
/*
  Check for convergence every 20 iterations.
  Find the average absolute change in elements of M.
  Maximum iterations is 5000.
*/

    if ( ( done != 1 ) )
    { 
      diff = 0.0;
      for ( i = 1; i <= nodeedge; i++ )
      {
        for ( j = 1; j <= nodeedge; j++ )
        {
          ediff = M[i][j] - mold[i][j];
          if ( ediff < 0.0 ) 
          {
            ediff = - ediff;
          }
          diff = diff + ediff;
          mold[i][j] = M[i][j];
        }
      }
      diff = diff / ( ( double ) ( nodeedge * nodeedge ) );
/*
  IN = sum of DIFF over all processes.
*/
      MPI_Allreduce ( &diff, &in, 1, MPI_DOUBLE, MPI_SUM, comm_prod );
/*
      if ( in < ( double ) nproc * 0.001 ) 
      {
        done = 1;
      }
*/
    }


    // output dir for adios
    char filename [256];
	char *filepath = getenv("BP_DIR");
	int ys,xs;
	ys = 0;
	xs = 0;
	xs = (rank % nblock) * nodeedge;
//	ys = (nblock - rank / 3 -1) * nodeedge;


	/* decaf put */
	if (1)
	{
		pConstructData container;

		// lammps gathered all positions to rank 0
		//if (decaf->prod_comm()->rank() == 0)
			// debug
			//         for (int i = 0; i < 10; i++)         // print first few atoms
			//           fprintf(stderr, "%.3lf %.3lf %.3lf\n",
			// x[3 * i], x[3 * i + 1], x[3 * i + 2]);

		VectorFliedd data(&M[1][1], nodeedge * nodeedge, SIZE_ONE);

		container->appendData("pos", data,
								  DECAF_NOFLAG, DECAF_PRIVATE,
								  DECAF_SPLIT_DEFAULT, DECAF_MERGE_DEFAULT);
		if (rank == 0)
		{
			printf("[nmoments]: after laplace append data into decaf\n");
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
			printf("after laplace put data into decaf\n");
		}
	}


	step ++;
  }
/* 
  Send results to task 0.
*/ 
  for ( i = 0; i < nodeedge; i++ )
  {
    for ( j = 0; j < nodeedge; j++ )
    {
      send[i][j] = M[i+1][j+1];
    }
  }

for(int i = 0; i < nodeedge + 2; ++i) {
    delete [] mold[i];
}
delete [] mold;
for(int i = 0; i < nodeedge; ++i) {
    delete [] send[i];
}
delete [] send;
/*  count = nodeedge * nodeedge;

  MPI_Gather ( &send, count, MPI_DOUBLE, &MM, count, MPI_DOUBLE, 0, 
    MPI_COMM_WORLD );

  printf ( "  ITERATE gathered updated results to process 0.\n" );
*/
  return;
}
/******************************************************************************/

void setcomm ( int rank, int comm[] )

/******************************************************************************/
/*
  Purpose:

    SETCOMM determines the active communication directions.

  Discussion:

    In this picture, we're assuming the RESULTS array is split among 
    four processes, numbered 0 through 3 and arranged as suggested by the 
    following:

        0  |  1
     ------+-------
        2  |  3

    Then process 0 must communicate with processes 1 and 2 only,
    so its COMM array would be { 0, 1, 1, 0 }.

  Modified:

    14 November 2011

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Input, int RANK, the rank of the process.

    Output, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.
*/
{
  int i;
    comm[0] = 0;
    comm[2] = 0;
    if(rank == 0)
        comm[3] = 0;
    else
        comm[3] = 1;
    if(rank == nproc - 1)
        comm[1] = 0;
    else
        comm[1] = 1;


  return;
}
/******************************************************************************/

void setex ( double ex[], double** M, int which )

/******************************************************************************/
/*
  Purpose:

    SETEX pulls off the edge values of M to send to another task.

  Modified:

    14 November 2011

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Output, double EX[NODEEDGE], the values to be exchanged.

    Input, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process. 

    Input, int WHICH, 0, 1, 2, or 3, indicates the edge from which
    the data is to be copied.
*/                  
{
  int i;

  switch ( which ) 
  {
    case 0:
    {
      for ( i = 1; i <= nodeedge; i++) 
      {
        ex[i-1] = M[1][i];
      }
      break;
    }
    case 1:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[i][nodeedge];
      }
      break;
    }
    case 2:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[nodeedge][i];
      }
      break;
    }
    case 3:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[i][1];
      }
      break;
    }
  }
  return;
}
/******************************************************************************/

void unpack ( double** M, int where, double in[] )

/******************************************************************************/
/*
  Purpose:

    UNPACK puts the vector of new edge values into the edges of M.

  Modified:

    14 November 2011

  Author:

    Sequential C version by Robb Newman.
    MPI C version by Xianneng Shen.

  Parameters:

    Output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Input, int WHERE, 0, 1, 2, or 3, indicates the edge to which the 
    data is to be applied.

    Input, int IN[nodeedge], the boundary data.
*/
{
  int i;

  if ( where == 0 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[0][i+1] = in[i]; 
    }
  }
  else if ( where == 1 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[i+1][nodeedge+1] = in[i];
    }
  }
  else if ( where == 2 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[nodeedge+1][i+1] = in[i];
    }
  }
  else if ( where == 3 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[i+1][0] = in[i];
    }
  }

  return;
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

    prod(workflow, nsteps, infile);

           
    return 0;
}


