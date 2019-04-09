#include <mpi.h>

//void init_analyze(int nsteps, int size_one);
void analyze(int nsteps, int ii, double *data,  int size_one, /*uint64_t bounds[6],*/ int nlocal, double **msd);
//void finalize_analyze();
