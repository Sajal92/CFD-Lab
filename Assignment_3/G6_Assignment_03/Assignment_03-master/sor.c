#include "sor.h"
#include "math.h"
#include <mpi.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  double *bufSend ;
	double *bufRecv;
	MPI_Status status;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;

	/* Set left & right global domain boundaries according to Neumann boundary conditions */
	if(rank_l == MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Only receive/send data from/to right */
	{
		for(j = 1; j <= y_dim; j++) {
			P[0][j] = P[1][j];
		}
	}
	else if(rank_l != MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /* Only send/receive data to/from left */
	{
		for(j = 1; j <= y_dim; j++) {
			P[x_dim+1][j] = P[x_dim][j];
		}
	}
	else if(rank_l == MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /* No bordering processes */
	{
		for(j = 1; j <= y_dim; j++) {
			P[0][j] = P[1][j];
			P[x_dim+1][j] = P[x_dim][j];
		}
	}

	/* Set top & bottom global domain boundaries according to Neumann boundary conditions */
	if(rank_t == MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /* Only receive/send data from/to bottom */
	{
		for(i = 1; i <= x_dim; i++) {
			P[i][y_dim+1] = P[i][y_dim];
		}
	}
	else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* Only send/receive data to/from top */
	{
		for(i = 1; i <= x_dim; i++) {
			P[i][0] = P[i][1];
		}
	}
	else if(rank_t == MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* No bordering processes */
	{
		for(i = 1; i <= x_dim; i++) {
			P[i][0] = P[i][1];
			P[i][y_dim+1] = P[i][y_dim];
		}
	}

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }


  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  *res= rloc;

}
