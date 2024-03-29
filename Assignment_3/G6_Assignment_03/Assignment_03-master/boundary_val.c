#include "boundary_val.h"

void boundaryvalues(int imax,int jmax,double **U,double **V, int rank_l, int rank_r, int rank_b, int rank_t)
{
	int i, j;

	/* For bottom boundary */
	if (rank_b == MPI_PROC_NULL){
		for(i = 1; i < imax+1; i++){
			U[i][0]	= -U[i][1];
			V[i][0]	= 0;
		}
	}

	/* For top boundary */
	if (rank_t == MPI_PROC_NULL){
		for(i = 1; i < imax+1; i++){
			U[i][jmax+1] = 2.0-U[i][jmax];
			V[i][jmax]	= 0;
		}
	}

	/* for left boundary*/
	if (rank_l == MPI_PROC_NULL){
		for(j = 1; j < jmax+1; j++){
			U[0][j] = 0;
			V[0][j]	= -V[1][j];
		}
	}

	/* for right wall */
	if (rank_r == MPI_PROC_NULL){
		for(j = 1; j < jmax+1; j++){
			U[imax][j] = 0;
			V[imax+1][j]	= -V[imax][j];
		}
	}
}

