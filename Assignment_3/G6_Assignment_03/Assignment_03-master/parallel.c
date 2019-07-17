#include "parallel.h"
#include <mpi.h>


void Program_Message(char *txt)
/* produces x_dim stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces x_dim stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce x_dim text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}



void init_parallel (int iproc,int jproc,int imax,int jmax,
int *myrank,int *il,int *ir,int *jb,int *jt,
int *rank_l,int *rank_r,int *rank_b,int *rank_t,
int *omg_i,int *omg_j,int num_proc)
{
  
	int i_per_iproc, j_per_jproc, i_rem, j_rem;

	/* Set sub-domain coordinates */
	*omg_i = ((*myrank) % iproc) + 1;
	*omg_j = (((*myrank) - (*omg_i) + 1) / iproc) + 1;

	/* Compute il, ir for each sub-domain*/
	i_per_iproc = (imax / iproc);
	j_per_jproc = (jmax / jproc);
	i_rem = (imax % iproc);
	j_rem = (jmax % jproc);

	/* for il ir*/
	if ((*omg_i) == 1)   /* to rank zero assign the remainder*/
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + 1;
	}
	else
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + i_rem + 1;
	}
	*ir = ((*omg_i) * i_per_iproc) + i_rem;

	/* for jb jt*/
	if ((*omg_j) == 1)
	{
	  *jb = (((*omg_j) -1) * j_per_jproc) + 1;
	}
	else
	{
	  *jb = (((*omg_j) -1) * j_per_jproc) + j_rem + 1;
	}
	*jt = ((*omg_j) * j_per_jproc) + j_rem;

    /* Assign rank of neighbour to each sub-domain: rank_l, rank_r, rank_b, rank_t*/
	/* Left boundary*/
	if ((*il) == 1)
	{
	  *rank_l = MPI_PROC_NULL;
	}
	else
	{
	  *rank_l = (*myrank) - 1;
	}

	/* Right boundary*/
	if ((*ir) == imax)
	{
	  *rank_r = MPI_PROC_NULL;
	}
	else
	{
	  *rank_r = (*myrank) + 1;
	}

	/* Bottom boundary*/
	if ((*jb) == 1)
	{
	  *rank_b = MPI_PROC_NULL;
	}
	else
	{
	  *rank_b = (*myrank) - iproc;
	}

	/* Top boundary*/
	if ((*jt) == jmax)
	{
	  *rank_t = MPI_PROC_NULL;
	}
	else
	{
	  *rank_t = (*myrank) + iproc;
}

}



void pressure_comm(double **P,int il,int ir,int jb,int jt
int rank_l,int rank_r,int rank_b,int rank_t,
double *bufSend,double *bufRecv, MPI_Status *status, int chunk)
{
	MPI_Status status;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb +1;

	bufSend = (double*)malloc(y_dim*sizeof(double));
	bufRecv = (double*)malloc(y_dim*sizeof(double));

	if (rank_l != MPI_PROC_NULL) // send to left
	{
		int counter=0;
		for(j=1; j<=y_dim; j++)
		{
			bufSend[counter] = P[1][j];
			counter++;
		}
		//bufSend[0] = P[1][1];
		MPI_Send( bufSend, y_dim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
	}


	if (rank_r != MPI_PROC_NULL) //recieve from right and then send to right
	{
		MPI_Recv( bufRecv, y_dim, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, &status );
		counter=0;
		for(j=1; j<=y_dim; j++)
		{
			P[x_dim+1][j] = bufRecv[counter];
			counter++;
		}
		//bufRecv[0] = P[ir+1][1];

		counter=0;
		for(j=1; j<=y_dim; j++)
		{
			bufSend[counter] = P[x_dim][j];
			counter++;
		}
		//bufSend[0] = P[ir][1];
		MPI_Send( bufSend, y_dim, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
	}


	if (rank_l != MPI_PROC_NULL) //recieving from left 
	{
		MPI_Recv(bufRecv, y_dim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, &status );
		counter=0;
		for(j=1; j<=y_dim; j++)
		{
			P[0][j] = bufRecv[counter]; 
			counter++;
		}
		//bufRecv[0] = P[0][1];

	}

	//destroy first
	free(bufSend);
	free(bufRecv);

	bufSend = (double*)malloc(x_dim*sizeof(double));
	bufRecv = (double*)malloc(x_dim*sizeof(double));

	if (rank_t != MPI_PROC_NULL) //send to the top
	{
		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			bufSend[counter] = P[i][y_dim];
			counter++;
		}
		MPI_Send( bufSend, x_dim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
	}


	if (rank_b != MPI_PROC_NULL) //recieving from the bottom and then send to the bottom
	{
		MPI_Recv( bufRecv, x_dim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, &status );
		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			P[i][0] =bufRecv[counter];
			counter++;
		}

		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			bufSend[counter] = P[i][1];
			counter++;
		}
		MPI_Send( bufSend, x_dim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
	}

	if (rank_t != MPI_PROC_NULL) //recieving from the top
	{
		MPI_Recv( bufRecv, x_dim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, &status );
		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			P[i][y_dim+1] = bufRecv[counter];
			counter++;
		}
	}
	free(bufSend);
	free(bufRecv);
}




void uv_comm (double**U,double**V,int il,int ir, int jb,int jt,
int rank_l,int rank_r,int rank_b,int rank_t,
double *bufSend,double *bufRecv, MPI_Status *status, int chunk)
{
	`
	int x_dim = ir - il + 1;
	int y_dim = jt - jb +1;

	bufSend = (double*)malloc(2*y_dim*sizeof(double));
	bufRecv = (double*)malloc(2*y_dim*sizeof(double));

	MPI_Status status;

	if (rank_l != MPI_PROC_NULL) //sending to left
	{
		int counter=0;
		for(j=1; j<=y_dim; j++)
		{
			bufSend[counter] = U[2][j];
			counter++;
		}

		for(j=1; j<=y_dim; j++)
		{
			bufSend[counter] = V[1][j];
			counter++;
		}
		MPI_Send( bufSend, 2*y_dim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );

	}


	if (rank_r != MPI_PROC_NULL) //recieving from right and then sending to right
	{
		MPI_Recv( bufRecv, 2*y_dim, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, &status );
		counter=0;
		for(j=1; j<=y_dim; j++)
		{
			U[x_dim+1][j] = bufRecv[counter];
			counter++;
		}

		for(j=1; j<=y_dim; j++)
		{
			V[x_dim+1][j] = bufRecv[counter];
			counter++;
		}
		counter=0;
		for(j=1; j<=y_dim; j++)
		{
			bufSend[counter] = U[x_dim-1][j];
			counter++;
		}
		for(j=1; j<=y_dim; j++)
		{
			bufSend[counter] = V[x_dim][j];
			counter++;
		}
		MPI_Send( bufSend, 2*y_dim, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
	}


	if (rank_l != MPI_PROC_NULL) //recieving from left
	{
		MPI_Recv(bufRecv, 2*y_dim, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, &status );

		counter=0;
		for(j=1; j<=y_dim; j++)
		{
			U[0][j] = bufRecv[counter] ;
			counter++;
		}
		for(j=1; j<=y_dim; j++)
		{
			bufRecv[counter] = V[0][j];
			counter++;
		}

	}

	//destroy first
	free(bufSend);
	free(bufRecv);

	bufSend = (double*)malloc(2*x_dim*sizeof(double));
	bufRecv = (double*)malloc(2*x_dim*sizeof(double));

	if (rank_t != MPI_PROC_NULL) //send to top
	{
		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			bufSend[counter] = U[i][y_dim];
			counter++;
		}

		for(i=1; i<=x_dim; i++)
		{
			bufSend[counter] = V[i][y_dim-1];
			counter++;
		}

		MPI_Send( bufSend, 2*x_dim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
	}


	if (rank_b != MPI_PROC_NULL) //recieve from bottom and then send to the bottom
	{
		MPI_Recv( bufRecv, 2*x_dim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, &status );
		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			U[i][0] = bufRecv[counter];
			counter++;
		}

		for(i=1; i<=x_dim; i++)
		{
			V[i][0] = bufRecv[counter];
			counter++;
		}


		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			bufSend[counter] = U[i][1];
			counter++;
		}

		for(i=1; i<=x_dim; i++)
		{
			bufSend[counter] = V[i][2];
			counter++;
		}
		MPI_Send( bufSend, 2*x_dim, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );

	}

	if (rank_t != MPI_PROC_NULL) //recieving from top
	{

		MPI_Recv( bufRecv, 2*x_dim, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, &status );

		counter=0;
		for(i=1; i<=x_dim; i++)
		{
			U[i][y_dim+1] = bufRecv[counter];
			counter++;
		}

		for(i=1; i<=x_dim; i++)
		{
			V[i][y_dim+1] = bufRecv[counter];
			counter++;
		}

	}
	free(bufSend);
	free(bufRecv);

}


