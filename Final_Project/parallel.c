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



void init_parallel (int iproc,int jproc,int kproc,int imax,int jmax,int kmax,
					int *myrank,int *il,int *ir,int *jb,int *jt,int *kf, int *kb,
					int *rank_l,int *rank_r,int *rank_b,int *rank_t, int *rank_f, int *rank_kb,
					int *omg_i,int *omg_j,int *omg_k,int num_proc)
{
  
	int i_per_iproc, j_per_jproc,k_per_kproc, i_rem, j_rem,k_rem;

	/* Set sub-domain coordinates */
	*omg_i = ((*myrank) % iproc) + 1;
	*omg_j = (((*myrank) - (*omg_i) + 1) / iproc) + 1;
	//*omg_k = 



	/* Compute il, ir for each sub-domain*/
	i_per_iproc = (imax / iproc);
	j_per_jproc = (jmax / jproc);
	k_per_kproc = (kmax / kproc);
	i_rem = (imax % iproc);
	j_rem = (jmax % jproc);
	k_rem = (kmax % kproc);



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
	/* for jb jt*/
	if ((*omg_k) == 1)
	{
	  *kf = (((*omg_k) -1) * k_per_kproc) + 1;
	}
	else
	{
	  *kf = (((*omg_k) -1) * k_per_kproc) + k_rem + 1;
	}
	*kb = ((*omg_k) * k_per_kproc) + k_rem;



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
	/* Front boundary*/
	if ((*kf) == 1)
	{
	  *rank_f = MPI_PROC_NULL;
	}
	else
	{
	  *rank_f = (*myrank) - kproc;
	}
	/* Back boundary*/
	if ((*kb) == kmax)
	{
	  *rank_kb = MPI_PROC_NULL;
	}
	else
	{
	  *rank_kb = (*myrank) + kproc;
	}

}

void flag_comm(	double ***flag,
					int il, int ir, 
					int jb, int jt,
					int kf, int kb, 
					int rank_l, int rank_r,
					int rank_b, int rank_t,
					int rank_f, int rank_kb, 
					double *bufSend,double *bufRecv,
					MPI_Status *status, int chunk)
/* This method exchanges pressure values between processes that treat adjacent sub-domains */ 
{
	int counter;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;
	int z_dim = kb - kf + 1;

	/* ---------------- flag: LEFT & RIGHT COMMUNICATION BEGINS ------------- */

	/* Step a) Send to left - receive from right */
	bufSend = (double*)malloc((y_dim*z_dim)*sizeof(double));
	bufRecv = (double*)malloc((y_dim*z_dim)*sizeof(double));

	if (rank_l != MPI_PROC_NULL) // Send to left
	{
		counter=0;
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = flag[1][j][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (y_dim*z_dim), MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
	}


	if (rank_r != MPI_PROC_NULL) // Receive from right
	{
		MPI_Recv( bufRecv, (y_dim*z_dim), MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				flag[x_dim+1][j][k] = bufRecv[counter];
				counter++;
			}
		}

		/* Step b) Send to right and receive from left */
		counter=0;  // Send to right
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = flag[x_dim][j][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (y_dim*z_dim), MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
	}


	if (rank_l != MPI_PROC_NULL) // Receive from left 
	{
		MPI_Recv(bufRecv, (y_dim*z_dim), MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				flag[0][j][k] = bufRecv[counter]; 
				counter++;
			}
		}
	}

	//destroy first
	free(bufSend);
	free(bufRecv);
	/* ---------------- flag: LEFT & RIGHT COMMUNICATION ENDS ------------- */



	/* ---------------- flag: TOP & BOTTOM COMMUNICATION BEGINS ------------- */
	bufSend = (double*)malloc((x_dim*z_dim)*sizeof(double));
	bufRecv = (double*)malloc((x_dim*z_dim)*sizeof(double));


	/* Step c) Send to top - receive from bottom */
	if (rank_t != MPI_PROC_NULL) // send to the top
	{
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = flag[i][y_dim][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*z_dim), MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
	}


	if (rank_b != MPI_PROC_NULL) // Recieve from the bottom 	
	{
		MPI_Recv( bufRecv, (x_dim*z_dim), MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				flag[i][0][k] =bufRecv[counter];
				counter++;
			}
		}

		/* Step d) Send to the bottom - receive from the top */
		counter=0; // Send to the bottom
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = flag[i][1][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*z_dim), MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
	}

	if (rank_t != MPI_PROC_NULL) // Receive from the top
	{
		MPI_Recv( bufRecv, (x_dim*z_dim), MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				flag[i][y_dim+1][k] = bufRecv[counter];
				counter++;
			}
		}
	}
	free(bufSend);
	free(bufRecv);
	/* ----------------flag: TOP & BOTTOM COMMUNICATION ENDS ------------- */

	/* ---------------- flag: FRONT & BACK COMMUNICATION BEGINS ------------- */
	bufSend = (double*)malloc((x_dim*y_dim)*sizeof(double));
	bufRecv = (double*)malloc((x_dim*y_dim)*sizeof(double));


	/* Step e) Send to front - receive from back */
	if (rank_f != MPI_PROC_NULL) // send to the front
	{
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				bufSend[counter] = flag[i][j][z_dim];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*y_dim), MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD );
	}


	if (rank_kb != MPI_PROC_NULL) // Recieve from the back	
	{
		MPI_Recv( bufRecv, (x_dim*y_dim), MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				flag[i][j][0] =bufRecv[counter];
				counter++;
			}
		}

		/* Step f) Send to the back - receive from the front */
		counter=0; // Send to the back
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				bufSend[counter] = flag[i][j][1];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*y_dim), MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD );
	}

	if (rank_f != MPI_PROC_NULL) // Receive from the front
	{
		MPI_Recv( bufRecv, (x_dim*y_dim), MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				flag[i][j][z_dim+1] = bufRecv[counter];
				counter++;
			}
		}
	}
	free(bufSend);
	free(bufRecv);
	/* ----------------flag: FRONT & BACK COMMUNICATION ENDS ------------- */
}



void pressure_comm(	double ***P,
					int il, int ir, 
					int jb, int jt,
					int kf, int kb, 
					int rank_l, int rank_r,
					int rank_b, int rank_t,
					int rank_f, int rank_kb, 
					double *bufSend,double *bufRecv,
					MPI_Status *status, int chunk)
/* This method exchanges pressure values between processes that treat adjacent sub-domains */ 
{
	int counter;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;
	int z_dim = kb - kf + 1;

	/* ---------------- Pressure: LEFT & RIGHT COMMUNICATION BEGINS ------------- */

	/* Step a) Send to left - receive from right */
	bufSend = (double*)malloc((y_dim*z_dim)*sizeof(double));
	bufRecv = (double*)malloc((y_dim*z_dim)*sizeof(double));

	if (rank_l != MPI_PROC_NULL) // Send to left
	{
		counter=0;
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = P[1][j][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (y_dim*z_dim), MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD );
	}


	if (rank_r != MPI_PROC_NULL) // Receive from right
	{
		MPI_Recv( bufRecv, (y_dim*z_dim), MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				P[x_dim+1][j][k] = bufRecv[counter];
				counter++;
			}
		}

		/* Step b) Send to right and receive from left */
		counter=0;  // Send to right
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = P[x_dim][j][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (y_dim*z_dim), MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD );
	}


	if (rank_l != MPI_PROC_NULL) // Receive from left 
	{
		MPI_Recv(bufRecv, (y_dim*z_dim), MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int j=1; j<=y_dim; j++){
			for(int k=1; k<=z_dim; k++){
				P[0][j][k] = bufRecv[counter]; 
				counter++;
			}
		}
	}

	//destroy first
	free(bufSend);
	free(bufRecv);
	/* ---------------- Pressure: LEFT & RIGHT COMMUNICATION ENDS ------------- */



	/* ---------------- Pressure: TOP & BOTTOM COMMUNICATION BEGINS ------------- */
	bufSend = (double*)malloc((x_dim*z_dim)*sizeof(double));
	bufRecv = (double*)malloc((x_dim*z_dim)*sizeof(double));


	/* Step c) Send to top - receive from bottom */
	if (rank_t != MPI_PROC_NULL) // send to the top
	{
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = P[i][y_dim][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*z_dim), MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD );
	}


	if (rank_b != MPI_PROC_NULL) // Recieve from the bottom 	
	{
		MPI_Recv( bufRecv, (x_dim*z_dim), MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				P[i][0][k] =bufRecv[counter];
				counter++;
			}
		}

		/* Step d) Send to the bottom - receive from the top */
		counter=0; // Send to the bottom
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				bufSend[counter] = P[i][1][k];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*z_dim), MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD );
	}

	if (rank_t != MPI_PROC_NULL) // Receive from the top
	{
		MPI_Recv( bufRecv, (x_dim*z_dim), MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int k=1; k<=z_dim; k++){
				P[i][y_dim+1][k] = bufRecv[counter];
				counter++;
			}
		}
	}
	free(bufSend);
	free(bufRecv);
	/* ----------------Pressure: TOP & BOTTOM COMMUNICATION ENDS ------------- */

	/* ---------------- Pressure: FRONT & BACK COMMUNICATION BEGINS ------------- */
	bufSend = (double*)malloc((x_dim*y_dim)*sizeof(double));
	bufRecv = (double*)malloc((x_dim*y_dim)*sizeof(double));


	/* Step e) Send to front - receive from back */
	if (rank_f != MPI_PROC_NULL) // send to the front
	{
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				bufSend[counter] = P[i][j][z_dim];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*y_dim), MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD );
	}


	if (rank_kb != MPI_PROC_NULL) // Recieve from the back	
	{
		MPI_Recv( bufRecv, (x_dim*y_dim), MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				P[i][j][0] =bufRecv[counter];
				counter++;
			}
		}

		/* Step f) Send to the back - receive from the front */
		counter=0; // Send to the back
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				bufSend[counter] = P[i][j][1];
				counter++;
			}
		}
		MPI_Send( bufSend, (x_dim*y_dim), MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD );
	}

	if (rank_f != MPI_PROC_NULL) // Receive from the front
	{
		MPI_Recv( bufRecv, (x_dim*y_dim), MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD, status );
		counter=0;
		for(int i=1; i<=x_dim; i++){
			for(int j=1; j<=y_dim; j++){
				P[i][j][z_dim+1] = bufRecv[counter];
				counter++;
			}
		}
	}
	free(bufSend);
	free(bufRecv);
	/* ----------------Pressure: FRONT & BACK COMMUNICATION ENDS ------------- */
}




void uv_comm (	double ***U, double ***V, double ***W,
				int il,int ir,
				int jb, int jt, 
				int kf, int kb,
				int rank_l, int rank_r,
				int rank_b, int rank_t,
				int rank_f, int rank_kb,
				double *bufSend, double *bufRecv, 
				MPI_Status *status, int chunk)
/* This method exchanges velocity values between processes that treat adjacent sub-domains */ 
{
	int size;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;
	int z_dim = kb - kf + 1;
	int count;
	count =0; 

	/* ---------------- U-V-W Velocity: LEFT & RIGHT COMMUNICATION BEGINS ------------- */

	/* Send to the left, receive from the right */
	/* Send to the right, receive from the left */
	if(rank_l != MPI_PROC_NULL || rank_r != MPI_PROC_NULL){
		size = (2 * y_dim) + 1;

		if(rank_l != MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Perform both left-right transfers */
		{
			
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy left values to send */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[2][j][k];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[1][j][k];	
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = W[1][j][k];
					count++;
				}
			}

			/* Send left values, receive right values */
			MPI_Sendrecv(bufSend, count, MPI_DOUBLE, rank_l, 1, bufRecv, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received right values */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
		       			U[x_dim +1 ][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){	
					V[x_dim + 1][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){	
					W[x_dim + 1][j][k] = bufRecv[count];
					count++;
				}
			}
			count = 0;
			/* Copy right values to send */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[x_dim-1][j][k];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[x_dim][j][k];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = W[x_dim][j][k];
					count++;
				}
			}

			/* Send right values, receive left values */
			MPI_Sendrecv(bufSend, count, MPI_DOUBLE, rank_r, 1, bufRecv, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received left values */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					U[0][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= size; j++){
				for(int k=1; k<=z_dim; k++){
					V[0][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= size; j++){
				for(int k=1; k<=z_dim; k++){
					W[0][j][k] = bufRecv[count];
					count++;
				}
			}
			
			free(bufSend);
			free(bufRecv);
			count=0;
		}
		else if(rank_l == MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Receive data from right and  send data to right */
		{
			
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			
			count = 0;
			/* Copy right values to send */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[x_dim-1][j][k];
					count++;
				}
			}
			for(int j = 1 ; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[x_dim][j][k];
					count++;
				}
			}
			

			/* Send right values */
			MPI_Send(bufSend, count, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
			
			/* Receive right values */
			MPI_Recv(bufRecv, count, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received right values */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					U[x_dim + 1][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					V[x_dim + 1][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					W[x_dim + 1][j][k] = bufRecv[count];
					count++;
				}
			}
			free(bufRecv);
			free(bufSend);
			count = 0;
		}
		else if(rank_l != MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /*Send data to left and Receive data from left*/
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			count = 0;
			/* Copy left values to send */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[2][j][k];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[1][j][k];
					count++;
				}
			}
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = W[1][j][k];
					count++;
				}
			}

			/* Send left values */
			MPI_Send(bufSend, count, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);

			/* Receive left values */
			MPI_Recv(bufRecv, count, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received left values */
			for(int j = 1; j <= y_dim; j++){
				for(int k=1; k<=z_dim; k++){
					U[0][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= size; j++){
				for(int k=1; k<=z_dim; k++){
					V[0][j][k] = bufRecv[count];
					count++;
				}
			}
			for(int j = 1; j <= size; j++){
				for(int k=1; k<=z_dim; k++){
					W[0][j][k] = bufRecv[count];
					count++;
				}
			}
			free(bufRecv);
			free(bufSend);
			count = 0;
		}
	}
	/* ---------------- U-V-W Velocity: LEFT & RIGHT COMMUNICATION ENDS ------------- */



	/* ---------------- U-V-W Velocity: TOP & BOTTOM COMMUNICATION BEGINS ------------- */
	/* Send to the top, receive from the bottom */
	/* Send to the bottom, receive from the top */
	if(rank_t != MPI_PROC_NULL || rank_b != MPI_PROC_NULL)
	{
		size = (2 * x_dim) + 1;

		if(rank_t != MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /* Perform both top-bottom transfers */
		{
			count = 0;
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy top values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[i][y_dim-1][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[i][y_dim][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = W[i][y_dim][k];
					count++;
				}
			}

			/* Send values to top, receive values from bottom*/
			MPI_Sendrecv(bufSend, count, MPI_DOUBLE, rank_t, 1, bufRecv, size, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);

			/* Copy received bottom values */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					V[i][0][k] = bufRecv[count];	
					count++;
				}
			}
			for(int i = 1; i <= size; i++){
				for(int k=1; k<=z_dim; k++){
					U[i][0][k] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= size; i++){
				for(int k=1; k<=z_dim; k++){
					W[i][0][k] = bufRecv[count];
					count++;
				}
			}
			count = 0;
			/* Copy bottom values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[i][1][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[i][1][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = W[i][1][k];
					count++;
				}
			}

			/* Send bottom values, receive top values */
			MPI_Sendrecv(bufSend, count, MPI_DOUBLE, rank_b, 1, bufRecv, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received top values */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					V[i][y_dim+1][k] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					U[i][y_dim + 1][k] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					W[i][y_dim + 1][k] = bufRecv[count];
					count++;
				}
			}

			free(bufSend);
			free(bufRecv);
			count = 0;
		}
		else if(rank_t == MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /*Receive data from bottom and  send data to bottom*/
		{
			count = 0;
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			
			/* Copy bottom values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[i][1][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[i][1][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = W[i][1][k];
					count++;
				}
			}

			/* Send bottom values */
			MPI_Send(bufSend, count, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);
			/* Receive bottom values */
			MPI_Recv(bufRecv, count, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received bottom values */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					V[i][0][k] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= size; i++){
				for(int k=1; k<=z_dim; k++){
					U[i][0][k] = bufRecv[count];
					count++;
				}
			}
			count = 0;
			free(bufRecv);
			free(bufSend);
		}
		else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /*Send data to top and receive data from top */
		{
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			count = 0;
			/* Copy top values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = V[i][y_dim-1][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[i][y_dim][k];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					bufSend[count] = U[i][y_dim][k];
					count++;
				}
			}

			/* Send top values */
			MPI_Send(bufSend, count, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);

			/* Receive top values */
			MPI_Recv(bufRecv, count, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received top values */
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					V[i][y_dim+1][k] = bufRecv[count];
					count++;
				}		
			}
			for(int i = 1; i <= x_dim; i++){
				for(int k=1; k<=z_dim; k++){
					U[i][y_dim + 1][k] = bufRecv[count];
					count++;
				}
			}
			free(bufRecv);
			free(bufSend);
			count = 0;
		}
	}


/* ---------------- U-V-W Velocity: TOP & BOTTOM COMMUNICATION ENDS ------------- */










/* ---------------- U-V-W Velocity: Front & Back COMMUNICATION BEGINS ------------- */

	/* Send to the front, receive from the back */
	/* Send to the back, receive from the front */
	if(rank_f != MPI_PROC_NULL || rank_kb != MPI_PROC_NULL)	{
		size = (2 * y_dim) + 1;
		count = 0;
		if(rank_f != MPI_PROC_NULL && rank_kb != MPI_PROC_NULL)  /* Perform both front-back transfers */
		{
			
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy left values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = U[i][j][2];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = V[i][j][1];	
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = U[i][j][1];
					count++;
				}
			}

			/* Send front values, receive back values */
			MPI_Sendrecv(bufSend, count, MPI_DOUBLE, rank_f, 1, bufRecv, size, MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received back values */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
		       			U[i][j][z_dim +1] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){	
					V[i][j][z_dim +1] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){	
					W[i][j][z_dim +1] = bufRecv[count];
					count++;
				}
			}
			count = 0;
			/* Copy back values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = U[i][j][z_dim-1];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = V[i][j][z_dim];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = W[i][j][z_dim];
					count++;
				}
			}

			/* Send back values, receive front values */
			MPI_Sendrecv(bufSend, count, MPI_DOUBLE, rank_kb, 1, bufRecv, size, MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received front values */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					U[i][j][0] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					V[i][j][0] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					W[i][j][0] = bufRecv[count];
					count++;
				}
			}
			
			free(bufSend);
			free(bufRecv);
			count=0;
		}
		else if(rank_f == MPI_PROC_NULL && rank_kb != MPI_PROC_NULL)  /* Receive data from back and  send data to back */
		{
			
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			
			count = 0;
			/* Copy right values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = U[i][j][x_dim-1];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = V[i][j][x_dim];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = W[i][j][x_dim];
					count++;
				}
			}
			

			/* Send back values */
			MPI_Send(bufSend, count, MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD);
			
			/* Receive back values */
			MPI_Recv(bufRecv, count, MPI_DOUBLE, rank_kb, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received back values */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					U[i][j][z_dim + 1] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					V[i][j][z_dim + 1] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					W[i][j][z_dim + 1] = bufRecv[count];
					count++;
				}
			}
			free(bufRecv);
			free(bufSend);
			count = 0;
		}
		else if(rank_f != MPI_PROC_NULL && rank_kb == MPI_PROC_NULL)  /*Send data to front and Receive data from front*/
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			count = 0;
			/* Copy front values to send */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = U[i][j][2];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = V[i][j][1];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					bufSend[count] = W[i][j][1];
					count++;
				}
			}

			/* Send front values */
			MPI_Send(bufSend, count, MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD);

			/* Receive front values */
			MPI_Recv(bufRecv, count, MPI_DOUBLE, rank_f, 1, MPI_COMM_WORLD, status);
			count = 0;
			/* Copy received front values */
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					U[i][j][0] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					V[i][j][0] = bufRecv[count];
					count++;
				}
			}
			for(int i = 1; i <= x_dim; i++){
				for(int j=1; j<=y_dim; j++){
					W[i][j][0] = bufRecv[count];
					count++;
				}
			}
			free(bufRecv);
			free(bufSend);
			count = 0;
		}
	}	/* ---------------- U-V-W Velocity: front & back COMMUNICATION ENDS ------------- */
}
	

