#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void Program_Message (char *txt);
/* produces a stderr text output  */


void Programm_Sync (char *txt);
/* produces a stderr textoutput and synchronize all processes */


void Programm_Stop (char *txt);
/* all processes will produce a text output, be synchronized and finished */


void init_parallel (int iproc,int jproc,int kproc,int imax,int jmax,int kmax,
					int *myrank,int *il,int *ir,int *jb,int *jt,int *kf, int *kt,
					int *rank_l,int *rank_r,int *rank_b,int *rank_t, int *rank_f, int *rank_kb,
					int *omg_i,int *omg_j,int *omg_k,int num_proc);

void flag_comm(	double ***flag,
					int il, int ir, 
					int jb, int jt,
					int kf, int kb, 
					int rank_l, int rank_r,
					int rank_b, int rank_t,
					int rank_f, int rank_kb, 
					double *bufSend,double *bufRecv,
					MPI_Status *status, int chunk);

void pressure_comm(	double ***P,
					int il, int ir, 
					int jb, int jt,
					int kf, int kb, 
					int rank_l, int rank_r,
					int rank_b, int rank_t,
					int rank_f, int rank_kb, 
					double *bufSend,double *bufRecv,
					MPI_Status *status, int chunk);
/* this method exchanges pressure values between processes that treat
 adjacent sub-domains*/
/* Called within SOR before actually calculating the residual */


void uv_comm (	double ***U, double ***V, double ***W,
				int il,int ir,
				int jb, int jt, 
				int kf, int kb,
				int rank_l, int rank_r,
				int rank_b, int rank_t,
				int rank_f, int rank_kb,
				double *bufSend, double *bufRecv, 
				MPI_Status *status, int chunk);
/* this method exchanges the velocity values U und V between the processes
treating adjacent subdomains*/
/* Called at the end of calculate_uv */
