#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "parallel.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */	
int main(int argn, char** args){
			double xlength;           /* length of the domain x-dir.*/
                    	double ylength;           /* length of the domain y-dir.*/
                    	int  imax;                /* number of cells x-direction*/
                    	int  jmax;                /* number of cells y-direction*/
                    	double dx;                /* length of a cell x-dir. */
                    	double dy;                /* length of a cell y-dir. */
                   	double t ;                /* gravitation y-direction */
                    	double t_end;             /* end time */
                    	double dt;                /* time step */
                    	double tau;               /* safety factor for time step*/
                    	double dt_value;          /* time interval for writing visualization data in afile*/
                    	int n;                    /* current time iteration step*/
                    	int itermax;              /* max. number of iterations  */
                    	int it;                   /* SOR iteration counter*/
                    	double res;               /*residual norm of the pressure equation*/
                    	double eps;               /* accuracy bound for pressure*/
                    	double alpha;             /* uppwind differencing factor*/
                    	double omg;               /* relaxation factor */
                    	double Re;                /* reynolds number   */
                    	double UI;                /* velocity x-direction */
                    	double VI;                /* velocity y-direction */
                   	double PI;                /* pressure */
                    	double GX;                /* gravitation x-direction */
                    	double GY;
		  	const char *szFileName="cavity100.dat";       /* name of the file */

			/* MPI variables */
			int iproc = 0;
			int jproc = 0;
			int myrank = 0;
			int il = 0;
			int ir = 0;
			int jb = 0;
			int jt = 0;
			int rank_l = 0;
			int rank_r = 0;
			int rank_b = 0;
			int rank_t = 0;
			int omg_i = 0;
			int omg_j = 0;
			int num_proc = 0;
			int comm_size = 0;
			int x_dim, y_dim;

			/* Initialize MPI and begin parallel processes */
			MPI_Init(&argn, &args);
			MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
			MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

		        //Reading the program configuration file using read_parameters()
			if(myrank==0){
		 		int param = read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
				
				param++; 		// Just using param so that C does not throw an error message

				if(param != 1){
					Programm_Stop("Input parameters potentially corrupt!");
					return -1;
				}
				MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
				if(comm_size != (iproc*jproc)){
					printf("The number of processes entered via 'mpirun -np %u' does not equal the size specified by the input file (iproc * jproc = %u)\n", comm_size, iproc*jproc);
					return -1;
				}
			}

			/* Broadcast parameters to the remaining processes */
			MPI_Bcast(&Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&UI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&VI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&GX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&GY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&imax, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&jmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);


			/* Initialize process dependent variables */
			init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l,
			    &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);

			x_dim = ir - il + 1;
			y_dim = jt - jb + 1;
	
			double **P = matrix(il-1,ir+1,jb-1,jt+1);
			double **U = matrix(il-2,ir+1,jb-1,jt+1); //Dynamically allocating memmory for matrices U,V,P, RS, F and G
			double **V = matrix(il-1,ir+1,jb-2,jt+1); 
		 	double **F = matrix(il-2,ir+1,jb-1,jt+1); 
			double **G = matrix(il-1,ir+1,jb-2,jt+1);    	
		 	double **RS = matrix(il,ir,jb,jt);

			/* Initialize matrices for velocity, pressure, rhs, etc. */
			init_uvp(UI, VI, PI, x_dim, y_dim, U, V, P);

		 	
		 	t=0;
			n=0;
			int n1 = 0;
			int chunk = 0;
			int residue;
			double* bufSend = 0;
			double* bufRecv = 0;
			MPI_Status status;
			//init_uvp(UI,VI,PI,imax,jmax,U,V,P);    //Initializing U, V and P
			
			while (t<t_end)
			{
				
				boundaryvalues(x_dim, y_dim, U, V, rank_l, rank_r, rank_b, rank_t);		 //Setting the boundary values for the next time step.
				calculate_fg(Re,GX,GY,alpha,dt,dx,dy,x_dim,y_dim,U,V,F,G);			 //Determining the values of F and G (diffusion and confection).
				calculate_rs(dt,dx,dy,x_dim,y_dim,F,G,RS);					 //Calculating the right hand side of the pressure equation.
				
				it=0;
				res = 1;
			
				while (it<itermax && res>eps)  //Iterating the pressure poisson equation until the residual becomes smaller than eps or the maximal number of iterations is performed. 
				{
					sor(omg,dx,dy,imax,jmax,P,RS,&res,il,ir,jb,jt,rank_l,rank_r,rank_b,rank_t);      				//Within the iteration loop the operation sor() is used.

 			/* Communicate between processes regarding pressure boundaries */
			pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);
			  /* Sum the squares of all local residuals then square root that sum for global residual */
				
				MPI_Allreduce(&res, &residue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				residue = sqrt((residue)/(imax*jmax));
					it=it+1;
				}
				
				calculate_uv(dt,dx,dy,x_dim,y_dim,U,V,F,G,P);   //Calculating the velocity at the next time step.

				
				uv_comm(U, V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);

				calculate_dt(Re,tau,&dt,dx,dy,x_dim,y_dim,U,V); 

				if(t>=n1*dt_value){
					write_vtkFile("file", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
					n1 = n1 + 1;
				}
				
				t=t+dt;
				n=n+1;
				
			}

			/* Deallocate heap memory */
			free_matrix(U, 0, x_dim+1, 0, y_dim+1);
			free_matrix(V, 0, x_dim+1, 0, y_dim+1);
			free_matrix(P, 0, x_dim+1, 0, y_dim+1);
			free_matrix(F, 0, x_dim+1, 0, y_dim+1);
			free_matrix(G, 0, x_dim+1, 0, y_dim+1);
			free_matrix(RS, 0, x_dim+1, 0, y_dim+1);
			
  return -1;
}
