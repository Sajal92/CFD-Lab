#include "helper.h"
#include "visual.h"
#include "init.h"
#include "stdio.h"
#include "sor.h"
#include "uvp.h"
#include "boundary_val.h"
#include "string.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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
			double zlength;		  /* length of the domain z-dir.*/
                    	int  imax;                /* number of cells x-direction*/
                    	int  jmax;                /* number of cells y-direction*/
			int kmax;		  /* number of cells z-direction*/
                    	double dx;                /* length of a cell x-dir. */
                    	double dy;                /* length of a cell y-dir. */
			double dz;		  /* length of a cell z-dir. */
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
			double WI;		  /* velocity z-direction */ 
                   	double PI;                /* pressure */
                    	double GX;                /* gravitation x-direction */
                    	double GY;		  /* gravitation y-direction */
			double GZ;		  /* gravitation z-direction */
			double velIN;		  /* velocity at the inflow */
			double velMW[3]; 	  /*the moving wall velocity is a vector*/

			//char *problemOutput = malloc(strlen(szFileName) + 10);	/* Total length of output filename*/
			//char buffer[4] = {0};		/* Buffer name to include the thread rank and timestamp */
			
		  	char *szFileName=args[1]; /* name of the file */
			char *filename = malloc(strlen(szFileName) + 5);
			strcpy(filename, szFileName);
	                strcat(filename, ".dat");
			char *problemGeometry = malloc(strlen(szFileName) + 60);
			
			/* Append name for output vtk file */
			/*strcpy(problemOutput, szFileName);
			strcat(problemOutput, "_output.");
			sprintf(buffer, "%i", myrank);
			strcat(problemOutput, buffer);*/
			

		        //Reading the program configuration file using read_parameters()
	 		int param = read_parameters(filename, &Re, &UI, &VI, &WI, &PI, &GX, &GY, &GZ, &t_end, &xlength, &ylength, &zlength, 								&dt, &dx, &dy, &dz, &imax,&jmax, &kmax, &alpha, &omg, &tau, &itermax, &eps, 							&dt_value, problemGeometry, &velIN, &velMW[0]);
		        
			param++; // Just using param so that C does not throw an error message
			
			//Initializing .vtk folder
			char sol_folder[80];
			sprintf( sol_folder,"Results_%s",szFileName);

		    		mkdir(sol_folder, 0700);


			char sol_directory[80];
			sprintf( sol_directory,"Results_%s/sol", szFileName);
		   	double ***U = matrix2(0,imax+1,0,jmax+1,0,kmax+1);//Dynamically allocating memmory for matrices U,V,W,P,RS,F,G and H
		  	double ***V = matrix2(0,imax+1,0,jmax+1,0,kmax+1); 
			double ***W = matrix2(0,imax+1,0,jmax+1,0,kmax+1); 
		    	double ***P = matrix2(0,imax+1,0,jmax+1,0,kmax+1);
		 	double ***RS =matrix2(1, imax, 1, jmax, 1, kmax );
		 	double ***F = matrix2(0, imax, 1, jmax, 1, kmax );
		 	double ***G = matrix2(1, imax, 0, jmax, 1, kmax );
		 	double ***H = matrix2(1, imax, 1, jmax, 0, kmax );


			int ***Flag = imatrix2(0,imax+1,0,jmax+1,0,kmax+1);
		 	t=0;
			n=0;
			int n1 = 0;
		
			init_flag(problemGeometry, imax, jmax, kmax, Flag);//Setting flags for all the cells
			
			init_uvwp(UI, VI, WI, PI, Flag,imax, jmax, kmax, U, V, W, P); //Initializing U, V, W and P
			
			printf("Entering while loop \n");
			while (t<t_end)
			{
				//Calculating dt

				calculate_dt(Re, tau, &dt, dx, dy, dz, imax, jmax, kmax, U, V, W); 
				
				//Setting the boundary values for the next time step.

				boundaryvalues(imax, jmax, kmax, U, V, W, P, F, G, H, Flag, velIN, velMW);
				

 				//Determining the values of F,G and H (diffusion and confection).
				calculate_fgh(Re, GX, GY, GZ, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, Flag);

		
				//Calculating the right hand side of the pressure equation.		 
				calculate_rs(dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS,Flag);
				 
		
				it=0;
				res = 1;

				while (it<itermax && res>eps)  
				{

					sor(omg, dx, dy, dz, imax, jmax, kmax, P, RS, &res, Flag);     				
					it=it+1;
				}

				calculate_uvw(dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, P, Flag);

	
				reset_obstacles(U, V,W, P, Flag, imax, jmax,kmax);				

				printf("PROCESS: main.c -> SOR iterations = %d ,residual = %f dt = %f time = %f\n", it-1, res,dt, t);

				if(t>=n1*dt_value){
				write_vtkFile(sol_directory, n, xlength, ylength, zlength, imax, jmax, kmax, dx, dy, dz, U, V, W, P,Flag);
					n1 = n1 + 1;

				}
				
				t=t+dt;
				n=n+1;
				
			}

			
  return -1;
}
