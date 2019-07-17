#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "string.h"


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
int main(int argc, char* argv[]){
			int select=0;
			char *problem = argv[1];
			char *filename = malloc(strlen(problem) + 5);
			strcpy(filename, problem);
	                strcat(filename, ".dat");

			if(strcmp(problem, "karman_vortex")==0) select=1;
			if(strcmp(problem, "step_flow")==0) select=2;
			if(strcmp(problem, "natural_convection")==0) select=3;
			if(strcmp(problem, "fluid_trap")==0) select=4;
			if(strcmp(problem, "rb_convection")==0) select=5;
			char* geometry = (char*)(malloc(sizeof(char)*6));
			printf("%d", select);


			//define parameter variables
			double Re;                /* reynolds number   */
			double UI;                /* velocity x-direction */
			double VI;                /* velocity y-direction */
			double PI;                /* pressure */
			double GX;                /* gravitation x-direction */
			double GY;                /* gravitation y-direction */
			double t_end;             /* end time */
			double xlength;           /* length of the domain x-dir.*/
			double ylength;           /* length of the domain y-dir.*/
			double dt;                /* time step */
			double dx;                /* length of a cell x-dir. */
			double dy;                /* length of a cell y-dir. */
			int  imax;                /* number of cells x-direction*/
			int  jmax;                /* number of cells y-direction*/
			double alpha;             /* uppwind differencing factor*/
			double omg;               /* relaxation factor */
			double tau;               /* safety factor for time step*/
			int  itermax;             /* max. number of iterations for pressure per time step  */
			double eps;               /* accuracy bound for pressure*/
			double dt_value;           /* time for output */
			double Pr;
			double TI;
			double T_h;
			double T_c;
			double beta;

			//Read and assign the parameter values from file
			read_parameters(filename, &imax, &jmax, &xlength, &ylength, &dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax,&GX, &GY, &Re, &Pr, &UI, &VI, &PI, &TI, &T_h, &T_c,  
					&beta, &dx, &dy, geometry);

			//temperature =1 => Flag for heat transfer problems to include temperature equations for solving
			int temperature = 1;
			if(((select==1)||(select==2))){
				if( (Pr!=0)||(TI!=0)||(T_h!=0)||(T_c!=0)||(beta!=0) ){
					char szBuff[80];
					sprintf( szBuff, "Input file incompatible. Please check .dat file. \n");
					ERROR( szBuff );
				}
				else  temperature = 0;
			}


			//Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
			printf("Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap... \n");
			double **P = matrix(0, imax-1, 0, jmax-1);
			double **U = matrix(0, imax-1, 0, jmax-1);
			double **V = matrix(0, imax-1, 0, jmax-1);
			double **F = matrix(0, imax-1, 0, jmax-1);
			double **G = matrix(0, imax-1, 0, jmax-1);
			double **RS = matrix(0, imax-1, 0, jmax-1);
			int **flag = imatrix(0, imax-1, 0, jmax-1);
			double **T;
			if(temperature){	
				T = matrix(0, imax-1, 0, jmax-1);
			}
			printf("Matrices allocated... \n \n");

			//Flag Initialization
			init_flag(geometry, imax, jmax, flag);

			//Initialize the U, V and P
			init_uvp(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag, temperature);


			//Make solution folder
			struct stat st = {0};
			char sol_folder[80];
			sprintf( sol_folder,"Results_%s",problem);
			if (stat(sol_folder, &st) == -1) {
		    		mkdir(sol_folder, 0700);
			}

			char sol_directory[80];
			sprintf( sol_directory,"Results_%s/sol", problem);
			
			// Algorithm starts from here
			printf("Alogrithm started........\n");
			double t=0; int n=0; int n1=0;
			    
			while (t < t_end) {
	
				calculate_dt(Re,tau,&dt,dx,dy,imax,jmax, U, V, Pr, temperature);
				printf("t = %f ,dt = %f, ",t,dt);
							
				boundaryvalues(imax, jmax, U, V, flag);

				if(temperature){
					calculate_temp(T, Pr, Re, imax, jmax, dx, dy, dt, alpha, U, V, flag, TI, T_h, T_c, select);
				}

				spec_boundary_val(imax, jmax, U, V, flag);

				calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G,flag, beta, T, temperature);
									
				calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);
											
				int it = 0;
				double res = 10.0;

				do{
					sor(omg,dx,dy,imax,jmax,P,RS,&res,flag);
						++it;
				} while(it<itermax && res>eps);


				printf("SOR iterations = %d ,residual = %f \n", it-1, res);

				if((it==itermax)&&(res>eps)){
					printf("WARNING: SOR Iteration limit reached before convergence. \n");
				}


				calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);


				reset_obstacles(U, V, P, T, flag, imax, jmax,temperature);
	

				if ((t >= n1*dt_value)&&(t!=0.0)){
					write_vtkFile(sol_directory ,n ,xlength ,ylength ,imax-2 ,jmax-2,dx ,dy ,U ,V ,P,T,temperature);
					printf("Writing Solutions at %f seconds in the file \n",n1*dt_value);
				    	n1=n1+ 1;
				    	continue;
				}
			t =t+ dt;
			n = n+ 1;
			}


			printf("Algorithm successfully executed...\n \n");
			printf("Freeing dynamically allocated memory...\n");
			    //Free memory
			free_matrix( P, 0, imax-1, 0, jmax-1);
			free_matrix( U, 0, imax-1, 0, jmax-1);
			free_matrix( V, 0, imax-1, 0, jmax-1);
			free_matrix( F, 0, imax-1, 0, jmax-1);
			free_matrix( G, 0, imax-1, 0, jmax-1);
			free_matrix(RS, 0, imax-1, 0, jmax-1);
			free_imatrix(flag, 0, imax-1, 0, jmax-1);
			if(temperature) { 
				free_matrix(T, 0, imax-1, 0, jmax-1);
			}
			free(geometry);
			free(problem);

			printf("End \n");
			return -1;
    
}

