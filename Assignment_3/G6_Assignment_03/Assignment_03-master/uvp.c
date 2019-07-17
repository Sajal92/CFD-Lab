#include<stdio.h>
#include "uvp.h"
#include "math.h"
#include <mpi.h>

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{
		double a, b,c,d,du2x2,du2y2,du2dx,duvy,dv2y2,dv2x2,dv2dy,duvx;
		for (int j=1;j<=jmax;j++){
			F[0][j]=U[0][j];                                                      //Boundary conditions for F
			F[imax][j]=U[imax][j];
		}

		for(int i=1;i<=imax;i++){
			G[i][0]=V[i][0];                                                      // Boundary conditions for G
			G[i][jmax]=V[i][jmax];
		}
		
		for(int i=1; i<imax; i++){
			for (int j=1; j<=jmax; j++){
				du2x2 = (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);              //Central difference scheme for second derivatives
				du2y2 = (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);	      

				a = (U[i][j]+U[i+1][j])/2;
				b = (U[i-1][j]+U[i][j])/2;

				du2dx = (a*a-b*b)/dx+	(alpha/dx)*(fabs(a)*((U[i][j]-U[i+1][j])/2) - fabs(b)*((U[i-1][j]-U[i][j])/2));//Modified Donor Cell method for (d(u^2)/dx)

				duvy = (1/(4*dy))*((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1]) - (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))	+	(alpha/(4*dy))*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1]) - fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]));                                //Modified Donor Cell method on convective term (d(uv)/dy)

				F[i][j] = U[i][j] + dt*((1/Re)*(du2x2+du2y2) - du2dx - duvy + GX);
			}
		}

		for(int i=1; i<=imax; i++){
			for (int j=1; j<jmax; j++){
				dv2y2 = (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);              //Central difference scheme for second derivatives   
				dv2x2 = (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);	     

				c = (V[i][j]+V[i][j+1])/2;
				d = (V[i][j-1]+V[i][j])/2;

				dv2dy = (c*c-d*d)/dy	+	(alpha/dy)*(fabs(c)*((V[i][j]-V[i][j+1])/2) - fabs(d)*((V[i][j-1]-V[i][j])/2));//Modified Donor cell method for d(v^2)/dy  
				duvx = (1/(4*dx))*((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j]) - (U[i-1][j] + U[i-1][j+1])*(V[i-1][j]+V[i][j]))		+	(alpha/(4*dx))*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j]) - fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]));                      //Modified Donor cell method for d(uv)/dx

				G[i][j] = V[i][j] + dt*((1/Re)*(dv2x2+dv2y2) - duvx - dv2dy + GY);
			}
		}
	
}



void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{	
	for (int i=1; i<=imax; i++){
		for (int j=1; j<=jmax; j++){
			RS[i][j] = (1/(dt*dx))*(F[i][j] - F[i-1][j])+(1/(dt*dy))*(G[i][j] - G[i][j-1]);//Updating the right hand side of pressure poisson equation
		}
	}
	
}



void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
	
	double temp1,temp2,temp3,temp4;
    	double U1=fabs(U[0][0]);
	double V1=fabs(V[0][0]);
	double U_max, V_max;
	
        for(int c=0 ; c <=imax ; c++ ){
  	      for(int d= 0 ; d <=jmax ; d++ ){
        	 if ( fabs(U[c][d]) > fabs(U1) ){
          		  U1= U[c][d];							       //Calculating Umax
      		 }
		 if ( fabs(V[c][d]) > fabs(V1) ){
          		  V1= V[c][d];							    //Calculating Vmax						    				
                 }
  	       }
	}

	MPI_Allreduce(&U1, &U_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&V1, &V_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

	temp1=((Re/2)*(dx*dx)*(dy*dy))/((dx*dx)+(dy*dy)); 
        temp2=(dx)/fabs(U_max);
        temp3=(dy)/fabs(V_max);
        temp4=fmin(temp1,temp2);
        *dt=tau*fmin(temp4,temp3);							       // Calculating dt with tau belongs to [0,1]
    
}




void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
	// we may have to change limits of i and j
	for (int i=1; i<=imax-1; i++){
		for (int j=1; j<=jmax; j++){
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j] - P[i][j]);			//Updating the U components of velocity
		}
	}
	
	for (int i=1; i<=imax; i++){
		for (int j=1; j<=jmax-1; j++){
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j+1] - P[i][j]);			//updating the V components of velocity
		}
	}
}



