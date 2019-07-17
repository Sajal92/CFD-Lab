#include <math.h>
#include "uvp.h"
#include"helper.h"
#include"boundary_val.h"
#include <stdio.h>

void calculate_dt(
		double Re,
		double tau,
		double *dt,
		double dx,
		double dy,
		double dz,
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W
)
{
	
	double temp1,temp2,temp3,temp4,temp5,temp6;
    	double U1=fabs(U[0][0][0]);
	double V1=fabs(V[0][0][0]);
	double W1=fabs(W[0][0][0]);
        for(int c=0 ; c <=imax+1; c++ ){
  	      for(int d= 0 ; d <=jmax+1 ; d++ ){
		for(int e = 0;e<=kmax+1; e++){
        	 if ( fabs(U[c][d][e]) > fabs(U1) ){
          		  U1= U[c][d][e];  //Calculating Umax
      		 }
		 if ( fabs(V[c][d][e]) > fabs(V1) ){
          		  V1= V[c][d][e];	//Calculating Vmax						    				
                 }
		 if ( fabs(W[c][d][e]) > fabs(W1) ){
          		  W1= V[c][d][e];	//Calculating Wmax						    				
                 }
  	       }
	}
	}
	temp1=0.5*Re*(1.0/(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz))); 
        temp2=(dx)/fabs(U1);
        temp3=(dy)/fabs(V1);
	temp4=(dz)/fabs(W1);
        temp5=fmin(temp1,temp2);
	temp6=fmin(temp3,temp4);
        *dt=tau*fmin(temp5,temp6);							       // Calculating dt with tau belongs to [0,1]
    
}

void calculate_fgh(
		double Re,
		double GX,
		double GY,
		double GZ,
		double alpha,
		double dt,
		double dx,
		double dy,
		double dz,
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***F,
		double ***G,
		double ***H,
		int ***flag
){



	double uijk, vijk, wijk, uujk, udjk, uiuk, uidk, uiju, uijd, viuk, vidk, viju, vijd, vujk, vdjk, wiju, wijd, wujk, wdjk, wiuk, widk, uduk, udju, vudk, vidu, wiud,wujd;
	double Dx = 1/dx, Dy = 1/dy, Dz = 1/dz;
	double ReInv = 1.0/Re;
	int i, j, k;
	for (i=1; i<imax+1; i++){
		for  (j=1; j<jmax+1; j++){
			for  (k=1; k<kmax+1; k++){
				if ((flag[i][j][k] & (1<<0)))
					{
					uijk = U[i][j][k];   vijk = V[i][j][k];   wijk = W[i][j][k];
					uujk = U[i+1][j][k]; viuk = V[i][j+1][k]; wiju = W[i][j][k+1];
					udjk = U[i-1][j][k]; vidk = V[i][j-1][k]; wijd = W[i][j][k-1];
					uiuk = U[i][j+1][k]; viju = V[i][j][k+1]; wujk = W[i+1][j][k];
					uidk = U[i][j-1][k]; vijd = V[i][j][k-1]; wdjk = W[i-1][j][k];
					uiju = U[i][j][k+1]; vujk = V[i+1][j][k]; wiuk = W[i][j+1][k];
					uijd = U[i][j][k-1]; vdjk = V[i-1][j][k]; widk = W[i][j-1][k];
					uduk = U[i-1][j+1][k]; udju = U[i-1][j][k+1];
					vudk = V[i+1][j-1][k]; vidu = V[i][j-1][k+1];
					wiud = W[i][j+1][k-1]; wujd = W[i+1][j][k-1];
					if ((flag[i+1][j][k] & (1<<0))){ //calculation is only on edges separating two fluid cells.
						F[i][j][k] = uijk + dt*(ReInv*((uujk-2*uijk+udjk)*Dx*Dx + (uiuk-2*uijk+uidk)*Dy*Dy + (uiju-2*uijk+uijd)*Dz*Dz) -
								0.25*Dx*(pow((uijk+uujk),2)-pow((udjk+uijk),2) + alpha*(fabs(uijk+uujk)*(uijk-uujk)-fabs(udjk+uijk)*(udjk-uijk))) -
								0.25*Dy*((vijk+vujk)*(uijk+uiuk)-(vidk+vudk)*(uidk+uijk) + alpha*(fabs(vijk+vujk)*(uijk-uiuk)-fabs(vidk+vudk)*(uidk-uijk))) -
								0.25*Dz*((wijk+wujk)*(uijk+uiju)-(wijd+wujd)*(uijd+uijk) + alpha*(fabs(wijk+wujk)*(uijk-uiju)-fabs(wijd+wujd)*(uijd-uijk))) +
								GX);
					}
					if ((flag[i][j+1][k] & (1<<0))){
						G[i][j][k] = vijk + dt*(ReInv*((vujk-2*vijk+vdjk)*Dx*Dx + (viuk-2*vijk+vidk)*Dy*Dy + (viju-2*vijk+vijd)*Dz*Dz) -
								0.25*Dy*(pow((vijk+viuk),2)-pow((vidk+vijk),2) + alpha*(fabs(vijk+viuk)*(vijk-viuk)-fabs(vidk+vijk)*(vidk-vijk))) -
								0.25*Dx*((uijk+uiuk)*(vijk+vujk)-(udjk+uduk)*(vdjk+vijk) + alpha*(fabs(uijk+uiuk)*(vijk-vujk)-fabs(udjk+uduk)*(vdjk-vijk))) -
								0.25*Dz*((wijk+wiuk)*(vijk+viju)-(wijd+wiud)*(vijd+vijk) + alpha*(fabs(wijk+wiuk)*(vijk-viju)-fabs(wijd+wiud)*(vijd-vijk))) +
								GY);
					}
					if ((flag[i][j][k+1] & (1<<0))){
						H[i][j][k] = wijk + dt*(ReInv*((wujk-2*wijk+wdjk)*Dx*Dx + (wiuk-2*wijk+widk)*Dy*Dy + (wiju-2*wijk+wijd)*Dz*Dz) -
								0.25*Dz*(pow((wijk+wiju),2)-pow((wijd+wijk),2) + alpha*(fabs(wijk+wiju)*(wijk-wiju)-fabs(wijd+wijk)*(wijd-wijk))) -
								0.25*Dx*((uijk+uiju)*(wijk+wujk)-(udjk+udju)*(wdjk+wijk) + alpha*(fabs(uijk+uiju)*(wijk-wujk)-fabs(udjk+udju)*(wdjk-wijk))) -
								0.25*Dy*((vijk+viju)*(wijk+wiuk)-(vidk+vidu)*(widk+wijk) + alpha*(fabs(vijk+viju)*(wijk-wiuk)-fabs(vidk+vidu)*(widk-wijk))) +
								GZ);
					}
					
     			}

	else {
			
					if (B_F(flag[i][j][k])){ //U neighb. is fluid
						H[i][j][k] = W[i][j][k]; 
						//H[i][j][k+1]=W[i][j][k+1]
					} else if (B_B(flag[i][j][k])){  //D neighb. is fluid
						H[i][j][k-1] = W[i][j][k-1]; 
						
					}
					if (B_O(flag[i][j][k])){ // O neigh. is fluid
						F[i][j][k] = U[i][j][k]; 
					} else if (B_W(flag[i][j][k])) { // W neighb. is fluid
						F[i-1][j][k] = U[i-1][j][k]; 
					}
					if ( (B_N(flag[i][j][k] ))) { //N neigh. is fluid
						G[i][j][k] = V[i][j][k]; 
					} else if ( B_S(flag[i][j][k] )) { //S neigh. is fluid
						G[i][j-1][k] = V[i][j-1][k]; 
					}
					}
				}
			H[i][j][0]    = W[i][j][0];
			H[i][j][kmax] = W[i][j][kmax];
		}
		for (k=1; k<kmax+1; k++){
			G[i][0][k]   = V[i][0][k];
			G[i][jmax][k] = V[i][jmax][k];
		}		
	}
	for (j=1; j<jmax+1; j++){
		for (k=1; k<kmax+1; k++){
			F[0][j][k] = U[0][j][k];
			F[imax][j][k] = U[imax][j][k];
		}
	}
}
void calculate_uvw(
		double dt,
		double dx,
		double dy,
		double dz,
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***F,
		double ***G,
		double ***H,
		double ***P,
		int ***flag
){
	//updating the U components of velocity
	for(int i=1; i<imax; i++){
		for(int j=1; j<jmax+1; j++){
			for(int k=1; k<kmax+1; k++){
				if( ((flag[i][j][k]&(1<<0))&flag[i+1][j][k]) || ( (flag[i+1][j][k] & (1<<3)) && (flag[i][j][k]&(1<<0))) ){
					U[i][j][k] = F[i][j][k] - ((dt/dx) * (P[i+1][j][k] - P[i][j][k]));
				}
			}
		}
	}
	//updating the V components of velocity
	for(int i=1; i<imax+1; i++){
		for(int j=1; j<jmax; j++){
			for(int k=1; k<kmax+1; k++){
				if( ((flag[i][j][k]&(1<<0))&flag[i][j+1][k]) || ( (flag[i][j+1][k] & (1<<3)) && (flag[i][j][k]&(1<<0))) ){
					V[i][j][k] = G[i][j][k] - ((dt/dy) * (P[i][j+1][k] - P[i][j][k]));
				}
			}
		}
	}
	//updating the W components of velocity
	for(int i=1; i<imax+1; i++){
		for(int j=1; j<jmax+1; j++){
			for(int k=1; k<kmax; k++){
				if( ((flag[i][j][k]&(1<<0))&flag[i][j][k+1]) || ( (flag[i][j][k+1] & (1<<3)) && (flag[i][j][k]&(1<<0))) ){
					W[i][j][k] = H[i][j][k] - ((dt/dz) * (P[i][j][k+1] - P[i][j][k]));
				}
			}
		}
	}

}

void calculate_rs(
		double dt,
		double dx,
		double dy,
		double dz,
		int imax,
		int jmax,
		int kmax,
		double ***F,
		double ***G,
		double ***H,
		double ***RS,
		int ***flag
)
{	
	for (int i=1; i<imax+1; i++){
		for (int j=1; j<jmax+1; j++){
			for(int k=1; k<kmax+1; k++){
			   //Updating the right hand side of pressure poisson equation
		
			   RS[i][j][k] = (1/(dt*dx))*(F[i][j][k] - F[i-1][j][k])+(1/(dt*dy))*(G[i][j][k] - G[i][j-1][k])+(1/(dt*dz))*(H[i][j][k] - H[i][j][k-1]); 
				
		    }	
		}
	}
	
}

void reset_obstacles(double ***U, double ***V, double ***W, double ***P, int ***flag, int imax, int jmax, int kmax){
	for (int i = 0; i<=imax+1;i++){
		for (int j=0;j<=jmax+1; j++){
			for (int k=0;k<=kmax+1; k++){
			if(flag[i][j][k]&( (1<<1)|(1<<2)) ){
					U[i][j][k] = 0;
					V[i][j][k] = 0;
					W[i][j][k] = 0;
					P[i][j][k] = 0;
				}
			}	

		}	
	}
}
