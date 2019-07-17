#include "sor.h"
#include <math.h>
#include"boundary_val.h"
#include <stdio.h>

void sor(
		double omg,
		double dx,
		double dy,
		double dz,
		int    imax,
		int    jmax,
		int    kmax,
		double ***P,
		double ***RS,
		double *res,
		int ***flag
) {
	double rloc;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));

	double FluidCells = 0.0;
	double tmp = 0;

 //set pressure at boundary- no-slip/free-slip/inflow/outflow
for(int i = 0; i <= imax+1; i++) {
	
			for(int j = 0; j <= jmax+1; j++) {
				for(int k = 0; k <= kmax+1; k++) {
				   if((flag[i][j][k] & ((1<<1)|(1<<4)|(1<<2)|(1<<12)))) {	
					if( B_O(flag[i][j][k])) { P[i][j][k]  = P[i+1][j][k]; }
					if( B_W(flag[i][j][k])) { P[i][j][k]  = P[i-1][j][k]; }

					if( B_N(flag[i][j][k])) { P[i][j][k]  = P[i][j+1][k]; }
					if( B_S(flag[i][j][k])) { P[i][j][k]  = P[i][j-1][k]; }

					if( B_F(flag[i][j][k])) { P[i][j][k]  = P[i][j][k+1]; }
					if( B_B(flag[i][j][k])) { P[i][j][k]  = P[i][j][k-1]; }

					if( B_NO(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j+1][k])*0.5; }
					if( B_NW(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j+1][k])*0.5; }
					if( B_SO(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j-1][k])*0.5; }
					if( B_SW(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j-1][k])*0.5; }
					if( B_NF(flag[i][j][k])) { P[i][j][k] = (P[i][j][k+1] + P[i][j+1][k])*0.5; }
					if( B_NB(flag[i][j][k])) { P[i][j][k] = (P[i][j][k-1] + P[i][j+1][k])*0.5; }
					if( B_SF(flag[i][j][k])) { P[i][j][k] = (P[i][j][k+1] + P[i][j-1][k])*0.5; }
					if( B_SB(flag[i][j][k])) { P[i][j][k] = (P[i][j][k-1] + P[i][j-1][k])*0.5; }
					if( B_OF(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j][k+1])*0.5; }
					if( B_WF(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j][k+1])*0.5; }
					if( B_OB(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j][k-1])*0.5; }
					if( B_WB(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j][k-1])*0.5; }

					if( B_NOF(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j+1][k] + P[i][j][k+1])*1.0/3.0; }
					if( B_NWF(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j+1][k] + P[i][j][k+1])*1.0/3.0; }
					if( B_SOF(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j-1][k] + P[i][j][k+1])*1.0/3.0; }
					if( B_SWF(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j-1][k] + P[i][j][k+1])*1.0/3.0; }
					if( B_NOB(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j+1][k] + P[i][j][k-1])*1.0/3.0; }
					if( B_NWB(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j+1][k] + P[i][j][k-1])*1.0/3.0; }
					if( B_SOB(flag[i][j][k])) { P[i][j][k] = (P[i+1][j][k] + P[i][j-1][k] + P[i][j][k-1])*1.0/3.0; }
					if( B_SWB(flag[i][j][k])) { P[i][j][k] = (P[i-1][j][k] + P[i][j-1][k] + P[i][j][k-1])*1.0/3.0; }

					}

				   if((flag[i][j][k] & ((1<<3)))) {	P[i][j][k] = 0;  }
					
					}
				}
			}
			


  /* SOR iteration */
	for(int i = 1; i <= imax; i++) {
		for(int j = 1; j<=jmax; j++) {
			for(int k = 1; k<=kmax; k++) {
				if((flag[i][j][k] & (1<<0))){


					P[i][j][k] = (1.0-omg)*P[i][j][k]+
					                coeff*(
					                      (P[i+1][j  ][k  ]+P[i-1][j  ][k  ])/(dx*dx) +
					                      (P[i  ][j+1][k  ]+P[i  ][j-1][k  ])/(dy*dy) +
					                      (P[i  ][j  ][k+1]+P[i  ][j  ][k-1])/(dz*dz)

					                    - RS[i  ][j  ][k  ]
					                       );


					FluidCells++;
				}
			}
		}
	}


	/* compute the residual */
	rloc = 0;
	for(int i = 1; i <= imax; i++) {
		for(int j = 1; j <= jmax; j++) {
			for(int k = 1; k <= kmax; k++){
				if((flag[i][j][k] & (1<<0))){
					tmp =  (P[i+1][j][k]-2.0*P[i][j][k]+P[i-1][j][k])/(dx*dx) + (P[i][j+1][k]-2.0*P[i][j][k]+P[i][j-1][k])/(dy*dy) + (P[i][j][k+1]-2.0*P[i][j][k]+P[i][j][k-1])/(dz*dz) - RS[i][j][k];
					rloc += tmp*tmp;

				}
			}
		}
	}
	rloc = rloc/FluidCells;
	rloc = sqrt(rloc);
	/* set residual */
	*res = rloc;

}
