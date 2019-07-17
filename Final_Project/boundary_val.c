#include "boundary_val.h"
#include <stdio.h>

void spec_boundary_val(int imax,int jmax, int kmax, double ***U,double ***V, double***W,int ***flag, double velIN)
{
	for(int i = 0; i<=imax+1; ++i){
		for(int j = 0; j<=jmax+1; ++j){
			for(int k=0; k<=kmax+1; ++k){
			if(flag[i][j][k] &(1<<4)){
			
			if (B_O(flag[i][j][k])) {	U[i][j][k] = velIN; }
			if (B_W(flag[i][j][k])) { 	U[i-1][j][k] = -velIN; }
			if (B_N(flag[i][j][k])) { 	V[i][j][k] = velIN;	}
			if (B_S(flag[i][j][k])) {	V[i][j-1][k] = -velIN;printf("CHecking inflow \n\n\n\n"); }
			if (B_F(flag[i][j][k])) { 	W[i][j][k] = velIN; }
			if (B_B(flag[i][j][k])) { 	W[i][j][k-1] = -velIN; }
			}
		}
	}

}
}
	
int B_O(int flag){
	return ( (flag & (1<<8)) && ~( flag & ((1<<5) | (1<<6) | (1<<7) | (1<<9) | (1<<10)) ) );
}

int B_W(int flag){
	return ( (flag & (1<<7)) && ~( flag & ((1<<5) | (1<<6) | (1<<8) | (1<<9) | (1<<10)) ) );
}

int B_N(int flag){
	return ( (flag & (1<<5)) && ~( flag & ((1<<8) | (1<<6) | (1<<7) | (1<<9) | (1<<10)) ) );
}

int B_S(int flag){
	return ( (flag & (1<<6)) && ~( flag & ((1<<5) | (1<<8) | (1<<7) | (1<<9) | (1<<10)) ) );
}

int B_F(int flag){
	return ( (flag & (1<<9)) && ~( flag & ((1<<5) | (1<<6) | (1<<7) | (1<<8) | (1<<10)) ) );
}

int B_B(int flag){
	return ( (flag & (1<<10)) && ~( flag & ((1<<5) | (1<<6) | (1<<7) | (1<<9) | (1<<8)) ) );
}

int B_NO(int flag){
	return ( (flag & (1<<8)) &&  (flag & (1<<5)) && ~( flag & ((1<<6) | (1<<7) | (1<<9) | (1<<10)) ) );
}

int B_NW(int flag){
	return ( (flag & (1<<7)) &&  (flag & (1<<5)) && ~( flag & ((1<<6) | (1<<8) | (1<<9) | (1<<10)) ) );
}

int B_SO(int flag){
	return ( (flag & (1<<8)) &&  (flag & (1<<6)) && ~( flag & ((1<<5) | (1<<7) | (1<<9) | (1<<10)) ) ); 
}

int B_SW(int flag){
	return ( (flag & (1<<7)) &&  (flag & (1<<6)) && ~( flag & ((1<<5) | (1<<8) | (1<<9) | (1<<10)) ) );
}

int B_NF(int flag){
	return ( (flag & (1<<5)) &&  (flag & (1<<9)) && ~( flag & ((1<<6) | (1<<8) | (1<<7) | (1<<10)) ) );
}

int B_NB(int flag){
	return ( (flag & (1<<5)) &&  (flag & (1<<10)) && ~( flag & ((1<<6) | (1<<8) | (1<<7) | (1<<9)) ) );
}

int B_SF(int flag){
	return ( (flag & (1<<6)) &&  (flag & (1<<9)) && ~( flag & ((1<<5) | (1<<8) | (1<<7) | (1<<10)) ) );
}

int B_SB(int flag){
	return ( (flag & (1<<6)) &&  (flag & (1<<10)) && ~( flag & ((1<<5) | (1<<8) | (1<<7) | (1<<9)) ) );
}

int B_OF(int flag){
	return ( (flag & (1<<8)) &&  (flag & (1<<9)) && ~( flag & ((1<<6) | (1<<5) | (1<<7) | (1<<10)) ) );
}

int B_WF(int flag){
	return ( (flag & (1<<7)) &&  (flag & (1<<9)) && ~( flag & ((1<<6) | (1<<8) | (1<<5) | (1<<10)) ) );
}

int B_OB(int flag){
	return ( (flag & (1<<8)) &&  (flag & (1<<10)) && ~( flag & ((1<<6) | (1<<5) | (1<<7) | (1<<9)) ) );
}

int B_WB(int flag){
	return ( (flag & (1<<7)) &&  (flag & (1<<10)) && ~( flag & ((1<<6) | (1<<8) | (1<<5) | (1<<9)) ) );
}

int B_NOF(int flag){
	return ( (flag & (1<<5)) &&  (flag & (1<<8)) && (flag & (1<<9)) && ~( flag & ((1<<6) | (1<<7) | (1<<10) ) ) );
}

int B_NOB(int flag){
	return ( (flag & (1<<5)) &&  (flag & (1<<8)) && (flag & (1<<10)) && ~( flag & ((1<<6) | (1<<7) | (1<<9) ) ) );
}

int B_NWF(int flag){
	return ( (flag & (1<<5)) &&  (flag & (1<<7)) && (flag & (1<<9)) && ~( flag & ((1<<6) | (1<<8) | (1<<10) ) ) );
}

int B_NWB(int flag){
	return ( (flag & (1<<5)) &&  (flag & (1<<7)) && (flag & (1<<10)) && ~( flag & ((1<<6) | (1<<8) | (1<<9) ) ) );
}

int B_SOF(int flag){
	return ( (flag & (1<<6)) &&  (flag & (1<<8)) && (flag & (1<<9)) && ~( flag & ((1<<5) | (1<<7) | (1<<10) ) ) );
}

int B_SOB(int flag){
	return ( (flag & (1<<6)) &&  (flag & (1<<8)) && (flag & (1<<10)) && ~( flag & ((1<<5) | (1<<7) | (1<<9) ) ) );
}

int B_SWF(int flag){
	return ( (flag & (1<<6)) &&  (flag & (1<<7)) && (flag & (1<<9)) && ~( flag & ((1<<5) | (1<<8) | (1<<10) ) ) );
}

int B_SWB(int flag){
	return ( (flag & (1<<6)) &&  (flag & (1<<7)) && (flag & (1<<10)) && ~( flag & ((1<<5) | (1<<8) | (1<<9) ) ) );
}




void boundaryvalues(
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***P,
		double ***F,
		double ***G,
		double ***H,  
				int ***flag,
				double velIN,
				double *velMW
)
{
	for(int i = 1; i<=imax+1; ++i){
		for(int j = 1; j<=jmax+1; ++j){
		   for (int k = 1; k<=kmax+1; ++k){

			switch(flag[i][j][k]& ((1<<1)|(1<<2)|(1<<3)|(1<<12)|(1<<4)) ) {

				case 1<<1: //No slip conditions

					if ( B_O(flag[i][j][k]) ){

						U[i  ][j  ][k  ]   = 0.0;
						V[i  ][j-1][k  ] = -V[i+1][j-1][k ];
						V[i  ][j  ][k  ] = -V[i+1][j  ][k ];
						W[i  ][j  ][k-1] = -W[i+1][j  ][k-1];
						W[i  ][j  ][k  ] = -W[i+1][j  ][k ];
					}

					if ( B_W(flag[i][j][k]) ) {

						U[i-1][j  ][k  ] = 0.0;
						V[i  ][j-1][k  ] = -V[i-1][j-1][k  ];
						V[i  ][j  ][k  ] = -V[i-1][j  ][k  ];
						W[i  ][j  ][k-1] = -W[i-1][j  ][k-1];
						W[i  ][j  ][k  ] = -W[i-1][j  ][k  ];
						}

					if ( B_N(flag[i][j][k]) ) {

						V[i  ][j  ][k  ] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j+1][k  ];
						U[i  ][j  ][k  ] = -U[i  ][j+1][k  ];
						W[i  ][j  ][k-1] = -W[i  ][j+1][k-1];
						W[i  ][j  ][k  ] = -W[i  ][j+1][k  ];	
					}

					if ( B_S(flag[i][j][k]) ){
						V[i  ][j-1][k  ] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j-1][k  ];
						U[i  ][j  ][k  ] = -U[i  ][j-1][k  ];
						W[i  ][j  ][k-1] = -W[i  ][j-1][k-1];
						W[i  ][j  ][k  ] = -W[i  ][j-1][k  ];	
					}

					if ( B_F(flag[i][j][k]) ){
						W[i  ][j  ][k  ] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j  ][k+1];
						U[i  ][j  ][k  ] = -U[i  ][j  ][k+1];
						V[i  ][j-1][k  ] = -V[i  ][j-1][k+1];
						V[i  ][j  ][k  ] = -V[i  ][j  ][k+1];	
					}

					if ( B_B(flag[i][j][k]) ){
						W[i  ][j  ][k-1] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j  ][k-1];
						U[i  ][j  ][k  ] = -U[i  ][j  ][k-1];
						V[i  ][j-1][k  ] = -V[i  ][j-1][k-1];
						V[i  ][j  ][k  ] = -V[i  ][j  ][k-1];	
					}

					if ( B_NO(flag[i][j][k])  ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j+1][k];
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i+1][j-1][k];
						W[i][j][k] = -(W[i][j+1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j+1][k-1] + W[i+1][j][k-1]) * 0.5;
					}

					if ( B_NW(flag[i][j][k])  ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j+1][k];
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i-1][j-1][k];
						W[i][j][k] = -(W[i][j+1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j+1][k-1] + W[i-1][j][k-1]) * 0.5;
					} 

					if ( B_NF(flag[i][j][k])  ){
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i][j-1][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i][j+1][k-1];
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
					} 

					if ( B_NB(flag[i][j][k])  ){
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i][j-1][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i][j+1][k];
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
					} 

					if ( B_SO(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j-1][k];
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i+1][j][k];
						W[i][j][k] = -(W[i][j-1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j-1][k-1] + W[i+1][j][k-1]) * 0.5;
					}

					if ( B_SW(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j-1][k];
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i-1][j][k];
						W[i][j][k] = -(W[i][j-1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j-1][k-1] + W[i-1][j][k-1]) * 0.5;
					}

					if ( B_SF(flag[i][j][k]) ){
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i][j-1][k-1];
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
					}

					if ( B_SB(flag[i][j][k]) ){
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i][j-1][k];
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
					}

					if ( B_OF(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i+1][j][k-1];
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
					}

					if ( B_WF(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i-1][j][k-1];
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
					}

					if ( B_OB(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i+1][j][k];
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
					}

					if ( B_WB(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i-1][j][k];
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
					}

					if ( B_NOF(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i+1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NWF(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i-1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NOB(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i+1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_NWB(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i-1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_SOF(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i+1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if (B_SWF(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i-1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if ( B_SOB(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i+1][j][k] + W[i][j-1][k]) * 0.5;
						}
					if ( B_SWB(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i-1][j][k] + W[i][j-1][k]) * 0.5;
						}

				break;

				case 1<<4: //Inflow
					if ( B_O(flag[i][j][k]) ){
						U[i  ][j  ][k  ] = 0.0;
						V[i  ][j-1][k  ] = -V[i+1][j-1][k ];
						V[i  ][j  ][k  ] = -V[i+1][j  ][k ];
						W[i  ][j  ][k-1] = -W[i+1][j  ][k-1];
						W[i  ][j  ][k  ] = -W[i+1][j  ][k ];
					}

					if ( B_W(flag[i][j][k]) ) {
						U[i-1][j  ][k  ] = 0.0;
						V[i  ][j-1][k  ] = -V[i-1][j-1][k  ];
						V[i  ][j  ][k  ] = -V[i-1][j  ][k  ];
						W[i  ][j  ][k-1] = -W[i-1][j  ][k-1];
						W[i  ][j  ][k  ] = -W[i-1][j  ][k  ];
						}

					if ( B_N(flag[i][j][k]) ) {
						V[i  ][j  ][k  ] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j+1][k  ];
						U[i  ][j  ][k  ] = -U[i  ][j+1][k  ];
						W[i  ][j  ][k-1] = -W[i  ][j+1][k-1];
						W[i  ][j  ][k  ] = -W[i  ][j+1][k  ];	
					}

					if ( B_S(flag[i][j][k]) ){
						V[i  ][j-1][k  ] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j-1][k  ];
						U[i  ][j  ][k  ] = -U[i  ][j-1][k  ];
						W[i  ][j  ][k-1] = -W[i  ][j-1][k-1];
						W[i  ][j  ][k  ] = -W[i  ][j-1][k  ];	
					}

					if ( B_F(flag[i][j][k]) ){
						W[i  ][j  ][k  ] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j  ][k+1];
						U[i  ][j  ][k  ] = -U[i  ][j  ][k+1];
						V[i  ][j-1][k  ] = -V[i  ][j-1][k+1];
						V[i  ][j  ][k  ] = -V[i  ][j  ][k+1];	
					}

					if ( B_B(flag[i][j][k]) ){
						W[i  ][j  ][k-1] = 0.0;
						U[i-1][j  ][k  ] = -U[i-1][j  ][k-1];
						U[i  ][j  ][k  ] = -U[i  ][j  ][k-1];
						V[i  ][j-1][k  ] = -V[i  ][j-1][k-1];
						V[i  ][j  ][k  ] = -V[i  ][j  ][k-1];	
					}

					if ( B_NO(flag[i][j][k])  ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j+1][k];
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i+1][j-1][k];
						W[i][j][k] = -(W[i][j+1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j+1][k-1] + W[i+1][j][k-1]) * 0.5;
					}

					if ( B_NW(flag[i][j][k])  ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j+1][k];
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i-1][j-1][k];
						W[i][j][k] = -(W[i][j+1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j+1][k-1] + W[i-1][j][k-1]) * 0.5;
					} 

					if ( B_NF(flag[i][j][k])  ){
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i][j-1][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i][j+1][k-1];
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
					} 

					if ( B_NB(flag[i][j][k])  ){
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -V[i][j-1][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i][j+1][k];
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
					} 

					if ( B_SO(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j-1][k];
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i+1][j][k];
						W[i][j][k] = -(W[i][j-1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j-1][k-1] + W[i+1][j][k-1]) * 0.5;
					}

					if ( B_SW(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j-1][k];
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i-1][j][k];
						W[i][j][k] = -(W[i][j-1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = -(W[i][j-1][k-1] + W[i-1][j][k-1]) * 0.5;
					}

					if ( B_SF(flag[i][j][k]) ){
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i][j-1][k-1];
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
					}

					if ( B_SB(flag[i][j][k]) ){
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -V[i][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i][j-1][k];
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
					}

					if ( B_OF(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i+1][j][k-1];
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
					}

					if ( B_WF(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -W[i-1][j][k-1];
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
					}

					if ( B_OB(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -U[i-1][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i+1][j][k];
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
					}

					if ( B_WB(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -U[i][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -W[i-1][j][k];
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
					}

					if ( B_NOF(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i+1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NWF(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i-1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NOB(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i+1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_NWB(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = -(V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i-1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_SOF(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i+1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if (B_SWF(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = -(W[i-1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if ( B_SOB(flag[i][j][k]) ){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = -(U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i+1][j][k] + W[i][j-1][k]) * 0.5;
						}
					if ( B_SWB(flag[i][j][k]) ){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = -(U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = -(V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = -(W[i-1][j][k] + W[i][j-1][k]) * 0.5;
						}

				break;


				case 1<<2: //Free Slip conditions
					if ( B_O(flag[i][j][k])){
						U[i  ][j  ][k  ] = 0.0;
						V[i  ][j-1][k  ] = V[i+1][j-1][k ];
						V[i  ][j  ][k  ] = V[i+1][j  ][k ];
						W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
						W[i  ][j  ][k  ] = W[i+1][j  ][k ];
						}
					if ( B_W(flag[i][j][k])){
						U[i-1][j  ][k  ] = 0.0;
						V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
						V[i  ][j  ][k  ] = V[i-1][j  ][k  ];
						W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
						W[i  ][j  ][k  ] = W[i-1][j  ][k  ];
						}
					if ( B_N(flag[i][j][k])){
						V[i  ][j][k  ]   = 0.0;
						U[i-1][j][k  ]   = U[i-1][j+1][k  ];
						U[i  ][j][k  ]   = U[i  ][j+1][k  ];
						W[i  ][j][k-1]   = W[i  ][j+1][k-1];
						W[i  ][j][k  ]   = W[i  ][j+1][k  ];
						}
					if ( B_S(flag[i][j][k])){
						V[i  ][j-1][k  ] = 0.0;
						U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
						U[i  ][j  ][k  ] = U[i  ][j-1][k  ];
						W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
						W[i  ][j  ][k  ] = W[i  ][j-1][k  ];
						}
					if ( B_F(flag[i][j][k])){
						W[i][j][k]   = 0.0;
						U[i-1][j][k] = U[i-1][j][k+1];
						U[i][j][k]   = U[i][j][k+1];
						V[i][j-1][k] = V[i][j-1][k+1];
						V[i][j][k]   = V[i][j][k+1];
						}
					if ( B_B(flag[i][j][k])){
						W[i][j][k-1] = 0.0;
						U[i-1][j][k] = U[i-1][j][k-1];
						U[i][j][k]   = U[i][j][k-1];
						V[i][j-1][k] = V[i][j-1][k-1];
						V[i][j][k]   = V[i][j][k-1];
						}

					if ( B_NO(flag[i][j][k])){
						U[i][j][k]   = 0.0;
						U[i-1][j][k] = U[i-1][j+1][k];
						V[i][j][k]   = 0.0;
						V[i][j-1][k] = V[i+1][j-1][k];
						W[i][j][k]   = (W[i][j+1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j+1][k-1] + W[i+1][j][k-1]) * 0.5;
						}
					if ( B_NW(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k]   = U[i][j+1][k];
						V[i][j][k]   = 0.0;
						V[i][j-1][k] = V[i-1][j-1][k];
						W[i][j][k]   = (W[i][j+1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j+1][k-1] + W[i-1][j][k-1]) * 0.5;
						}
					if ( B_NF(flag[i][j][k])){
						V[i][j][k]   = 0.0;
						V[i][j-1][k] = V[i][j-1][k+1];
						W[i][j][k]   = 0.0;
						W[i][j][k-1] = W[i][j+1][k-1];
						U[i][j][k]   = (U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						}
					if ( B_NB(flag[i][j][k])){
						V[i][j][k] = 0.0;
						V[i][j-1][k] = V[i][j-1][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = W[i][j+1][k];
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						}
					if ( B_SO(flag[i][j][k])){
						U[i][j][k]   = 0.0;
						U[i-1][j][k] = U[i-1][j-1][k];
						V[i][j-1][k] = 0.0;
						V[i][j][k]   = V[i+1][j][k];
						W[i][j][k]   = (W[i][j-1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j-1][k-1] + W[i+1][j][k-1]) * 0.5;
						}
					if ( B_SW(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k]   = U[i][j-1][k];
						V[i][j-1][k] = 0.0;
						V[i][j][k]   = V[i-1][j][k];
						W[i][j][k]   = (W[i][j-1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j-1][k-1] + W[i-1][j][k-1]) * 0.5;
						}
					if ( B_SF(flag[i][j][k])){
						V[i][j-1][k] = 0.0;
						V[i][j][k]   = V[i][j][k+1];
						W[i][j][k]   = 0.0;
						W[i][j][k-1] = W[i][j-1][k-1];
						U[i][j][k]   = (U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						}
					if ( B_SB(flag[i][j][k])){
						V[i][j-1][k] = 0.0;
						V[i][j][k]   = V[i][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k]   = W[i][j-1][k];
						U[i][j][k]   = (U[i][j-1][k] + U[i][j][k-1]) * 0.5;
			  			U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						}
					if ( B_OF(flag[i][j][k])){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = U[i-1][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = W[i+1][j][k-1];
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						}
					if ( B_WF(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = U[i][j][k+1];
						W[i][j][k] = 0.0;
						W[i][j][k-1] = W[i-1][j][k-1];
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						}
					if ( B_OB(flag[i][j][k])){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = U[i-1][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = W[i+1][j][k];
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						}
					if ( B_WB(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = U[i][j][k-1];
						W[i][j][k-1] = 0.0;
						W[i][j][k] = W[i-1][j][k];
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						}

					if ( B_NOF(flag[i][j][k])){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = (W[i+1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NWF(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = (W[i-1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if (B_NOB(flag[i][j][k])){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = (W[i+1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_NWB(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						V[i][j][k] = 0.0;
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = (W[i-1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_SOF(flag[i][j][k])){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = (W[i+1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if ( B_SWF(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = (U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = 0.0;
						W[i][j][k-1] = (W[i-1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if ( B_SOB(flag[i][j][k])){
						U[i][j][k] = 0.0;
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = (W[i+1][j][k] + W[i][j-1][k]) * 0.5;
						}
					if ( B_SWB(flag[i][j][k])){
						U[i-1][j][k] = 0.0;
						U[i][j][k] = (U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						V[i][j-1][k] = 0.0;
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = 0.0;
						W[i][j][k] = (W[i-1][j][k] + W[i][j-1][k]) * 0.5;
						}			

					break;

				case 1<<3: // outflow 
					if ( B_O(flag[i][j][k])){
						U[i][j][k]   = U[i+1][j][k];

						V[i][j-1][k] = V[i+1][j-1][k];
						V[i][j][k]   = V[i+1][j][k];

						W[i][j][k-1] = W[i+1][j][k-1];
						W[i][j][k]   = W[i+1][j][k];
						}
					if ( B_W(flag[i][j][k])){ 
						U[i-1][j][k] = U[i-2][j][k];

						V[i][j-1][k] = V[i-1][j-1][k];
						V[i][j][k]   = V[i-1][j][k];

						W[i][j][k-1] = W[i-1][j][k-1];
						W[i][j][k]   = W[i-1][j][k];
						}
					if ( B_N(flag[i][j][k])){
						V[i][j][k] = V[i][j+1][k];

						U[i-1][j][k] = U[i-1][j+1][k];
						U[i][j][k] = U[i][j+1][k];

						W[i][j][k-1] = W[i][j+1][k-1];
						W[i][j][k] = W[i][j+1][k];
						}
					if ( B_S(flag[i][j][k])){
						V[i][j-1][k] = V[i][j-2][k];
						U[i-1][j][k] = U[i-1][j-1][k];
						U[i][j][k] = U[i][j-1][k];
						W[i][j][k-1] = W[i][j-1][k-1];
						W[i][j][k] = W[i][j-1][k];

						}
					if ( B_F(flag[i][j][k])){
						W[i][j][k] = W[i][j][k+1];
						U[i-1][j][k] = U[i-1][j][k+1];
						U[i][j][k] = U[i][j][k+1];
						V[i][j-1][k] = V[i][j-1][k+1];
						V[i][j][k] = V[i][j][k+1];
						}
					if ( B_B(flag[i][j][k])){
						W[i][j][k-1] = W[i][j][k-2];
						U[i-1][j][k] = U[i-1][j][k-1];
						U[i][j][k] = U[i][j][k-1];
						V[i][j-1][k] = V[i][j-1][k-1];
						V[i][j][k] = V[i][j][k-1];
						}

					if ( B_NO(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = U[i-1][j+1][k];
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = V[i+1][j-1][k];
						W[i][j][k] = (W[i][j+1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j+1][k-1] + W[i+1][j][k-1]) * 0.5;
						}
					if ( B_NW(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = U[i][j+1][k];
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = V[i-1][j-1][k];
						W[i][j][k] = (W[i][j+1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j+1][k-1] + W[i-1][j][k-1]) * 0.5;
						}
					if ( B_NF(flag[i][j][k])){
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = V[i][j-1][k+1];
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = W[i][j+1][k-1];
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						}
					if ( B_NB(flag[i][j][k])){
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = V[i][j-1][k-1];
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = W[i][j+1][k];
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						}
					if ( B_SO(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = U[i-1][j-1][k];
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = V[i+1][j][k];
						W[i][j][k] = (W[i][j-1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j-1][k-1] + W[i+1][j][k-1]) * 0.5;
						}
					if ( B_SW(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = U[i][j-1][k];
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = V[i-1][j][k];
						W[i][j][k] = (W[i][j-1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = (W[i][j-1][k-1] + W[i-1][j][k-1]) * 0.5;
						}
					if ( B_SF(flag[i][j][k])){
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = V[i][j][k+1];
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = W[i][j-1][k-1];
						U[i][j][k] = (U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						}
					if ( B_SB(flag[i][j][k])){
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = V[i][j][k-1];
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = W[i][j-1][k];
						U[i][j][k] = (U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						}
					if ( B_OF(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = U[i-1][j][k+1];
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = W[i+1][j][k-1];
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						}
					if ( B_WF(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = U[i][j][k+1];
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = W[i-1][j][k-1];
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						}
					if ( B_OB(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = U[i-1][j][k-1];
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = W[i+1][j][k];
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						}
					if ( B_WB(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = U[i][j][k-1];
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = W[i-1][j][k];
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						}

					if ( B_NOF(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = (W[i+1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NWF(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = (W[i-1][j][k-1] + W[i][j+1][k-1]) * 0.5;
						}
					if ( B_NOB(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = (U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = (V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = (W[i+1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_NWB(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = (U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						V[i][j][k] = V[i][j+1][k];
						V[i][j-1][k] = (V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = (W[i-1][j][k] + W[i][j+1][k]) * 0.5;
						}
					if ( B_SOF(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = (W[i+1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if ( B_SWF(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = (U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						W[i][j][k] = W[i][j][k+1];
						W[i][j][k-1] = (W[i-1][j][k-1] + W[i][j-1][k-1]) * 0.5;
						}
					if ( B_SOB(flag[i][j][k])){
						U[i][j][k] = U[i+1][j][k];
						U[i-1][j][k] = (U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = (V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = (W[i+1][j][k] + W[i][j-1][k]) * 0.5;
						}
					if ( B_SWB(flag[i][j][k])){
						U[i-1][j][k] = U[i-2][j][k];
						U[i][j][k] = (U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						V[i][j-1][k] = V[i][j-2][k];
						V[i][j][k] = (V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						W[i][j][k-1] = W[i][j][k-2];
						W[i][j][k] = (W[i-1][j][k] + W[i][j-1][k]) * 0.5;
						}
				break;

				case 1<<12:   //moving wall
						
					if ( B_O(flag[i][j][k])){ // wall is O. its moving direction is +/-y
						U[i][j][k] = 0.0;
						V[i][j-1][k] = 2.0*velMW[1] - V[i+1][j-1][k];
						V[i][j][k]   = 2.0*velMW[1] - V[i+1][j][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i+1][j][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i+1][j][k];
						}
					if ( B_W(flag[i][j][k])){ // wall is W. its moving direction is +/-y
						U[i-1][j][k] = 0.0;
						V[i][j-1][k] = 2.0*velMW[1] - V[i-1][j-1][k];
						V[i][j][k]   = 2.0*velMW[1] - V[i-1][j][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i-1][j][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i-1][j][k];
						}
					if ( B_N(flag[i][j][k])){ // wall is N. its moving direction is +/-z
						V[i][j][k] = 0.0;
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j+1][k];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j+1][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i][j+1][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i][j+1][k];
						}
					if ( B_S(flag[i][j][k])){ // wall is S. its moving direction is +/-z
						V[i][j-1][k] = 0.0;
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j-1][k];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j-1][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i][j-1][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i][j-1][k];
						}
					if ( B_F(flag[i][j][k])){ // wall is U. its moving direction is +/-x
						W[i][j][k]   = 0.0;
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j][k+1];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j][k+1];
						V[i][j-1][k] = 2.0*velMW[2] - V[i][j-1][k+1];
						V[i][j][k]   = 2.0*velMW[2] - V[i][j][k+1];
						}
					if ( B_B(flag[i][j][k])){ // wall is D. its moving direction is +/-x
						W[i][j][k-1] = 0.0;
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j][k-1];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j][k-1];
						V[i][j-1][k] = 2.0*velMW[1] - V[i][j-1][k-1];
						V[i][j][k]   = 2.0*velMW[1] - V[i][j][k-1];
						}

					if ( B_NO(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+y,-x,0) or (-y,+x,0)
						U[i][j][k]   = 2.0*velMW[0] - U[i][j+1][k];
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j+1][k];
						V[i][j][k]   = 2.0*velMW[1] - V[i+1][j][k];
						V[i][j-1][k] = 2.0*velMW[1] - V[i+1][j-1][k];
						W[i][j][k]   = - (W[i][j+1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = - (W[i][j+1][k-1] + W[i+1][j][k-1]) * 0.5;
						}
					if ( B_NW(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+y,+x,0) or (-y,-x,0)
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j+1][k];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j+1][k];
						V[i][j][k]   = 2.0*velMW[1] - V[i-1][j][k];
						V[i][j-1][k] = 2.0*velMW[1] - V[i-1][j-1][k];
						W[i][j][k] = - (W[i][j+1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = - (W[i][j+1][k-1] + W[i-1][j][k-1]) * 0.5;
						}
					if ( B_NF(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (0,-z,+y) or (0,+z,-y)
						V[i][j][k] = 2.0*velMW[1] - V[i][j][k+1];
						V[i][j-1][k] = 2.0*velMW[1] - V[i][j-1][k+1];
						W[i][j][k] = 2.0*velMW[2] - W[i][j+1][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i][j+1][k-1];
						U[i][j][k] = - (U[i][j+1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = - (U[i-1][j+1][k] + U[i-1][j][k+1]) * 0.5;
						}
					if ( B_NB(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (0,+z,+y) or (0,-z,-y)
						V[i][j][k]   = 2.0*velMW[1] - V[i][j][k-1];
						V[i][j-1][k] = 2.0*velMW[1] - V[i][j-1][k-1];
						W[i][j][k-1] = 2.0*velMW[2] - W[i][j+1][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i][j+1][k];
						U[i][j][k]   = -(U[i][j+1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = -(U[i-1][j+1][k] + U[i-1][j][k-1]) * 0.5;
						}
					if ( B_SO(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+y,+x,0) or (-y,-x,0)
						U[i][j][k]   = 2.0*velMW[0] - U[i][j-1][k];
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j-1][k];
						V[i][j-1][k] = 2.0*velMW[1] - V[i+1][j-1][k];
						V[i][j][k]   = 2.0*velMW[1] - V[i+1][j][k];
						W[i][j][k]   = - (W[i][j-1][k] + W[i+1][j][k]) * 0.5;
						W[i][j][k-1] = - (W[i][j-1][k-1] + W[i+1][j][k-1]) * 0.5;
						}
					if ( B_SW(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+y,-x,0) or (-y,+x,0)
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j-1][k];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j-1][k];
						V[i][j-1][k] = 2.0*velMW[1] - V[i-1][j-1][k];
						V[i][j][k]   = 2.0*velMW[1] - V[i-1][j][k];
						W[i][j][k]   = - (W[i][j-1][k] + W[i-1][j][k]) * 0.5;
						W[i][j][k-1] = - (W[i][j-1][k-1] + W[i-1][j][k-1]) * 0.5;
						}
					if ( B_SF(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (0,+z,+y) or (0,-z,-y)
						V[i][j-1][k] = 2.0*velMW[1] - V[i][j-1][k+1];
						V[i][j][k]   = 2.0*velMW[1] - V[i][j][k+1];
						W[i][j][k]   = 2.0*velMW[2] - W[i][j-1][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i][j-1][k-1];
						U[i][j][k]   = - (U[i][j-1][k] + U[i][j][k+1]) * 0.5;
						U[i-1][j][k] = - (U[i-1][j-1][k] + U[i-1][j][k+1]) * 0.5;
						}
					if ( B_SB(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (0,-z,+y) or (0,+z,-y)
						V[i][j-1][k] = 2.0*velMW[1] - V[i][j-1][k-1];
						V[i][j][k]   = 2.0*velMW[1] - V[i][j][k-1];
						W[i][j][k-1] = 2.0*velMW[2] - W[i][j-1][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i][j-1][k];
						U[i][j][k]   = - (U[i][j-1][k] + U[i][j][k-1]) * 0.5;
						U[i-1][j][k] = - (U[i-1][j-1][k] + U[i-1][j][k-1]) * 0.5;
						}
					if ( B_OF(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+z,0,-x) or (-z,0,+x)
						U[i][j][k]   = 2.0*velMW[0] - U[i][j][k+1];
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j][k+1];
						W[i][j][k]   = 2.0*velMW[2] - W[i+1][j][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i+1][j][k-1];
						V[i][j][k]   = - (V[i+1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = - (V[i+1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						}
					if ( B_WF(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+z,0,+x) or (-z,0,-x)
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j][k+1];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j][k+1];
						W[i][j][k]   = 2.0*velMW[2] - W[i-1][j][k];
						W[i][j][k-1] = 2.0*velMW[2] - W[i-1][j][k-1];
						V[i][j][k]   = - (V[i-1][j][k] + V[i][j][k+1]) * 0.5;
						V[i][j-1][k] = - (V[i-1][j-1][k] + V[i][j-1][k+1]) * 0.5;
						}
					if ( B_OB(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+z,0,+x) or (-z,0,-x)
						U[i][j][k]   = 2.0*velMW[0] - U[i][j][k-1];
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j][k-1];
						W[i][j][k-1] = 2.0*velMW[2] - W[i+1][j][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i+1][j][k];
						V[i][j][k]   = - (V[i+1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = - (V[i+1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						}
					if ( B_WB(flag[i][j][k])){ // circular moving direction of (O/W, N/S, U/D) is (+z,0,-x) or (-z,0,+x)
						U[i-1][j][k] = 2.0*velMW[0] - U[i-1][j][k-1];
						U[i][j][k]   = 2.0*velMW[0] - U[i][j][k-1];
						W[i][j][k-1] = 2.0*velMW[2] - W[i-1][j][k-1];
						W[i][j][k]   = 2.0*velMW[2] - W[i-1][j][k];
						V[i][j][k]   = - (V[i-1][j][k] + V[i][j][k-1]) * 0.5;
						V[i][j-1][k] = - (V[i-1][j-1][k] + V[i][j-1][k-1]) * 0.5;
						}

				break;


			}
		}
	}
}
spec_boundary_val(imax,jmax,kmax,U,V, W,flag,velIN);
}


