#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

#include<math.h>
/**
 * The boundary values of the problem are set.
 */
void spec_boundary_val(int imax,int jmax, int kmax, double ***U,double ***V, double***W,int ***flag, double velIN);

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
);

int B_O(int flag);

int B_W(int flag);

int B_N(int flag);

int B_S(int flag);

int B_F(int flag);

int B_B(int flag);

int B_NO(int flag);

int B_NW(int flag);

int B_SO(int flag);

int B_SW(int flag);

int B_NF(int flag);

int B_NB(int flag);

int B_SF(int flag);

int B_SB(int flag);

int B_OF(int flag);

int B_WF(int flag);

int B_OB(int flag);

int B_WB(int flag);

int B_NOF(int flag);

int B_NOB(int flag);

int B_NWF(int flag);

int B_NWB(int flag);

int B_SOF(int flag);

int B_SOB(int flag);

int B_SWF(int flag);

int B_SWB(int flag);


#endif
