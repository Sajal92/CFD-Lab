
#include "helper.h"
#include "init.h"
#include <stdio.h>

int read_parameters( const char *szFileName,       /* name of the file */
		double *Re,                /* reynolds number   */
		double *UI,                /* velocity x-direction */
		double *VI,                /* velocity y-direction */
		double *WI,                /* velocity z-direction */
		double *PI,                /* pressure */
		double *GX,                /* gravitation x-direction */
		double *GY,                /* gravitation y-direction */
		double *GZ,                /* gravitation z-direction */
		double *t_end,             /* end time */
		double *xlength,           /* length of the domain x-dir.*/
		double *ylength,           /* length of the domain y-dir.*/
		double *zlength,           /* length of the domain z-dir.*/
		double *dt,                /* time step */
		double *dx,                /* length of a cell x-dir. */
		double *dy,                /* length of a cell y-dir. */
		double *dz,                /* length of a cell z-dir. */
		int  *imax,                /* number of cells x-direction*/
		int  *jmax,                /* number of cells y-direction*/
		int  *kmax,                /* number of cells z-direction*/
		double *alpha,             /* uppwind differencing factor*/
		double *omg,               /* relaxation factor */
		double *tau,               /* safety factor for time step*/
		int  *itermax,             /* max. number of iterations for pressure per time step */
		double *eps,               /* accuracy bound for pressure*/
		double *dt_value,          /* time for output */
		char *problemGeometry,     /*problem to solve*/
		double *velIN,             /*velocity of inflow*/
		double *velMW		   /*velocity of the Wall*/
)
{
	READ_DOUBLE( szFileName, *xlength );
	READ_DOUBLE( szFileName, *ylength );
	READ_DOUBLE( szFileName, *zlength );

	READ_INT   ( szFileName, *imax );
	READ_INT   ( szFileName, *jmax );
	READ_INT   ( szFileName, *kmax );

	READ_DOUBLE( szFileName, *dt    );
	READ_DOUBLE( szFileName, *t_end );
	READ_DOUBLE( szFileName, *tau   );

	READ_DOUBLE( szFileName, *dt_value );
	READ_INT   ( szFileName, *itermax );

	READ_DOUBLE( szFileName, *eps   );
	READ_DOUBLE( szFileName, *omg   );
	READ_DOUBLE( szFileName, *alpha );

	READ_DOUBLE( szFileName, *Re    );

	READ_DOUBLE( szFileName, *GX );
	READ_DOUBLE( szFileName, *GY );
	READ_DOUBLE( szFileName, *GZ );

	READ_DOUBLE( szFileName, *PI );
	READ_DOUBLE( szFileName, *UI );
	READ_DOUBLE( szFileName, *VI );
	READ_DOUBLE( szFileName, *WI );

	READ_STRING( szFileName, problemGeometry );
	READ_DOUBLE( szFileName, *velIN );

	double *velMWx = &velMW[0], *velMWy = &velMW[1], *velMWz = &velMW[2];
	READ_DOUBLE( szFileName, *velMWx );
	READ_DOUBLE( szFileName, *velMWy );
	READ_DOUBLE( szFileName, *velMWz );
	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);
	*dz = *zlength / (double)(*kmax);
	
	return 1;
}

void init_uvwp(
		double UI,
		double VI,
		double WI,
		double PI,
		int ***flag,
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***P
)
{
printf("Initialization of U,V,W,P,T ... \n");
	
	init_matrix2(U, 0, imax+1, 0, jmax+1, 0, kmax+1, UI);
	init_matrix2(V, 0, imax+1, 0, jmax+1, 0, kmax+1, VI);
	init_matrix2(W, 0, imax+1, 0, jmax+1, 0, kmax+1, WI);
	init_matrix2(P, 0, imax+1, 0, jmax+1, 0, kmax+1, PI);


	for(int i = 1; i <= imax+1; i++) {
		for(int j = 1; j <= jmax+1; j++) {
			for(int k = 1; k <= kmax+1; k++) {
				if((flag[i][j][k] &(1<<0)) == 0){
					U[i  ][j  ][k  ] = 0;
					U[i-1][j  ][k  ] = 0;
					V[i  ][j  ][k  ] = 0;
					V[i  ][j-1][k  ] = 0;
					W[i  ][j  ][k  ] = 0;
					W[i  ][j  ][k-1] = 0;
					P[i  ][j  ][k  ] = 0;
				}
			}
		}
	}
	printf("U,V,W,P,T matrices have been initialized... \n \n");
}


int  Fluid(int pic){
	if((pic == 4)) {
		return 1;
	}
		else {
		return 0;
	}
}



void assert_error()
{
	char szBuff[80];
	sprintf( szBuff, "Geometry is forbidden. Consider modifying .pgm file. \n");
	ERROR( szBuff );
}



//Avoids any forbidden configuration
/*void forbid_assert(int imax, int jmax, int kmax, int **pic)
{
int **pic1 = read_pgm(geometry);
 

int counter = 0;

    //Checking forbidden configuration
    for(int i=1; i<=imax ; i++)
    {

        for(int j=1; j<=jmax; j++)
        {
	   for (int k=1; k<=kmax; k++)
	   {

            counter = 0;
            if(pic1[i][j + k*(jmax+2)] != 4 && pic1[i][j + k*(jmax+2)] != 3 && pic1[i][j + k*(jmax+2)] != 2)
	    {
		    if((pic1[i+1][j + k*(jmax+2)] == 4) || (pic1[i+1][j + k*(jmax+2)] == 3) || (pic1[i+1][j + k*(jmax+2)] == 2))
		    {
		    counter++;   
		    }
		    if((pic1[i-1][j + k*(jmax+2)] == 4) || (pic1[i-1][j + k*(jmax+2)] == 3) || (pic1[i-1][j + k*(jmax+2)] == 2))
		    {
		    counter++;   
		    }
		    if((pic1[i][j+1 + k*(jmax+2)] == 4) || (pic1[i][j+1 + k*(jmax+2)] == 3) || (pic1[i][j+1 + k*(jmax+2)] == 2))
		    {
		    counter++;   
		    }
		    if((pic1[i][j-1 + k*(jmax+2)] == 4) || (pic1[i][j-1 + k*(jmax+2)] == 3) || (pic1[i][j-1 + k*(jmax+2)] == 2))
		    {
		    counter++;   
		    }
	    //include forbidden conditions for front and back
	    }

            if(counter > 3)
            {
            assert_error();
            }


        }
    }
}
}*/


void init_flag( char* geometry, int imax, int jmax, int kmax, int ***flag)
{
	printf("Flags are being set... \n");
	int **pic = read_pgm(geometry);
	int ***Q = imatrix2(0,imax+1,0,jmax+1,0,kmax+1);
	for (int i=0; i<=imax+1; i++){
		for (int j=0; j<=jmax+1; j++){
			for (int k=0; k<=kmax+1; k++) {
				Q[i][j][k]=pic[i][j+k*(jmax+2)];
			}
		}
	}

	//forbid_assert(imax, jmax, pic); //change this one
	
	for (int i=0; i<=imax+1; i++){
		for (int j=0; j<=jmax+1; j++){
			for (int k=0; k<=kmax+1; k++) {

			flag[i][j][k] = 0;

			switch(Q[i][j][k]){
				case 0: //no-slip
				flag[i][j][k] = 1<<1;
				break;

				case 1: //free-slip
				flag[i][j][k] = 1<<2;
				break;

				case 2: //outflow
				flag[i][j][k] = 1<<3;
				break;

				case 3: //inflow
				flag[i][j][k] = 1<<4;
				break;

				case 4: //fluid
				flag[i][j][k] = 1<<0;
				break;

				case 6: //moving_wall
				flag[i][j][k] = 1<<12;
				break;
			}

			if((Q[i][j][k] != 4)){ //set boundaries if not fluid

				if(i<imax+1 && Fluid(Q[i+1][j][k])){ //east
					flag[i][j][k] |= 1<<8;
				}
				if( i>0 && Fluid(Q[i-1][j][k])){ //west
					flag[i][j][k] |= 1<<7;
				}
				if(j<jmax+1 && Fluid(Q[i][j+1][k])){ //north
					flag[i][j][k] |= 1<<5;
				}
				if(j>0 && Fluid(Q[i][j-1][k])){ //south
					flag[i][j][k] |= 1<<6;
				}
				if(k<kmax+1 && Fluid(Q[i][j][k+1])){ //front
					flag[i][j][k] |= 1<<9;  
				}
				if(k>0 && Fluid(Q[i][j][k-1])){ //back
					flag[i][j][k] |= 1<<10;  
				}
			}
 

		}

	}

}	

	free_imatrix(pic, 0,imax+1,0,jmax+1);
	printf("Flags set according .pgm file...\n \n");

}
