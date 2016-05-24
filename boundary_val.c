#include "boundary_val.h"
void boundaryvalues(int imax,int jmax,double **U,double **V, double** P, double** G, double** F)
{
	// wl wr wt wb ..... object wall NO_SLIP.
/*

1 for no-slip 
2 for free-slip
3 for outflow
 
TODO 
change this function to accept flag field.
impliment spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V), for the inlet boundary cond. 

 if(wl==1){ }
else if(wl==2){ }
else(wl==3){ }

 if(wr==1){ }
else if(wr==2){ }
else(wr==3){ }

 if(wt==1){ }
else if(wt==2){ }
else(wt==3){ }

 if(wb==1){ }
else if(wb==2){ }
else(w==3){ }

object NO_SLIP BC  make use of flag field for fluid obstacle distinction.  




*/





	for(int i=1; i<=imax; i++)
	{
				
		U[i][0] = -U[i][1];
		U[i][jmax+1] = 2.0 - U[i][jmax];
        V[i][0] = -V[i][1];
		V[i][jmax+1] = -V[i][jmax];
        P[i][0] = P[i][1];
		P[i][jmax+1] = P[i][jmax];
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax]; 
    }
    
	for(int j=1; j<=jmax; j++)
	{
		U[0][j] = -U[1][j];
		U[imax+1][j] = -U[imax][j];
        V[0][j] = -V[1][j];
		V[imax+1][j] = -V[imax][j];
        P[0][j] = P[1][j];
		P[imax+1][j] = P[imax][j];
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
	}
}
