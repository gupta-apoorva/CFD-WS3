#include "boundary_val.h"

#DEFINE NO_SLIP 1
#DEFINE FREE_SLIP 2
#DEFINE OUTFLOW 3

void boundaryvalues(int imax,int jmax,double **U,double **V, double** P, double** G, double** F)
{
	// wl wr wt wb ..... object wall NO_SLIP.
//TODO 
//change this function to accept flag field.
//impliment spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V), for the inlet boundary cond. 

//object NO_SLIP BC  make use of flag field for fluid obstacle distinction.  
if (wl==NO_SLIP)
{
	for (int j = 1; j <jmax ; ++j)
	{
		V[0][j] = -V[1][j];
		U[0][j] = 0;
	}
}
else if (wl==FREE_SLIP)
{
	for for (int j = 1; j <jmax ; ++j)
	{
		V[0][j] = V[1][j];
		U[0][j] = 0;
	}
}
else if (wl==OUTFLOW))
{
	for (int j = 1; j <jmax ; ++j)
	{
		V[0][j] = V[1][j];
		U[0][j] = U[1][j];
	}
}

if (wr==NO_SLIP)
{
	for (int j = 1; j <jmax ; ++j)
	{
		V[imax+1][j] = -V[imax][j];
		U[imax][j] = 0;
	}
}
else if (wr==FREE_SLIP)
{
	for for (int j = 1; j <jmax ; ++j)
	{
		V[imax+1][j] = V[imax][j];
		U[imax][j] = 0;
	}
}
else if (wr==OUTFLOW))
{
	for (int j = 1; j <jmax ; ++j)
	{
		V[imax+1][j] = V[imax][j];
		U[imax][j] = U[imax-1][j];
	}
}

if (wt==NO_SLIP)
{
	for (int i = 1; i <imax ; ++i)
	{
		U[i][jmax+1] = -U[i][jmax];
		V[i][jmax] = 0;
	}
}
else if (wt==FREE_SLIP)
{
	for for (int i = 1; i <imax ; ++i)
	{
		U[i][jmax+1] = U[i][jmax];
		V[i][jmax] = 0;
	}
}
else if (wt==OUTFLOW))
{
	for (int i = 1; i <imax ; ++i)
	{
		U[i][jmax+1] = U[i][jmax];
		V[i][jmax] = V[i][jmax-1];
	}
}

if (wb==NO_SLIP)
{
	for (int i = 1; i <imax ; ++i)
	{
		U[i][0] = -U[i][1];
		V[i][0] = 0;
	}
}
else if (wb==FREE_SLIP)
{
	for for (int i = 1; i <imax ; ++i)
	{
		U[i][0] = U[i][1];
		V[i][0] = 0;
	}
}
else if (wb==OUTFLOW))
{
	for (int i = 1; i <imax ; ++i)
	{
		U[i][0 = U[i][1];
		V[i][0] = V[i][1];
	}
}




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
