#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include <string.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_FLAG(), output_uvp()
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




int main(int argn, char** args)
{

   double** U;
   double** V;
   double** P;
   double Re;               
   double UI;                
   double VI;               
   double PI;           
   double GX;                
   double GY;                
   double t_end;            
   double xlength;          
   double ylength;           
   double dt;                
   double dx;             
   double dy;               
   int  imax;                
   int  jmax;               
   double alpha;            
   double omg;               
   double tau;              
   int  itermax;             
   double eps;              
   double dt_value;          
   double** RS;
   double** F;
   double** G;
   int** FLAG;
   int** pgm = NULL;
   int wl;	
   int wr;
   int wt;
   int wb;
   int problemtype;
   double delta_p;
   double input_vel;
   int pType;
   
	
//setting the parameters
read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax,
                 &jmax, &alpha, &omg, &tau,&itermax, &eps, &dt_value, &wl, &wr, &wt, &wb,&pType,&delta_p,&input_vel);
printf(" pType =  %d \n", pType);

//if (strcmp(problem,"STEP") == 0)
printf("%d\n",problemtype );


pgm = read_pgm("mesh2.pgm");

// Creating the arrays U,V and P
  U = matrix ( 0 , imax+1 , 0 , jmax+1 );
  V = matrix ( 0 , imax+1 , 0 , jmax+1 );
  P = matrix ( 0 , imax+1 , 0 , jmax+1 );
        

// Creating arrays for right side of pressure poissons equation (RS) and F and G
  RS = matrix ( 0,imax+1,0,jmax+1);
  F = matrix (0,imax+1,0,jmax+1);
  G = matrix (0,imax+1,0,jmax+1);
  FLAG = imatrix (0,imax+1,0,jmax+1);

// Initializing the arrays U,V,P,RS,F and G
  init_uvp( UI, VI,PI,imax, jmax,U,V,P);
  init_matrix(RS,0,imax+1,0,jmax+1,0);
  init_matrix(F,0,imax+1,0,jmax+1,0);
  init_matrix(G,0,imax+1,0,jmax+1,0);
  init_imatrix(FLAG,0,imax+1,0,jmax+1,1);

  for (int i = 0; i < imax+2; ++i)
  {
  	for (int j = 0; j < jmax+2; ++j)
    {
  		if (pgm[i][j] == 0)
      {
  			FLAG[i][j] = 0;
  		}
  	}
  }

  for (int j = 0; j < jmax+2; ++j){
  	for (int i = 0; i < imax+2; ++i)
  		printf(" %d",FLAG[i][j]);
		printf("\n");
  		
  	}

  	int dummy[imax+4][jmax+4];

  	for (int j = 0; j < jmax+4; ++j)
    {
  	for (int i = 0; i < imax+4; ++i)
  		dummy[i][j] = 0;
  	}



  	for (int j = 0; j < jmax+2; ++j)
    {
  	for (int i = 0; i < imax+2; ++i)
  		dummy[i+1][j+1] = FLAG[i][j];
  	}

        for (int j = 0; j < jmax+4; ++j){
    for (int i = 0; i < imax+4; ++i)
      printf(" %d",dummy[i][j]);
    printf("\n");}

	for (int j = 0; j < jmax+2; ++j)
  {
		for (int i = 0; i < imax+2; ++i)
    {
			FLAG[i][j] = 16*dummy[i+1][j+1] + 8*dummy[i+2][j+1] + 4*dummy[i][j+1] + 2*dummy[i+1][j+2] + dummy[i+1][j]; // big problem here...look at the output...		
		}
	}

for (int j = 0; j < jmax+2; ++j){
  	for (int i = 0; i < imax+2; ++i)
  		printf(" %d",FLAG[i][j]);
		printf("\n");
  		
  	}


  double t=0;   // initialize the time
  int n = 0;    // number of time steps


while (t<t_end)
  {

      calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);

      boundaryvalues(imax, jmax, wl , wr, wt, wb , U, V , P, G, F, FLAG);

      spec_boundary_val (pType, imax, jmax, U, V,delta_p,input_vel, Re);
write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P);
      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G,FLAG);
      calculate_rs(dt,dx,dy, imax,jmax, F, G, RS,FLAG);
      int it = 0;
      double res = 1000;

      while(it<itermax && res > eps) 
          {
            sor(omg, dx,dy,imax,jmax, P, RS, &res,FLAG);
            it++; 
          }

      calculate_uv(dt,dx, dy,imax,jmax,U,V,F,G,P,FLAG);
      t = t+dt;
      n = n+1;
//write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P);
  }

//write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P);

free_matrix(U,0,imax+1,0,jmax+1);
free_matrix(V,0,imax+1,0,jmax+1);
free_matrix(P,0,imax+1,0,jmax+1);
free_matrix(RS,0,imax+1,0,jmax+1);
free_matrix(F,0,imax+1,0,jmax+1);
free_matrix(G,0,imax+1,0,jmax+1);
free_imatrix(FLAG,0,imax+1,0,jmax+1);



  return 0;
}
