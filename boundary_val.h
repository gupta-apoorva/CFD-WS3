#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int  wl,
  int wr,
  int wt,
  int wb,
  double **U,
  double **V,
  double** P,
  double** G, 
  double** F,
  int **FLAG
);

#endif
