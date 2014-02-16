#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-EPS)
/* removed other definition of EPS  */

#define EPS 1.0e-7
#define GOLD 1.618034
#define CGOLD 0.3819660

#include "mpi.h"
#include "nrutil.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265358979323846
#define K_1 1.0
#define K_3 1.0

#define DEBUG 0

int nx, ny;             /* array dimensions, taken as arguments to the program */
int h;                  /* array height per process */
int rad;                /* radius of particles */
int sep;                /* particle separation */
int id;                 /* run id */
int ba;                 /* boundary angle */

/* initialisation */
void initialise(double*** grid, int*** lock);
/* local energy calculation */
double locfunc(double **grid, int x, int y, int flag);
/* random number generator */
float ran2(long *idum);
/* total energy function */
double func(double **grid, int **lock, int flag);
/* monte carlo */
void monte(double **grid, int **lock, int maxiter, float t0);
/* grid printer */
void print(double **grid, char *fname);
/* derivative calculation - returns MINUS grad */
void dfunc(double **fgrid, int **flock, double **dgrid);
/* conjugate gradient function */
void conjgrad(double **fgrid, int **flock, int maxits, double *fret);
/* grid gatherer */
void grid2root(double **grid, double ***fgrid);
/* lock grid gatherer */
void lock2root(int **lock, int ***flock);
/* grid printer - single processor, specified range */
void prn(double **grid, int sy, int ey, int sx, int ex, char *fname);
/* prints contour map of grid */
void prncont(double **grid, char *fname);
/* reconstructs grid from output file */
void file2grid(double **grid, char *gname);
