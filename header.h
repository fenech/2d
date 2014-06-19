#ifndef HEADER_H
#define HEADER_H

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

static const double K_1 = 4.2e-12;
static const double K_3 = 5.3e-12;
static const double CELL_LEN = 5e-8;

#define PI M_PI
#define DEBUG 0

typedef struct {
    int cx, cy;
    int major, minor;
    float theta;
    char align[16];
} t_par;

int nx, ny;             /* array dimensions, taken as arguments to the program */
int h;                  /* array height per process */
int id;                 /* run id */

/* initialisation */
int initialise(double*** grid, int*** lock, const t_par * par, int sep, int ba, const char * suffix);
/* local energy calculation */
double locfunc(double **grid, int x, int y, int flag);
/* random number generator */
float ran2(long *idum);
/* total energy function */
double func(double **grid, int **lock, int flag);
/* monte carlo */
void monte(double **grid, int **lock, int maxiter, float t0, FILE * fp, const char * suffix);
/* grid printer */
void print(double **grid, char *fname);
/* derivative calculation - returns MINUS grad */
void dfunc(double **fgrid, int **flock, double **dgrid);
/* conjugate gradient function */
void conjgrad(double **fgrid, int **flock, int maxits, double *fret, FILE * fp, const char * suffix, const t_par * par, int sep);
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

#endif /* HEADER_H */
