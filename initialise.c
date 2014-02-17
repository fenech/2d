#include "header.h"

void initialise(double*** grid, int*** lock, const t_par * par, int sep)
{
    int i, j, k;                               /* counters */
    long idum;                                 /* random seed */
    int rank, numProcs;                        /* MPI variables */
    int r;                                     /* remainder used in grid allocation */
    float angle = 0;                           /* angle at top and bottom boundary */
    int bottom, top;                           /* grid boundaries indices */
    MPI_Status stat;                           
    double particle;                           
    FILE *fp = NULL;
    int *hs;
    char fname[20];
    float bangle = ba * PI / 180;               /* convert ba to radians */

    /*  char fname[10];*/
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    idum = -MPI_Wtime() - (rank + 1) - id;
    /* allocate grid and lock arrays */
    h = ny / numProcs;
    r = ny - h * numProcs;
    if (r > 0)
	for (i = 1; i <= r; ++i)
	    if (rank == i % (numProcs - 1))
		h++;  
    /* each process ranges from j = hs[rank-1] to j = hs[rank] */
    hs = ivector(0, numProcs - 1);
    MPI_Gather(&h, 1, MPI_INT, hs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(hs, numProcs, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 1; i < numProcs; ++i)
	hs[i] += hs[i-1];
    /* parts of grid unique to each process go from 1 to h
       need 2 rows either side of grid for local energy calculations */
    bottom = -1;
    top = h + 2;
    /* row 0 is lower boundary on rank 0 */
    if (rank == 0)
	++bottom;
    /* row h + 1 is upper boundary on last process */
    else if (rank == numProcs - 1)
	--top;    
    *grid = dmatrix(bottom, top, 1, nx);
    *lock = imatrix(bottom, top, 1, nx);
    /* each process initialises its unique parts */
    for (j = 1; j <= h; ++j) {
	for (i = 2; i <= nx - 1; ++i) {
	    (*grid)[j][i] = (ran2(&idum) * PI);
	    (*lock)[j][i] = 0;
	}
	/* lock left and right edges */
	(*grid)[j][1] = (*grid)[j][nx] = bangle;
	(*lock)[j][1] = (*lock)[j][nx] = 1;
    }
    /* first process initialises bottom row */
    if (rank == 0) {
	r = 0;
	angle = bangle;
    }
    /* last process initialises top row */
    if (rank == numProcs - 1) {
	r = h + 1;
	angle = bangle;
    }
    if (rank == 0 || rank == numProcs - 1)
	for (i = 1; i <= nx; ++i) {
	    (*grid)[r][i] = angle;
	    (*lock)[r][i] = 1;
	}
    if (rank == 0) {
	sprintf(fname, "par%dx%d_s%d", par[0].major, par[0].minor, sep); 
	fp = fopen(fname, "w");
    }
    
    for (k = 0; k < 2; ++k) {
	for (j = par[k].cy - par[k].major; j <= par[k].cy + par[k].major; ++j) {
	    for (i = par[k].cx - par[k].major; i <= par[k].cx + par[k].major; ++i) {
		double a = (cos(par[k].theta) * (i - par[k].cx) + sin(par[k].theta) * (j - par[k].cy)) *
		    (cos(par[k].theta) * (i - par[k].cx) + sin(par[k].theta) * (j - par[k].cy));
		double b = (sin(par[k].theta) * (i - par[k].cx) - cos(par[k].theta) * (j - par[k].cy)) *
		    (sin(par[k].theta) * (i - par[k].cx) - cos(par[k].theta) * (j - par[k].cy));
		if (a / (par[k].major * par[k].major) + b / (par[k].minor * par[k].minor) <= 1.0) {
		    if (i == par[k].cx) {
			if (strcmp(par[k].align, "para") == 0) particle = 0;
			else particle = PI / 2;
		    }
		    else {			
			particle = atan((float)(j - par[k].cy) / (i - par[k].cx));			
			if (strcmp(par[k].align, "para") == 0) particle += PI / 2;			
		    }
		    
		    if (j <= hs[rank] && j >= hs[rank-1]) {
			(*grid)[j-hs[rank-1]][i] = particle;
			(*lock)[j-hs[rank-1]][i] = 1;
		    }

		    if (rank == 0)
			fprintf(fp, "%f %f %f %f\n",
				(float)i - cos(particle) / 2,
				(float)j - sin(particle) / 2,
				cos(particle),
				sin(particle));
		}
	    }
	}
    }
	      
    if (rank == 0) fclose(fp);
  
    if (rank % 2 == 1) {
	/* odd process send downward first, then upward */
	MPI_Send(&((*grid)[1][1]), 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	if (rank < numProcs - 1) {
	    MPI_Send(&((*grid)[h-1][1]), 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    /* now receive from above, then below */
	    MPI_Recv(&((*grid)[h+1][1]), 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &stat);
	}
	MPI_Recv(&((*grid)[-1][1]), 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &stat);
    }
    if (rank % 2 == 0) {
	/* even processes */
	if (rank < numProcs - 1)
	    MPI_Recv(&((*grid)[h+1][1]), 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &stat);
	if (rank > 0)
	    MPI_Recv(&((*grid)[-1][1]), 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &stat);
	if (rank > 0)
	    MPI_Send(&((*grid)[1][1]), 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	if (rank < numProcs - 1)
	    MPI_Send(&((*grid)[h-1][1]), 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }
}
