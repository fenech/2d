#include "header.h"

inline double frank(double **grid, int ni[], int nj[])
{
    double xcell[3] = {
	grid[nj[1]][ni[0]],
	grid[nj[1]][ni[1]],
	grid[nj[1]][ni[2]]
    };
    double ycell[3] = {
	grid[nj[0]][ni[1]],
	grid[nj[1]][ni[1]],
	grid[nj[2]][ni[1]]
    };
    
    if (cos(xcell[2] - xcell[1]) < 0) xcell[2] += PI;
    if (cos(xcell[1] - xcell[0]) < 0) xcell[0] += PI;
    if (cos(ycell[2] - ycell[1]) < 0) ycell[2] += PI;
    if (cos(ycell[1] - ycell[0]) < 0) ycell[0] += PI;
		
    double dsindxr = (sin(xcell[2]) - sin(xcell[1])) / CELL_LEN,
	dsindxl = (sin(xcell[0]) - sin(xcell[1])) / CELL_LEN,
	dcosdxr = (cos(xcell[2]) - cos(xcell[1])) / CELL_LEN,
	dcosdxl = (cos(xcell[0]) - cos(xcell[1])) / CELL_LEN,
	dsindyu = (sin(ycell[2]) - sin(ycell[1])) / CELL_LEN,
	dsindyd = (sin(ycell[0]) - sin(ycell[1])) / CELL_LEN,
	dcosdyu = (cos(ycell[2]) - cos(ycell[1])) / CELL_LEN,
	dcosdyd = (cos(ycell[0]) - cos(ycell[1])) / CELL_LEN;
		
    double splay = ((dcosdxr + dsindyu) * (dcosdxr + dsindyu) + (dcosdxl + dsindyu) * (dcosdxl + dsindyu) 
		    + (dcosdxr + dsindyd) * (dcosdxr + dsindyd) + (dcosdxl + dsindyd) * (dcosdxl + dsindyd)) / 4;
    double bend = ((dsindxl - dcosdyu) * (dsindxl - dcosdyu) 
		   + (dsindxr - dcosdyu) * (dsindxr - dcosdyu)
		   + (dsindxl - dcosdyd) * (dsindxl - dcosdyd)
		   + (dsindxr - dcosdyd) * (dsindxr - dcosdyd)) / 4;
    
    return (K_1 / 2) * splay + (K_3 / 2) * bend;
}


double locfunc(double **grid, int x, int y, int flag)
{
    /* x and y are the coordinates of the modified point in the monte carlo simulation */
    /* ci and cj are lists of central points, around which the new derivatives must be calculated */
    /* these extend from x - 1 to x + 1, same for y */
    int ni[3], nj[3];
    double lF = 0;
    int i, j;
    int ci, cj;
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (flag == 1) {
	h = ny;
	numProcs = 1;
    }
  
    int c[5][2] = {
	{ x-1, y-1 },
	{ x-1, y }, 
	{ x, y },
	{ x+1, y },
	{ x+1, y+1 }
    };

    for (j = 0; j < 5; ++j) {
	ci = c[i][0];
	cj = c[i][1];
      
	ni[0] = ci - 1;
	ni[1] = ci;
	ni[2] = ci + 1;		
	nj[0] = cj - 1;
	nj[1] = cj;
	nj[2] = cj + 1;
				
	/* sort out any indices overflowing boundaries */
	for (i = 0; i <= 2; ++i) {
	    if (ni[i] < 1) ni[i] += nx;
	    if (ni[i] > nx) ni[i] -= nx;
	    if (rank == numProcs - 1 && nj[i] > h + 1) nj[i] = h + 1;
	    if (rank == 0 && nj[i] < 0) nj[i] = 0;
	}
	lF += frank(grid, ni, nj);
    }
  
    return lF;
}	

double func(double **grid, int **lock, int flag)
{
    double lF = 0.0, f = 0.0;
    int i;
    int ni[3], nj[3];
    int ci, cj;
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    /* in conjugate gradient only root computes energy for whole grid */
    if (flag == 1)
	h = ny;

    for (cj = 1; cj <= h; ++cj) {
	for (ci = 1; ci <= nx; ++ci) {
	    if (lock[cj][ci] == 0) {
		ni[0] = ci - 1;
		ni[1] = ci;
		ni[2] = ci + 1;		
		nj[0] = cj - 1;
		nj[1] = cj;
		nj[2] = cj + 1;
				
		/* sort out any indices overflowing boundaries */
		for (i = 0; i <= 2; ++i) {
		    if (ni[i] < 1) ni[i] += nx;
		    else if (ni[i] > nx) ni[i] -= nx;
		}

		lF += frank(grid, ni, nj);
	    }
	}
    }
    if (flag == 0) {
	/* monte carlo */
	MPI_Reduce(&lF, &f, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  
    else if (flag == 1) {
	/* conjugate gradient */
	f = lF;
    }
    return f;
}

void dfunc(double **fgrid, int **flock, double **dgrid)
{
    /* calculates the MINUS grad */
    int ni[3], nj[3];
    int i;
    int ci, cj;
    double xcell[3], ycell[3];
    double dsindxl, dcosdxl, dsindyd, dcosdyd;
    double dsindxr, dcosdxr, dsindyu, dcosdyu;
    double ds, db;
    
    for (cj = 1; cj <= ny; ++cj) {
	for (ci = 1; ci <= nx; ++ci) {
	    if (flock[cj][ci] == 0) {
		ni[0] = ci - 1;
		ni[1] = ci;
		ni[2] = ci + 1;		
		nj[0] = cj - 1;
		nj[1] = cj;
		nj[2] = cj + 1;
	
		/* sort out any indices overflowing boundaries */
		for (i = 0; i <= 2; ++i) {
		    if (ni[i] < 1)
			ni[i] += nx;
		    if (ni[i] > nx)
			ni[i] -= nx;
		    if (nj[i] > ny + 1)
			nj[i] = ny + 1;
		    if (nj[i] < 0)
			nj[i] = 0;
		}

		xcell[2] = fgrid[nj[1]][ni[2]];
		xcell[1] = fgrid[nj[1]][ni[1]];
		xcell[0] = fgrid[nj[1]][ni[0]];
		ycell[2] = fgrid[nj[2]][ni[1]];
		ycell[1] = fgrid[nj[1]][ni[1]];
		ycell[0] = fgrid[nj[0]][ni[1]];
		
		if (cos(xcell[2] - xcell[1]) < 0)
		    xcell[2] += PI;
		if (cos(xcell[1] - xcell[0]) < 0)
		    xcell[0] += PI;
		if (cos(ycell[2] - ycell[1]) < 0)
		    ycell[2] += PI;
		if (cos(ycell[1] - ycell[0]) < 0)
		    ycell[0] += PI;
	
		dsindxr = sin(xcell[2]) - sin(xcell[1]);
		dsindxl = sin(xcell[0]) - sin(xcell[1]);
		dcosdxr = cos(xcell[2]) - cos(xcell[1]);
		dcosdxl = cos(xcell[0]) - cos(xcell[1]);
		dsindyu = sin(ycell[2]) - sin(ycell[1]);
		dsindyd = sin(ycell[0]) - sin(ycell[1]);
		dcosdyu = cos(ycell[2]) - cos(ycell[1]);
		dcosdyd = cos(ycell[0]) - cos(ycell[1]);
		
		ds = (sin(xcell[1]) + cos(xcell[1])) * (dcosdxr - dsindyu + dcosdxl - dsindyd); 
		db = (-cos(xcell[1]) - sin(xcell[1])) * (dsindxl - dcosdyu + dsindxr - dcosdyd);

		dgrid[cj][ci] = -ds - db;
	    }
	}
    }  
}
