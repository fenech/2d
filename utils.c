#include "header.h"

void print(double **grid, char *fname)
{
  /* prints full grid to file */
  int i, j;
  double **fgrid;
  FILE *fp;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  grid2root(grid, &fgrid); 

  if (rank == 0) {
    fp = fopen(fname, "w");
    for (j = 0; j <= ny + 1; ++j){
      for (i = 1; i <= nx; ++i){
	fprintf(fp, "%f %f %f %f\n",
		(float)i - (cos(fgrid[j][i])) / 2,
		(float)j - (sin(fgrid[j][i])) / 2,
		cos(fgrid[j][i]),
		sin(fgrid[j][i]));
      }
    }
    fclose(fp);
  }
  if (rank == 0)
    free_dmatrix(fgrid, 0, ny + 1, 1, nx);
}

void grid2root(double **grid, double ***fgrid)
{
  /* gathers grid to root */
  int i;
  int rank, numProcs;
  int *recvcnts = NULL, *displs = NULL;
  MPI_Status stat;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (rank == 0) {  
    displs = ivector(0, numProcs - 1);
    recvcnts = ivector(0, numProcs - 1);
  }
  *fgrid = dmatrix(0, ny + 1, 1, nx);
  MPI_Gather(&h, 1, MPI_INT, recvcnts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (i = 0; i < numProcs; ++i) {
      displs[i] = 0;
      recvcnts[i] *= nx;
    }
    for (i = 1; i < numProcs; ++i)
      displs[i] = displs[i-1] + recvcnts[i-1];
  }
  
  MPI_Gatherv(&grid[1][1], h * nx, MPI_DOUBLE, &((*fgrid)[1][1]), recvcnts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  
  /* send top row */
  if (rank == numProcs - 1)
    MPI_Send(&grid[h+1][1], nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    /* receive top row */
    MPI_Recv(&((*fgrid)[ny+1][1]), nx, MPI_DOUBLE, numProcs - 1, 0, MPI_COMM_WORLD, &stat);
    /* copy bottom row */
    for (i = 1; i <= nx; ++i)
      (*fgrid)[0][i] = grid[0][i];
  }
  if (rank == 0) {
    free_ivector(recvcnts, 0, numProcs - 1);
    free_ivector(displs, 0, numProcs - 1);
  }
}

void lock2root(int **lock, int ***flock)
{
  /* gathers lock grid to root */
  int i;
  int rank, numProcs;
  int *recvcnts = NULL, *displs = NULL;
  MPI_Status stat;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  
  if (rank == 0) {
    recvcnts = ivector(0, numProcs - 1);
    displs = ivector(0, numProcs - 1);
  }
  *flock = imatrix(0, ny + 1, 1, nx);
  MPI_Gather(&h, 1, MPI_INT, recvcnts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    for (i = 0; i < numProcs; ++i) {
      displs[i] = 0;
      recvcnts[i] *= nx;
    }
    for (i = 1; i < numProcs; ++i)
      displs[i] = displs[i-1] + recvcnts[i-1];
  }
  
  MPI_Gatherv(&lock[1][1], h * nx, MPI_INT, &((*flock)[1][1]), recvcnts, displs, MPI_INT, 0, MPI_COMM_WORLD);
  
  /* send top row */
  if (rank == numProcs - 1)
    MPI_Send(&lock[h+1][1], nx, MPI_INT, 0, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    /* receive top row */
    MPI_Recv(&((*flock)[ny+1][1]), nx, MPI_INT, numProcs - 1, 0, MPI_COMM_WORLD, &stat);
    /* copy bottom row */
    for (i = 1; i <= nx; ++i)
      (*flock)[0][i] = lock[0][i];
  }
  if (rank == 0) {
    free_ivector(recvcnts, 0, numProcs - 1);
    free_ivector(displs, 0, numProcs - 1);
  }
}

void prn(double **grid, int sy, int ey, int sx, int ex, char *gname)
{
  /* prints director field grid from single processor
     arguments:
     grid name
     start and end x and y coordinates
     filename */
  int i, j;
  FILE *fp;
  fp = fopen(gname, "w");
  for (j = sy; j <= ey; ++j){
    for (i = sx; i <= ex; ++i){
      fprintf(fp, "%f %f %f %f\n",
	      (float) i - cos(grid[j][i]) / 2.,
	      (float) j - sin(grid[j][i]) / 2.,
	      cos(grid[j][i]),
	      sin(grid[j][i]));
    }
  }
  fclose(fp);
}

void prncont(double **grid, char *gname)
{
  /* prints contour map of a grid */
  int i, j;
  FILE *fp;
  fp = fopen(gname, "w");
  for (j = 1; j <= ny; ++j){
    for (i = 1; i <= nx; ++i){
      fprintf(fp, "%d %d %f \n",
	      i, j, grid[j][i]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);	
}

void file2grid(double **grid, char *gname)
{
  /* reconstructs a grid from file */
  int i, j;
  float x;
  FILE *fp;
  fp = fopen(gname, "r");
  for (j = 1; j <= ny; ++j) 
    for (i = 1; i <= nx; ++i) {
      fscanf(fp, "%f %f %lf %f", &x, &x, &grid[j][i], &x);
      grid[j][i] = acos(grid[j][i]);
    }
  fclose(fp);
}
