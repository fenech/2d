#include "header.h"

void monte(double **grid, int **lock, int maxiter, float t0) 
{
  unsigned long t;                   /* iteration counter */
  unsigned long accloc, acctot;      /* number of accepts */
  double beta;                       /* boltzmann factor */
  double gamma = PI / 2;             /* angular variation */
  long idum;                         /* random seed */
  int i, j;                          /* grid coodinate */
  double angle;                      /* original angle */
  double e, e2;                      /* old and new energy */
  float num, p = 0;                  /* random number and acceptance probability */
  int rank, numProcs;                /* processor ID and number of processes */
  int pass;                          /* 1 or 2 */
  int srow = 0, erow = 0;            /* start row and end row */
  char fname[30];
  MPI_Status stat;
  FILE *fp = NULL;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  idum = -MPI_Wtime() - (rank + 1) - id;
  if (rank == 0) {
    sprintf(fname, "logr%ds%da%d_%d", rad, sep, ba, id);
    fp = fopen(fname, "w");
  }

  accloc = acctot = 0;
  beta = t0;

  for (t = 1; t <= maxiter; ++t) {
    if (t % 10000 == 0) {
      sprintf(fname, "gridr%ds%da%d_%d", rad, sep, ba, id);
      print(grid, fname);
    }
    
    if (t % 280 == 0)    
      beta *= 1.01;

    /* gamma selected to keep acceptance ratio at 1/2 */
    if (t % 100 == 0) {
      MPI_Reduce(&accloc, &acctot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Bcast(&acctot, 1, MPI_INT, 0, MPI_COMM_WORLD);
      gamma *= 2 * (double)acctot / (100 * nx * ny);
      /* reset acceptance count */
      accloc = acctot = 0;
      if (gamma > PI / 2) {
	if (DEBUG == 1)
	  printf("gamma %f - resetting\n", gamma);
	gamma = PI / 2;
      }
      else if (gamma < PI / 1800) {
	gamma = PI / 1800;
	if (DEBUG == 1)
	  printf("gamma below 0.1 degrees\n");
      }
      e = func(grid, lock, 0);
      if (rank == 0 && DEBUG == 1) {
	fprintf(fp, "%lu %f %f %f %f\n", t, e, gamma, p, beta);
      }
    }
    
    for (pass = 1; pass <= 2; pass++) {
      if (pass == 1) {
	srow = 1;
	erow = (float)h / 2;
      }
      else if (pass == 2) {
	srow = erow + 1;
	erow = h;
      }
      for (j = srow; j <= erow; ++j) {
	for (i = 1; i <= nx; ++i) {
	  if (lock[j][i] == 0) {
	    angle = grid[j][i];
	    e = locfunc(grid, i, j, 0);
	    grid[j][i] += ((ran2(&idum)*2)-1)*gamma;
	    e2 = locfunc(grid, i, j, 0);
	    if (e2 - e > 0) {
	      num = ran2(&idum);
	      p = exp((e-e2)*beta);
	      if (num < p) {
		accloc++;
	      }
	      else
		grid[j][i] = angle;
	    }
	    else
	      accloc++;
	  }
	}
      }
      /* communication step */
      /* on pass 0, source is rank - 1, destination is rank + 1 */
      if (pass == 1) {
	if (rank % 2 == 1) {
	  MPI_Send(&grid[1][1], 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	  if (rank < numProcs - 1)
	    MPI_Recv(&grid[h+1][1], 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &stat);
	}
	if (rank % 2 == 0) {
	  if (rank < numProcs - 1)
	    MPI_Recv(&grid[h+1][1], 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &stat);
	  if (rank > 0)
	    MPI_Send(&grid[1][1], 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	}
      }
      if (pass == 2) {
	if (rank % 2 == 1) {
	  if (rank < numProcs - 1)
	    MPI_Send(&grid[h-1][1], 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	  MPI_Recv(&grid[-1][1], 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &stat);
	}
	if (rank % 2 == 0) {
	  if (rank > 0)
	    MPI_Recv(&grid[-1][1], 2 * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &stat);
	  if (rank < numProcs - 1)
	    MPI_Send(&grid[h-1][1], 2 * nx, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	}
      }
    }
  }
  if (rank == 0)
    fclose(fp);
}
