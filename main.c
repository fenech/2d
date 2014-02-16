#include "header.h"
	
int main(int argc, char **argv)
{	
  double fret;            /* Frank energy */
  double start, end;      /* for timing */
  double **grid;          /* 2D grid */
  int **lock;             /* locked cells */
  double **fgrid;         /* fullgrid */
  int **flock;            /* locked cells (fullgrid) */
  char fname[20];
  char rname[20];
  int iter = 0;
  long maxiter;
  int flag, rank;         /* MPI variables */
  float t0;               /* starting "temperature" */
  FILE *fp;            

  MPI_Init(&argc, &argv);
  MPI_Initialized(&flag);
  if (flag != 1)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 8) {
    printf("Usage: <exe> <x length> <y length> <monte carlo steps> <temp> <particle radius> <separation> <boundary angle> [id]\n");
    exit(1);
  }
 
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  maxiter = atoi(argv[3]);
  t0 = atof(argv[4]);
  rad = atoi(argv[5]);
  sep = atoi(argv[6]);
  ba = atoi(argv[7]);

  if (argc == 9)
    id = atoi(argv[8]);
  else 
    id = 1;

  initialise(&grid, &lock);

  fret = func(grid, lock, 0);
  sprintf(fname, "%di%.1f_%d", ny, t0, sep);
  print(grid, fname);
  
  sprintf(rname, "report%d", sep); 
  if (rank == 0)
    printf("Initial Frank Energy: %f\n", fret);
  
  start = MPI_Wtime();
  monte(grid, lock, maxiter, t0);
  end = MPI_Wtime();
  
  fret = func(grid, lock, 0);  

  grid2root(grid, &fgrid);
  lock2root(lock, &flock);
  prncont(fgrid, "testgrid");

  if (rank == 0)
    conjgrad(fgrid, flock, 10000, &fret);

  if (rank == 0) {
    fp = fopen(rname, "w");
    fprintf(fp, "%d %f\n", sep, fret);
    fclose(fp);
    printf("End Frank Energy:     %f\n", fret);
    printf("No. iterations:       %d\n", iter);
    printf("Time taken:           %f\n", end - start);
  }

  MPI_Finalize();
  return 0;
}
