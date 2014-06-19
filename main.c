#include "header.h"
	
int main(int argc, char **argv)
{	
    double fret;            /* Frank energy */
    double start, end;      /* for timing */
    double **grid;          /* 2D grid */
    int **lock;             /* locked cells */
    double **fgrid;         /* fullgrid */
    int **flock;            /* locked cells (fullgrid) */
    char fname[128] = "log";
    char gname[128] = "grid";
    char suffix[128];
    int iter = 0;
    long maxiter;
    int flag, rank, np;     /* MPI variables */
    float t0;               /* starting "temperature" */
    FILE *fp, *log_fp;
    t_par par[2];
    int sep, ba;

    MPI_Init(&argc, &argv);
    MPI_Initialized(&flag);
    if (flag != 1) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (argc < 11 || argc > 12) {
	if (rank == 0) printf("Usage: <%s> <x length> <y length> <monte carlo steps> <temp> <major axis> <minor axis> <align> <theta> <separation> <boundary angle> [id]\n", argv[0]);
	MPI_Finalize();
	exit(1);
    }
 
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    maxiter = atoi(argv[3]);
    t0 = atof(argv[4]);
    par[0].major = par[1].major = atoi(argv[5]);
    par[0].minor = par[1].minor = atoi(argv[6]);
    if (strcmp(argv[7], "para") != 0 && strcmp(argv[7], "perp") != 0) {
	if (rank == 0) printf("Alignment must be para or perp\n");
	MPI_Finalize();	
	exit(1);
    }
    strcpy(par[0].align, argv[7]);
    strcpy(par[1].align, argv[7]);
    par[0].theta = 0;
    par[1].theta = atof(argv[8]);
    sep = atoi(argv[9]);
    ba = atoi(argv[10]);

    if (argc == 12) id = atoi(argv[11]);
    else id = 1;

    par[0].cy = par[1].cy = ny / 2;
    par[0].cx = nx / 2 - par[0].major - sep / 2;
    par[1].cx = nx / 2 + par[1].major + sep / 2 - 1;
    
    sprintf(suffix, "r%dx%d_t%.0f_s%d_a%d_%d_%s", 
	    par[0].major, par[0].minor, par[1].theta, sep, ba, id, par[0].align);
    par[1].theta = PI * par[1].theta / 180.0;

    int success = initialise(&grid, &lock, par, sep, ba, suffix);
    int all_succeeded;
    MPI_Allreduce(&success, &all_succeeded, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (all_succeeded != np) {
	MPI_Finalize();
	return 0;
    }

    strcat(gname, suffix);
    print(grid, gname);

    fret = func(grid, lock, 0);
    if (rank == 0) {
	printf("Initial Frank Energy: %f\n", fret);
	strcat(fname, suffix);	
	log_fp = fopen(fname, "w");
    }

    start = MPI_Wtime();
    monte(grid, lock, maxiter, t0, log_fp, suffix);
    end = MPI_Wtime();

    grid2root(grid, &fgrid);
    lock2root(lock, &flock);
    prncont(fgrid, "testgrid");

    fret = func(grid, lock, 0);
    if (rank == 0) conjgrad(fgrid, flock, 10000, &fret, log_fp, suffix, par, sep);
    
    if (rank == 0) {		
	printf("End Frank Energy:     %f\n", fret);
	printf("No. iterations:       %d\n", iter);
	printf("Time taken:           %f\n", end - start);
    }

    MPI_Finalize();
    return 0;
}
