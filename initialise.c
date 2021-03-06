#include "header.h"

int initialise(double*** grid, int*** lock, const t_par * par, int sep, int ba, const char * suffix, random_key * key)
{
    int i, j, k, l;                            /* counters */
    int rank, numProcs;                        /* MPI variables */
    int r;                                     /* remainder used in grid allocation */
    float angle = 0;                           /* angle at top and bottom boundary */
    int bottom, top;                           /* grid boundaries indices */
    MPI_Status stat;
    double particle;
    FILE *fp = NULL;
    int *hs;
    char fname[128] = "par";
    float bangle = ba * PI / 180;               /* convert ba to radians */
    int ret_val = 1;

    strcat(fname, suffix);

    /*  char fname[10];*/
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    threefry4x64_ukey_t uk={{}};
    size_t len = strlen(fname);
    int seed = 0;
    for (i = 0; i < len; ++i) {
        seed += fname[i];
    }
    uk.v[0] = seed;
    uk.v[1] = rank;

    random_ctr rc = {{}};
    *key = threefry4x64keyinit(uk);

    /* allocate grid and lock arrays */
    h = ny / numProcs;
    r = ny - h * numProcs;
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
        rc.v[0] = j;
        for (i = 1; i <= nx; i += 4) {
            rc.v[1] = i;
            threefry4x64_ctr_t ran = threefry4x64(rc, *key);
            for (int c = 0; c < 4 && i + c <= nx; ++c) {
                (*grid)[j][i+c] = u01fixedpt_open_closed_64_53(ran.v[c]) * PI;
                (*lock)[j][i] = 0;
            }
        }
        /* lock left and right edges */
        // (*grid)[j][1] = (*grid)[j][nx] = bangle;
        // (*lock)[j][1] = (*lock)[j][nx] = 1;
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
    if (rank == 0 || rank == numProcs - 1) {
        for (i = 1; i <= nx; ++i) {
            (*grid)[r][i] = angle;
            (*lock)[r][i] = 1;
        }
    }
    if (rank == 0) {
        fp = fopen(fname, "w");
    }

    for (k = 0; k < 2; ++k) {
        t_par p = par[k];
        double f = sqrt(p.major * p.major - p.minor * p.minor);
        double foc[][2] = {
                { p.cx + f * cos(p.theta), p.cy + f * sin(p.theta) },
                { p.cx - f * cos(p.theta), p.cy - f * sin(p.theta) }
        };
        for (j = p.cy - p.major; j <= p.cy + p.major; ++j) {
            for (i = p.cx - p.major; i <= p.cx + p.major; ++i) {
                double a = (cos(p.theta) * (i - p.cx) + sin(p.theta) * (j - p.cy)) *
                        (cos(p.theta) * (i - p.cx) + sin(p.theta) * (j - p.cy));
                double b = (sin(p.theta) * (i - p.cx) - cos(p.theta) * (j - p.cy)) *
                        (sin(p.theta) * (i - p.cx) - cos(p.theta) * (j - p.cy));
                if (a / (p.major * p.major) + b / (p.minor * p.minor) <= 1.0) {
                    double p2foc[][2] = {
                            { foc[0][0] - i, foc[0][1] - j },
                            { foc[1][0] - i, foc[1][1] - j }
                    };
                    for (l = 0; l < 2; ++l) {
                        double mag = sqrt(p2foc[l][0] * p2foc[l][0] + p2foc[l][1] * p2foc[l][1]);
                        if (mag != 0.0) {
                            p2foc[l][0] /= mag;
                            p2foc[l][1] /= mag;
                        }
                        else {
                            // printf("zero magnitude vector\n");
                            // ret_val = 0;
                        }
                    }
                    particle = (atan2(p2foc[0][1], p2foc[0][0]) + atan2(p2foc[1][1], p2foc[1][0])) / 2.0;
                    if (strcmp(p.align, "para") == 0) particle += PI / 2.0;

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

    return ret_val;

}
