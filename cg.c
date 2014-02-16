#include "header.h"

void conjgrad(double **fgrid, int **lock, int maxits, double *fret)
{
  int i, j, k, l;                     /* counters */ 
  double **g, **h;                    /* direction matrices used in conjugate gradient */
  double **dgrid;                     /* grid of partial derivatives */
  double gamt, gamb;                  /* top and bottom part of fraction in evaluating gamma */
  double gam;                         /* gamma */
  double **gra, **grb, **grc;         /* three grid positions bracketing a minimum */
  double **grt;                       /* trial/dummy grid position */
  double fa, fb, fc, ft;              /* energy values at different grid positions */
  double fprev;                       /* records the energy from the previous iteration */
  double a2b, b2c;                    /* distances gra <-> grb and grb <-> grc */
  int strike = 0;
  double estart;
  FILE *fp;

  int pass;
  double **emap;
  /* FILE *fp2; */
  /*  char dname[20]; */
  char ename[20];
  char lname[30];
  char gname[30];

  fprev = estart = *fret;
  fa = fb = fc = ft = 0.;

  sprintf(gname, "finalr%ds%da%d_%d", rad, sep, ba, id);
  sprintf(lname, "logr%ds%da%d_%d", rad, sep, ba, id);
  fp = fopen(lname, "a");    

  emap = dmatrix(1, ny, 1, nx);
  g = dmatrix(1, ny, 1, nx);
  h = dmatrix(1, ny, 1, nx);
  dgrid = dmatrix(1, ny, 1, nx);
  gra = dmatrix(0, ny + 1, 1, nx);
  grb = dmatrix(0, ny + 1, 1, nx);
  grc = dmatrix(0, ny + 1, 1, nx);
  grt = dmatrix(0, ny + 1, 1, nx);

  for (j = 1; j <= ny; ++j) {         /* copy fgrid to gra */
    for (i = 1; i <= nx; ++i) {
      gra[j][i] = fgrid[j][i];
    }
  }
   
  for (i = 1; i <= nx; ++i) {         /* copy top and bottom rows of fgrid to all grids */
    gra[0][i] = grb[0][i] = grc[0][i] = grt[0][i] = fgrid[0][i];    
    gra[ny+1][i] = grb[ny+1][i] = grc[ny+1][i] = grt[ny+1][i] = fgrid[ny+1][i];
  }

  for (j = 1; j <= ny; ++j)           /* start with arbitrary g0 = h0 */  
    for (i = 1; i <= nx; ++i)
      if (lock[j][i] == 0)
	h[j][i] = g[j][i] = 1;

  for (pass = 1; pass <= 2; ++pass) {
    /* METHOD BEGINS HERE */
    for (k = 1; k <= maxits; ++k) {
      sprintf(ename, "emap%d_%d", pass, k);
      /*
      for (j = 1; j <= ny; ++j)
	for (i = 1; i <= nx; ++i)
	  emap[j][i] = locfunc(gra, i, j, 1);
      prncont(emap, 1, ny, 1, nx, ename);
      */
      if (DEBUG == 1)
	fprintf(fp, "#iteration %d\n", k);
      
      /*
	sprintf(gname, "gra%d", k);
	prn(gra, 1, ny, 1, nx, gname);
      */

      dfunc(gra, lock, dgrid);         /* calculate new g from derivative */

      /*
      sprintf(dname, "dgr%ds%d_%d", rad, sep, k);
      prncont(dgrid, 1, ny, 1, nx, dname);
      */
      gamt = gamb = 0.0;                /* calculate gamma */
      for (j = 1; j <= ny; ++j) {
	for (i = 1; i <= nx; ++i) {
	  gamt += (dgrid[j][i] - g[j][i]) * dgrid[j][i];
	  gamb += g[j][i] * g[j][i];
	}
      }
      gam = gamt / gamb;
      /* 
	 if (gam < EPS) {
	 fprintf(fp, "#gamma = 0. exiting\n");
	 prn(gra, 1, ny, 1, nx, gname);
	 fclose(fp);
	 return;
	 }
      */
      if (DEBUG == 1)
	fprintf(fp, "#gamma: %f\n", gam);
    
      for (j = 1; j <= ny; ++j) {       /* calculate h */
	for (i = 1; i <= nx; ++i) {
	  h[j][i] *= gam;
	  h[j][i] += dgrid[j][i];
	  h[j][i] /= 1;
	}
      }
      if (pass == 2) {
	for (j = 35; j <= 65; ++j) {
	  for (i = nx / 2 - 2 * rad - sep / 2; i <= nx / 2 - sep / 2; ++i)
	    h[j][i] = 0;
	  for (i = nx / 2 + sep / 2; i <= nx / 2 + 2 * rad + sep / 2; ++i)
	    h[j][i] = 0;
	}
      }
      /*
      sprintf(dname, "h%ds%d_%d", rad, sep, k);
      prncont(h, 1, ny, 1, nx, dname);
      */
      /* now minimise in direction of h */
      /* bracket a minimum */
      for (j = 1; j <= ny; ++j)
	for (i = 1; i <= nx; ++i)           
	  grb[j][i] = gra[j][i] + h[j][i];    /* 1st point is gra, grb is gra + h */
    
      /* evaluate energy at first and second points
	 to establish direction of decrease */
      fa = func(gra, lock, 1);         
      fb = func(grb, lock, 1);
        
      if (fb > fa) {
	for (j = 1; j <= ny; ++j)
	  for (i = 1; i <= nx; ++i) {
	    grt[j][i] = gra[j][i];
	    gra[j][i] = grb[j][i];
	    grb[j][i] = grt[j][i];
	  }
	ft = fb; fb = fa; fa = ft;
      }

      for (j = 1; j <= ny; ++j)
	for (i = 1; i <= nx; ++i)
	  grc[j][i] = grb[j][i] + GOLD * (grb[j][i] - gra[j][i]);
    
      fc = func(grc, lock, 1);
      /*
	sprintf(ename, "%de%d_%d", sep, k, index);        
	fp2 = fopen(ename, "w");
	for (l = 0; l <= 20; ++l) {
	for (j = 1; j <= ny; ++j)
	for (i = 1; i <= nx; ++i)
	grt[j][i] = gra[j][i] + l * h[j][i];
	fprintf(fp2, "%d %f\n", l, func(grt, lock, 1));
	}
	fclose(fp2);	
      */
      l = 0;
      if (DEBUG == 1)
	fprintf(fp, "#fa: %f, fb: %f, fc: %f\n", fa, fb, fc); 
      while (fb > fc) {     /* keep moving in h direction until we are out of the valley */
	l++;
	for (j = 1; j <= ny; ++j) 
	  for (i = 1; i <= nx; ++i) {
	    gra[j][i] = grb[j][i];
	    grb[j][i] = grc[j][i];
	    grc[j][i] = grb[j][i] + GOLD * (grb[j][i] - gra[j][i]);
	  }
	fb = func(grb, lock, 1);
	fc = func(grc, lock, 1);      
	if (DEBUG == 1)
	  fprintf(fp, "#hop %d - fa: %f, fb: %f, fc: %f\n", l, fa, fb, fc); 
      }
      if (DEBUG == 1)
	fprintf(fp, "#made %d hops\n", l);

      /* minimum should lie between gra, grb and grc
	 now perform line minimisation */

      a2b = b2c = 1.0;
      while (a2b + b2c > 3.0e-8) {
	/* work out distances gra <-> grb and grb <-> grc */
	a2b = b2c = 0.0;
	for (j = 1; j <= ny; ++j) {
	  for (i = 1; i <= nx; ++i) {
	    a2b += (grb[j][i] - gra[j][i]) * (grb[j][i] - gra[j][i]);
	    b2c += (grc[j][i] - grb[j][i]) * (grc[j][i] - grb[j][i]);
	  }
	}
	a2b = sqrt(a2b);
	b2c = sqrt(b2c);
      
	if (a2b > b2c) {
	  if (DEBUG == 1)
	    fprintf(fp, "#a2b larger, moving from grb towards gra\n");
	  for (j = 1; j <= ny; ++j) {
	    for (i = 1; i <= nx; ++i) {
	      /* if gra <-> grb is larger then move in that direction */
	      grt[j][i] = grb[j][i] + CGOLD * (gra[j][i] - grb[j][i]);
	    }
	  }	
	  ft = func(grt, lock, 1);
	  if (DEBUG == 1)
	    fprintf(fp, "#fa: %f, fb: %f, fc: %f ft: %f\n", fa, fb, fc, ft); 
	  if (ft < fb) {       /* if ft is lower we have a new minimum to bracket around*/
	    if (DEBUG == 1) {
	      fprintf(fp, "#ft lower than fb\n");
	      fprintf(fp, "#grc changed to previous grb\n");
	      fprintf(fp, "#grb changed to grt\n");
	    }
	    for (j = 1; j <= ny; ++j) {
	      for (i = 1; i <= nx; ++i) {
		grc[j][i] = grb[j][i];
		grb[j][i] = grt[j][i];
	      }
	    }
	    fc = fb;
	    fb = ft;
	  }
	  else {               /* if ft is higher */
	    for (j = 1; j <= ny; ++j) {
	      for (i = 1; i <= nx; ++i) {
		/* grt is the new left hand of the bracket */
		gra[j][i] = grt[j][i];
	      }
	    }
	    fa = ft;
	  }
	}	
	
	else 
	  if (DEBUG == 1) 
	    fprintf(fp, "#b2c larger, moving from grb towards grc\n");
	  for (j = 1; j <= ny; ++j) {
	    for (i = 1; i <= nx; ++i) {
	      /* if grb <-> grc is larger then move in that direction */
	      grt[j][i] = grb[j][i] + CGOLD * (grc[j][i] - grb[j][i]);
	    }
	  }
	  ft = func(grt, lock, 1);
	  if (DEBUG == 1)
	    fprintf(fp, "#fa: %f, fb: %f, fc: %f ft: %f\n", fa, fb, fc, ft);
	  if (ft < fb) {       /* if ft is lower we have a new minimum to bracket around*/
	    if (DEBUG == 1) {
	      fprintf(fp, "#ft lower than fb\n");
	      fprintf(fp, "#gra changed to previous grb\n");
	      fprintf(fp, "#grb changed to grt\n");
	    }
	    for (j = 1; j <= ny; ++j) {
	      for (i = 1; i <= nx; ++i) {
		gra[j][i] = grb[j][i];
		grb[j][i] = grt[j][i];
	      }
	    }
	    fa = fb;
	    fb = ft;
	  }
	  else {               /* if ft is higher */
	    for (j = 1; j <= ny; ++j) {
	      for (i = 1; i <= nx; ++i) {
		/* grt is the new right hand */
		grc[j][i] = grt[j][i];
	      }
	    }
	    fc = ft;
	  }
      }
      
    
    
      /* line minimisation complete, set gra to middle of the bracket */
      for (j = 1; j <= ny; ++j)
	for (i = 1; i <= nx; ++i) {
	  gra[j][i] = grb[j][i];
	  g[j][i] = dgrid[j][i];
	}
      
      if (DEBUG == 1)
	fprintf(fp, "#_________________\n");

      if(fabs(fb - fprev) <= EPS)
	strike++;
      else
	strike = 0;
    
      fprev = *fret = fb;

      if(strike == 50) {
	fprintf(fp, "#exiting after %d iterations\n", k);
	fprintf(fp, "#final energy:\n # %f\n", fb);
	prn(gra, 0, ny+1, 1, nx, gname);
	break;
      }
    }
  }
  fprintf(fp, "#final energy:\n # %f\n", fb);
  fclose(fp);
  prn(gra, 0, ny+1, 1, nx, gname);
}

