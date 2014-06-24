#ifndef CONVERT_H
#define CONVERT_H

#include <stdio.h>
#include "libaat.h"

/*#include "halotype.h"*/
/*
#ifndef LATLON_TYPE
typedef struct {
  double lt, lg;
} LatLon;
#define LATLON_TYPE
#endif
*/

typedef struct {
/* The following defines a record in ?basins.dta */
  char code[4];
  char type;
  char name[19];
  LatLon G, P;
  double xig, yjg, gsize; /* Center point and grid size */
  double ximin, ximax, yjmin, yjmax; /* Dimmension of SLOSH Grid 1,xmax 1,ymax */
  double elipty;
/* End definition of record. */

  /* These are w.r.t. local north-east coordinates?? */
  double Xgpt, Ygpt;
  double Xppt, Yppt;
  double Thtg2p, Delthg, Delrg, Rg2p, Abqab;
} bsndta_type;

typedef struct {
  LatLon T;
  double Stlatd, Ctlatd, Tscald;
} tanplane_type;

typedef struct {
  bsndta_type bsn;
  tanplane_type slsh;
} convert_type;

double Fort_atof (char *s, size_t w, size_t d);
double a2lat (char *s);
void cl2cnf (double lat, double *cnf, double *scale);
void cnf2cl (double cnf, double *clk, double *scale);
void sll2xy (LatLon In, double *x, double *y, tanplane_type *slsh);
void sxy2ll (double x, double y, LatLon *Out, tanplane_type *slsh);
void pq2xy (double Xi, double Yj, double *x, double *y, bsndta_type *bsn);
void xy2pq(double x, double y, double *p, double *q, bsndta_type *bsn);
void str2bsndta (char *s, bsndta_type *bsn, tanplane_type *slsh);
int load_bsndta (char *filename, char *abrev, bsndta_type *bsn,
                 tanplane_type *slsh);
void ltln2pq (LatLon DX, double *p, double *q, bsndta_type *bsn,
             tanplane_type *slsh);
void pq2ltln (double p, double q, LatLon *DX, bsndta_type *bsn,
             tanplane_type *slsh);
int InsideBasin (convert_type *ct, LatLon pt, double d);
int InsideBasin2 (convert_type *ct, LatLon pt, double d, int pMin, int pMax,
                  int qMin, int qMax);

#endif
