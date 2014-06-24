/*****************************************************************************
 ** <convert.c> :: Ansi C
 **   1/98  Arthur Taylor(AAT)     TDL
 **     Tested on MS-Windows 3.1 platform
 **
 ** Purpose:
 **     This module contains the code needed to convert from lat/lon on the
 **   Clarke spheroid to p/q co-ordinates internal to a SLOSH basin.  It also
 **   contains the necessary code to convert from p/q co-ordinates to lat/lon
 **   on the Clarke spheroid.
 **
 ** Files Needed:
 **   Source: <ltln2pq.c> <ltln2pq.h> <complex.h> <complex.h>
 **
 ** Global Variables: (S = Static/Private) (E = Extern/Public)
 **
 ** History:
 **   Original code in Jye Chen's ltln2pq.for and his pq2ltln.for
 **     Jye's code extracted parts of the ADTaylor map transformations, used
 **     to produce data for the SLOSH program, and also in the SLOSH program
 **     itself.
 **   1/28/98  AAT:  Finished debugging ltln->pq, added pq->ltln
 **   1/29/98  AAT:  Finished debugging pq->ltln, started commenting.
 ** Bugs/Desired Features:
 ****************************************************************************/
#include "convert.h"
/*
#include "haloinit.h"
*/
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* RADPDG is radians per degree. */
/* ECCEN is the Ecentricity of the Earth. */
/* A_RAD is the semi-major axis of the Earth in meters */
/* B_RAD is the semi-minor axis of the Earth in meters */
/* METPSM is meters per statute mile */
#define RADPDG 0.01745329251994328
#define ECCEN 0.082437112686627
#define A_RAD 6.378294E6
/* #define B_RAD 6.356584E6 */
#define METPSM 1609.344

/*****************************************************************************
 * <Fort_atof> :: Arthur Taylor TDL
 * Purpose:
 *     Converts a string to a double using FORTRAN conventions.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   s           (I) The string to extract the double from
 *   w           (I) The FORTRAN field width.
 *   d           (I) The FORTRAN decimal field.
 *
 * Returns: A double taken from s using FORTRAN conventions.
 * History: 1/29/98 AAT Commented.
 * Notes:
 *   Apparently fortran does not put zeros in the blanks if it is reading an
 *     integer into a F4.1 statement.  Also F4.5 is allowed.
 *     (I checked both lahey and microsoft DOS fortran compilers)
 ****************************************************************************/
double Fort_atof (char *s, size_t w, size_t d)
{
   char c;
   double ans;

   if (strlen (s) > w) {
      c = s[w];
      s[w] = '\0';
      if (strchr (s, '.') != NULL)
         ans = atof (s);
      else
         ans = atof (s) / pow (10, d);
      s[w] = c;
      return ans;
   } else {
      if (strchr (s, '.') != NULL)
         ans = atof (s);
      else
         ans = atof (s) / pow (10, d);
   }
   return ans;
}

/*****************************************************************************
 * <a2lat> :: Arthur Taylor TDL
 * Purpose:
 *     Converts a string containing deg min sec for latitude to a double.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   s           (I) The string to get the latitude from.
 *
 * Returns: A double with the latitude.
 * History: 1/29/98 AAT Commented.
 *   1/6/99 AAT Fixed to handle lack of 0's better. (and made a2lat == a2lon)
 * Notes:
 ****************************************************************************/
double a2lat (char *s)
{
   int deg, min, sec;
   double ans;
   char buffer[20], c, *tmp;

   c = s[12];
   s[12] = '\0';
   ans = atof (s);
   s[12] = c;
   sprintf (buffer, "%+09.5f", ans);

   tmp = strchr (buffer, '.');
   *tmp = '\0';
   deg = atoi (buffer);
   *tmp = '.';
   tmp++;
   c = *(tmp + 2);
   *(tmp + 2) = '\0';
   min = atoi (tmp);
   *(tmp + 2) = c;
   tmp += 2;
   c = *(tmp + 3);
   *(tmp + 3) = '\0';
   sec = atoi (tmp);
   *(tmp + 3) = c;

   if (deg < 0) {
      ans = deg - min / 60.0 - sec / (3600.0 * 10.0);
   } else {
      ans = deg + min / 60.0 + sec / (3600.0 * 10.0);
   }
   return (ans);
}

/*****************************************************************************
 * <cl2cnf> :: Arthur Taylor TDL
 * Purpose:
 *     Convert a latitude on the clarke sphere (regular earth latitude) to one
 *   on the conformal sphere (a perfect sphere gotten from the clarke sphere
 *   by pulling on the poles).
 *     Also returns local scale of the map between spheroids.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clk         (I) The clarke latitude.
 *   cnf         (O) The conformal latitude.
 *   scale       (O) The local scale of the map between spheroids.
 *
 * Returns: NULL
 * History: 1/29/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void cl2cnf (double clk, double *cnf, double *scale)
{
   double sin_clk, temp0, temp1, temp2;

   sin_clk = sin (clk * RADPDG);
   temp0 = ECCEN * sin_clk;
   temp1 = (1. - sin_clk) * pow ((1. + temp0), ECCEN);
   temp2 = (1. + sin_clk) * pow ((1. - temp0), ECCEN);
   *cnf = atan2 (temp2 - temp1, 2. * sqrt (temp1 * temp2)) / RADPDG;
   *scale =
      pow ((1 - temp0 * temp0), ((1 + ECCEN) * .5)) * 2. / (temp1 + temp2);
}

/*****************************************************************************
 * <cnf2cl> :: Arthur Taylor TDL
 * Purpose:
 *     Convert a latitude on the conformal sphere (a perfect sphere gotten
 *   from the clarke sphere by pulling on the poles) to the clarke sphere
 *   (regular earth latitude).
 *     Also returns local scale of the map between spheroids.  The scale
 *   should be a ratio of distance on conformal sphere to clarke sphere, and
 *   should be >= 1.0.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   cnf         (I) The conformal latitude.
 *   clk         (O) The clarke latitude.
 *   scale       (O) The local scale of the map between spheroids.
 *
 * Returns: NULL
 * History: 1/29/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void cnf2cl (double cnf, double *clk, double *scale)
{
   double sin_clk, temp0, temp1, temp2, temp_lat, grad;
   size_t i;

   *clk = cnf;
   for (i = 0; i < 3; i++) {
      sin_clk = sin (RADPDG * (*clk));
      temp0 = ECCEN * sin_clk;
      temp1 = (1. - sin_clk) * pow ((1. + temp0), ECCEN);
      temp2 = (1. + sin_clk) * pow ((1. - temp0), ECCEN);
      temp_lat = atan2 (temp2 - temp1, 2 * sqrt (temp1 * temp2)) / RADPDG;
      grad = 2. * (1. - ECCEN * ECCEN) / (temp1 + temp2) /
         pow ((1. - temp0 * temp0), (1. - ECCEN / 2.0));
      (*clk) = (*clk) + (cnf - temp_lat) / grad;
   }
   sin_clk = sin (RADPDG * (*clk));
   temp0 = ECCEN * sin_clk;
   temp1 = (1. - sin_clk) * pow ((1. + temp0), ECCEN);
   temp2 = (1. + sin_clk) * pow ((1. - temp0), ECCEN);
   *scale =
      pow ((1. - temp0 * temp0),
           ((1. + ECCEN) * 0.5)) * 2. / (temp1 + temp2);
}

/*****************************************************************************
 * <sll2xy> :: Arthur Taylor TDL
 * Purpose:
 *     Convert lat/lon (on Conformal sphere) to x,y on SLOSH grid where x,y
 *   are in statute miles, with the origin at the tangent point, y-axis to
 *   local north.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   In          (I) The lat/lon input pair
 *   x,y         (O) The result on the SLOSH grid.
 *   slsh        (I) The parameters defining the tangent plane.
 *
 * Returns: NULL
 * History: 1/29/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void sll2xy (LatLon In, double *x, double *y, tanplane_type * slsh)
{
   double sin_lat, cos_lat, cos_del_lon;
   double del_lon, t1, t2, t3, t4;

   sin_lat = sin (In.lat * RADPDG);
   cos_lat = cos (In.lat * RADPDG);
   del_lon = In.lon - slsh->T.lon;
   cos_del_lon = cos (del_lon * RADPDG);
   t1 = 1. + sin_lat * slsh->Stlatd + cos_lat * slsh->Ctlatd * cos_del_lon;
   t2 = -1 * cos_lat * sin (del_lon * RADPDG);
   t3 = sin_lat * slsh->Ctlatd - cos_lat * slsh->Stlatd * cos_del_lon;
   t4 = ((A_RAD / METPSM) * (2. / t1) / slsh->Tscald);
   *x = t4 * t2;
   *y = t4 * t3;
}

/*****************************************************************************
 * <sxy2ll> :: Arthur Taylor TDL
 * Purpose:
 *     Convert x,y on the SLOSH grid to lat/lon on Conformal Sphere.
 *   x,y coordinates are in statute miles, with the origin at the tangent
 *   point, y-axis to local north.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   x,y         (I) The x,y coordinates on the SLOSH grid.
 *   Out         (O) The resulting lat/lon pair.
 *   slsh        (I) The parameters defining the tangent plane.
 *
 * Returns: NULL
 * History: 1/29/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void sxy2ll (double x, double y, LatLon * Out, tanplane_type * slsh)
{
   double xp, yp, r, temp, temp_2;
   double sin_psi, cos_psi, sin_nu, cos_nu;

   xp = x * slsh->Tscald;
   yp = y * slsh->Tscald;
   r = sqrt (xp * xp + yp * yp);
   if (r > 0) {
      sin_nu = xp / r;
      cos_nu = yp / r;
   } else {
      sin_nu = 0;
      cos_nu = 1;
   }

   temp = 0.5 * r * METPSM / A_RAD;
   temp_2 = temp * temp;
   sin_psi = 2. * temp / (1. + temp_2);
   cos_psi = (1. - temp_2) / (1. + temp_2);

   xp = cos_psi * slsh->Ctlatd - sin_psi * slsh->Stlatd * cos_nu;
   yp = -1 * sin_psi * sin_nu;
   r = sqrt (xp * xp + yp * yp);
   if (r > 0) {
      Out->lon = slsh->T.lon + atan2 (yp, xp) / RADPDG;
   } else {
      Out->lon = slsh->T.lon;
   }

   temp = cos_psi * slsh->Stlatd + sin_psi * slsh->Ctlatd * cos_nu;
   Out->lat = atan2 (temp, r) / RADPDG;
}

/*****************************************************************************
 * <pq2xy> :: Arthur Taylor TDL
 * Purpose:
 *     Converts SLOSH pq coordinates to xy co-ordinates on the SLOSH tangent
 *   plane.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   p,q         (I) The p,q coordinates in the slosh grid.
 *   x,y         (O) The x,y coordinates on the SLOSH tangent plane.
 *   bsn         (I) The defining parameters of the basin.
 *
 * Returns: NULL
 * History: 1/30/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void pq2xy (double p, double q, double *x, double *y, bsndta_type * bsn)
{
   myComplex Z1, Z2;
   double xa, ya, x1, y1, cs, ss;

   xa = (p - bsn->xig) * bsn->Delrg;
   ya = (bsn->yjg - q) * bsn->Delrg;
   Z1 = myCexp(myCset(xa, ya));
   Z2 = myCmul_Real(myCadd(Z1, myCmul_Real(myCinv(Z1), bsn->Abqab)),
                    bsn->Rg2p);
   x1 = myCimag(Z2);
   y1 = myCreal(Z2);

   cs = cos (bsn->Thtg2p * RADPDG);
   ss = sin (bsn->Thtg2p * RADPDG);
   *x = bsn->Xppt + x1 * cs + y1 * ss;
   *y = bsn->Yppt - x1 * ss + y1 * cs;
}

/*****************************************************************************
 * <xy2pq> :: Arthur Taylor TDL
 * Purpose:
 *     Converts xy co-ordinates on the SLOSH tangent plane to SLOSH pq
 *   co-ordinates.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   p,q         (I) The p,q coordinates in the slosh grid.
 *   x,y         (O) The x,y coordinates on the SLOSH tangent plane.
 *   bsn         (I) The defining parameters of the basin.
 *
 * Returns: NULL
 * History: 1/30/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void xy2pq (double x, double y, double *p, double *q, bsndta_type * bsn)
{
   myComplex Z1, Z2, Z3;
   double xa, ya, cs, ss;

/* Convert p-pt to g-pt line from local north to y azis in x-y plane, east
 *   to x-axis (Mathematical Coordinates)
 *   using fmod (bsn->Thtg2p +270., 360.)
 */

   cs = cos (fmod (bsn->Thtg2p + 270., 360.) * RADPDG);
   ss = sin (fmod (bsn->Thtg2p + 270., 360.) * RADPDG);

/* rotate x-y axis */

   x = x - bsn->Xppt;
   y = y - bsn->Yppt;
   xa = x * cs - y * ss;
   ya = x * ss + y * cs;

   Z1 = myCmul_Real(myCset(xa, ya), 1. / bsn->Rg2p);
   if (bsn->Abqab != 0) {
      Z2 = myCsqrt(myCadd(myCmul(Z1, Z1), myCset(-4 * bsn->Abqab, 0)));
      if ((bsn->Abqab > 0) && (xa < 0)) {
         Z2 = myCmul(myCset(-1, 0), Z2);
      }
      Z1 = myCmul_Real(myCadd(Z1, Z2), 0.5);
   }
/*
   Z3 = myCmul_Real(myClog(Z1), 1./bsn->Delrg);
   *p = myCreal(Z3) + bsn->xig;
   *q = myCimag(Z3) + bsn->yjg;
*/
   Z3 = myClog(Z1);
   *p = myCreal(Z3) / bsn->Delrg + bsn->xig;
   *q = fmod(myCimag(Z3) + (bsn->yjg - 1) * bsn->Delrg + 2 * M_PI, 2 * M_PI)
         / bsn->Delrg + 1;
}

/*****************************************************************************
 * <str2bsndta> :: Arthur Taylor TDL
 * Purpose:
 *     Converts a NULL terminated string containing a line from the ?basins.dta
 *   file into the bsndta type.  Also converts the read in lat/lon from clarke
 *   spheroid to conformal sphere, and initializes the tangent plane.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   s           (I) The NULL terminated string (assumed lowercase)
 *   bsn         (O) The resulting bsndta type.
 *   slsh        (O) Holds the math for the SLOSH tangent plane.
 *
 * Returns: NULL
 * History: 1/29/98 AAT Commented.
 * Notes:
 *   Fortran Input line is:
 *   A3,A1,A18,A9,1X,A10,2X,A9,1X,A10,1X,F5.1,1X,F5.1,1X,F6.4,4F5.1,20X,F6.5
 ****************************************************************************/
void str2bsndta (char *s, bsndta_type * bsn, tanplane_type * slsh)
{
   double temp;

/* Parsing the character string and generating the bsndta_type. */
   strncpy (bsn->code, s + 0, 3);
   bsn->code[3] = '\0';
   bsn->type = (char) tolower (s[3]);
   if (bsn->type == ' ')
      bsn->type = 'p';
   strncpy (bsn->name, s + 4, 18);
   bsn->name[18] = '\0';
   bsn->G.lat = a2lat (s + 22); /* 1x */
   bsn->G.lon = a2lat (s + 32); /* 2x */
   bsn->P.lat = a2lat (s + 44); /* 1x */
   bsn->P.lon = a2lat (s + 54); /* 1x */
   bsn->xig = Fort_atof (s + 65, 5, 1); /* 1x */
   bsn->yjg = Fort_atof (s + 71, 5, 1); /* 1x */
   bsn->gsize = Fort_atof (s + 77, 6, 4);
   bsn->ximin = Fort_atof (s + 83, 5, 1);
   bsn->ximax = Fort_atof (s + 88, 5, 1);
   bsn->yjmin = Fort_atof (s + 93, 5, 1);
   bsn->yjmax = Fort_atof (s + 98, 5, 1); /* 20x */
   if (bsn->type == 'h')
      bsn->elipty = Fort_atof (s + 123, 6, 5);

/* Initializing the tangent plane. */
   cl2cnf (bsn->P.lat, &bsn->P.lat, &temp);
   cl2cnf (bsn->G.lat, &bsn->G.lat, &slsh->Tscald);
   slsh->T = bsn->G;
   slsh->Stlatd = sin (RADPDG * slsh->T.lat);
   slsh->Ctlatd = cos (RADPDG * slsh->T.lat);
   bsn->Xgpt = 0.;
   bsn->Ygpt = 0.;
   sll2xy (bsn->P, &bsn->Xppt, &bsn->Yppt, slsh);
}

/*****************************************************************************
 * <load_bsndta> :: Arthur Taylor TDL
 * Purpose:
 *     Loads the bsndta, and initializes several mathematical constants needed
 *   by the pq2ltlg, and ltlg2pq functions.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   filename    (I) The name of the file to read from (ebasins.dta,
 *                   basins.dta, or hbasins.dta)
 *   abrev       (I) The 3 letter abreviation for the basin.
 *   bsn         (O) The resulting bsndta type.
 *   slsh        (O) Holds the math for the SLOSH tangent plane.
 *
 * Returns: -1 if invalid file, or unable to find the basin, 0 otherwise
 * History: 1/30/98 AAT Commented.
 * Notes:
 *   Calling procedure should verrify that bsn->type == 'e/p/h', and verify
 *   the other fields.
 ****************************************************************************/
int load_bsndta (char *filename, char *abrev, bsndta_type * bsn,
                 tanplane_type * slsh)
{
   FILE *fp;
   int c1;
   char f_found = 0, buffer[1000];
   double temp1, temp2;

   strToLower (abrev);
   if ((fp = fopen (filename, "rt")) == NULL) {
      printf ("Unable to open %s\n", filename);
      return (-1);
   }
   while (((c1 = fgetc (fp)) != EOF) && (f_found == 0)) {
      ungetc (c1, fp);
      fgets (buffer, 1000, fp);
      buffer[strlen (buffer) - 1] = '\0';
      strToLower (buffer);
      if (strncmp (abrev, buffer, 3) == 0)
         f_found = 1;
   }
   fclose (fp);
   if (f_found == 0) {
      printf ("Unable to find Basin %s\n", abrev);
      return (-1);
   }
   str2bsndta (buffer, bsn, slsh);

/* Compute the basic SLOSH basin parameters. */
   temp1 = bsn->Xppt - bsn->Xgpt;
   temp2 = bsn->Yppt - bsn->Ygpt;
   bsn->Rg2p = sqrt (temp1 * temp1 + temp2 * temp2);
   bsn->Thtg2p = 180. + atan2 (temp1, temp2) / RADPDG;
   if (bsn->type == 'e') {
      bsn->Abqab = -1;
      bsn->Delrg =
         RADPDG * (90. - bsn->gsize) / (((int) (bsn->yjg + 0.001)) - 1.);
      bsn->Rg2p = .5 * bsn->Rg2p / cos (bsn->gsize * RADPDG);
/* For e-case only-- reset math origin as center point.  X-axis changed to
 *   be perpendicular to line of P-point to G-point.
 */
      bsn->Thtg2p = bsn->Thtg2p + 90.;
      bsn->Xppt = 0.;
      bsn->Yppt = 0.;
   } else if (bsn->type == 'h') {
      bsn->Abqab = (1. - bsn->elipty) / (1. + bsn->elipty);
      bsn->Rg2p = bsn->Rg2p / (1. + bsn->Abqab);
      if (bsn->name[0] == '+')
         bsn->Delrg = 360. * RADPDG / (bsn->yjmax - bsn->yjmin);
      else
         bsn->Delrg = bsn->gsize / bsn->Rg2p;
   } else {
      bsn->Abqab = 0;
      temp1 = bsn->gsize / (2. * bsn->Rg2p);
      bsn->Delrg = 2. * log (sqrt (1.0 + temp1 * temp1) + temp1);
   }
   bsn->Delthg = bsn->Delrg / RADPDG;
   return (0);
}

/*****************************************************************************
 * <ltln2pq> :: Arthur Taylor TDL
 * Purpose:
 *     Convert a lat/lon to a SLOSH Grid pq.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   DX          (I) The lat/lon to convert.
 *   p,q         (O) The resulting p/q value.
 *   bsn         (I) Holds several mathematical constants for this basin.
 *   slsh        (I) Holds the math for the SLOSH tangent plane.
 *
 * Returns: NULL
 * History: 1/30/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void ltln2pq (LatLon DX, double *p, double *q, bsndta_type * bsn,
              tanplane_type * slsh)
{
   double temp, x1, y1;

   cl2cnf (DX.lat, &DX.lat, &temp);
   sll2xy (DX, &x1, &y1, slsh);
   xy2pq (x1, y1, p, q, bsn);
}

/*****************************************************************************
 * <pq2ltln> :: Arthur Taylor TDL
 * Purpose:
 *     Convert a SLOSH Grid pq to lat/lon.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   p,q         (I) The SLOSH grid p/q value.
 *   DX          (O) The resulting lat/lon.
 *   bsn         (I) Holds several mathematical constants for this basin.
 *   slsh        (I) Holds the math for the SLOSH tangent plane.
 *
 * Returns: NULL
 * History: 1/30/98 AAT Commented.
 * Notes:
 ****************************************************************************/
void pq2ltln (double p, double q, LatLon * DX, bsndta_type * bsn,
              tanplane_type * slsh)
{
   double temp, x1, y1;

   pq2xy (p, q, &x1, &y1, bsn);
   sxy2ll (x1, y1, DX, slsh);
   cnf2cl (DX->lat, &DX->lat, &temp);
}

/*****************************************************************************
 * <convertLoadCmd> :: Arthur Taylor TDL  (halo_convert_Load)
 * Purpose:
 *     To Load a basin's mathematical formulas for converting from lat/lon to
 *   p/q and vice versa.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     filename  (I) The filename for this particular basin.  (ebasins.dta,
 *                   basins.dta, or hbasins.dta)
 *     abrev     (I) 3 letter abreviation for the basin (bos-boston...)
 *     type      (I) e/p/h - optional. (make sure basin is of corect type
 *                   even though it was in the file.
 *
 * Returns: TCL_ERROR or TCL_OK
 * History: 2/2/98 AAT Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int convertLoadCmd (ClientData clientData, Tcl_Interp * interp,
                           int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;

   if ((argc != 3) && (argc != 4)) {
      sprintf (interp->result, "usage: %.50s <filename> <abrev>  <type "
               "(optional)>", argv[0]);
      return TCL_ERROR;
   }
   if (load_bsndta (argv[1], argv[2], &(ct->bsn), &(ct->slsh)) != 0) {
      sprintf (interp->result, "Could not open %s or could not find %s",
               argv[1], argv[2]);
      return TCL_ERROR;
   }
   if (argc == 4) {
      if (ct->bsn.type != tolower (argv[3][0])) {
         sprintf (interp->result,
                  "ERROR: .dta file has different (e/p/h) than it"
                  " should for this basin.\n");
         return TCL_ERROR;
      }
   }
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <convertSaveLLXCmd> :: Arthur Taylor TDL  (halo_conSave_llxSet)
 * Purpose:
 *     To Load a basin's mathematical formulas for converting from lat/lon to
 *   p/q and vice versa, and to save to a .llx file
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     filename  (I) The filename for this particular basin.  (ebasins.dta,
 *                   basins.dta, or hbasins.dta)
 *     abrev     (I) 3 letter abreviation for the basin (bos-boston...)
 *     type      (I) e/p/h (make sure basin is of corect type even though it
 *                   was in the file and filename).
 *     filename2 (I) Where to save it to.
 *     Header    (I) Optional.
 *
 * Returns: TCL_ERROR or TCL_OK
 * History: 4/1999 Arthur Taylor (RSIS/TDL) Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int convertSaveLLXCmd (ClientData clientData, Tcl_Interp * interp,
                              int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   sInt4 imxb, jmxb;
   LatLon DX;
   float temp;
   int i, j;
   char header[161];
   TIO_type *tp;
   short int FID = 1;

   if ((argc != 5) && (argc != 6)) {
      sprintf (interp->result, "usage: %.50s <filename> <abrev> <type> "
               "<filename2> <header>", argv[0]);
      return TCL_ERROR;
   }
   if (load_bsndta (argv[1], argv[2], &(ct->bsn), &(ct->slsh)) != 0) {
      sprintf (interp->result, "Could not open %s or could not find %s",
               argv[1], argv[2]);
      return TCL_ERROR;
   }
   if (ct->bsn.type != tolower (argv[3][0])) {
      sprintf (interp->result,
               "ERROR: .dta file has different (e/p/h) than it"
               " should for this basin.\n");
      return TCL_ERROR;
   }

   if ((tp = tOpen (FID, argv[4], TFLAG_WRITE, TFLAG_MadeOnIntel)) == NULL) {
      sprintf (interp->result, "Difficulties opening %s\n", argv[4]);
      tClose (tp);
      return TCL_ERROR;
   }
   imxb = ct->bsn.ximax + 1;
   jmxb = ct->bsn.yjmax + 1;
   tWrite (&imxb, sizeof (sInt4), 1, tp);
   tWrite (&jmxb, sizeof (sInt4), 1, tp);
   header[0] = '\0';
   if (argc == 6) {
      strncpy (header, argv[5], 160);
   }
   for (i = strlen (header); i < 160; i++) {
      header[i] = ' ';
   }
   header[160] = '\0';
   tWrite (header, sizeof (char), 160, tp);

   for (j = 0; j < jmxb; j++) {
      for (i = 0; i < imxb; i++) {
         /* i+1, j+1 to convert to [1..imxb] system */
         pq2ltln (i + 1, j + 1, &DX, &(ct->bsn), &(ct->slsh));
         temp = DX.lat;
         tWrite (&(temp), sizeof (float), 1, tp);
      }
   }
   for (j = 0; j < jmxb; j++) {
      for (i = 0; i < imxb; i++) {
         /* i+1, j+1 to convert to [1..imxb] system */
         pq2ltln (i + 1, j + 1, &DX, &(ct->bsn), &(ct->slsh));
         temp = -1 * DX.lon;
         tWrite (&(temp), sizeof (float), 1, tp);
      }
   }
   tClose (tp);
   return TCL_OK;
}
#endif

/* return's the current loaded basin's p/q dimmensions. */
#ifdef INCLUDE_TCL
static int convertBasinDimCmd (ClientData clientData, Tcl_Interp * interp,
                               int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   sInt4 imxb, jmxb;

   if (argc != 1) {
      sprintf (interp->result, "usage: %.50s", argv[0]);
      return TCL_ERROR;
   }
   imxb = ct->bsn.ximax + 1;
   jmxb = ct->bsn.yjmax + 1;
   sprintf (interp->result, "%ld %ld", imxb, jmxb);
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <convertLoadSetCmd> :: Arthur Taylor TDL  (halo_conLoad_llxSet)
 * Purpose:
 *     To Load a basin's mathematical formulas for converting from lat/lon to
 *   p/q and vice versa, and to init the latlon grid, thus eliminating need
 *   for .llx files.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     filename  (I) The filename for this particular basin.  (ebasins.dta,
 *                   basins.dta, or hbasins.dta)
 *     abrev     (I) 3 letter abreviation for the basin (bos-boston...)
 *     type      (I) e/p/h (make sure basin is of corect type even though it
 *                   was in the file and filename).
 *     imin, imax (I) Min and max i values to display (compute entire region,
 *                   but only allow [imin...imax] to influence the lat/lon.
 *     jmin, jmax (I) Min and max j values to display (compute entire region,
 *                   but only allow [jmin...jmax] to influence the lat/lon.
 *
 * Returns: TCL_ERROR or TCL_OK
 *          Gives Tcl/Tk: The bounding box of the basin. (In the co-ord system
 *          that I am using in Tcl/Tk (ie positive west longitutde))
 * History: 5/4/98 AAT Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int convertLoadSetCmd (ClientData clientData, Tcl_Interp * interp,
                              int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   Tcl_CmdInfo info;
   basin_type *bt;
   sInt4 imxb, jmxb;
   LatLon DX, point;
   LatLon up, lw;
   int i, j;
   int imin, imax, jmin, jmax;

   if ((argc != 8)) {
      sprintf (interp->result, "usage: %.50s <filename> <abrev>  <type> "
               "<imin> <imax> <jmin> <jmax>", argv[0]);
      return TCL_ERROR;
   }
   if (load_bsndta (argv[1], argv[2], &(ct->bsn), &(ct->slsh)) != 0) {
      sprintf (interp->result, "Could not open %s or could not find %s",
               argv[1], argv[2]);
      return TCL_ERROR;
   }
   if (ct->bsn.type != tolower (argv[3][0])) {
      sprintf (interp->result,
               "ERROR: .dta file has different (e/p/h) than it"
               " should for this basin.\n");
      return TCL_ERROR;
   }
   imin = atoi (argv[4]);
   imax = atoi (argv[5]);
   jmin = atoi (argv[6]);
   jmax = atoi (argv[7]);
   Tcl_GetCommandInfo (interp, "halo_bsnLoadSurge", &info);
   if (info.isNativeObjectProc) {
      sprintf (interp->result,
               "ERROR: halo_bsnLoadSurge should not be an Obj Proc\n");
      return TCL_ERROR;
   } else {
      bt = (basin_type *) info.clientData;

      imxb = ct->bsn.ximax;
      jmxb = ct->bsn.yjmax;
      /* Check if we need to free memory so we can reallocate it. */
      /* Frees memory if i!=0, i!= imxb or if j!=0, j!=jmxb */
      if (((bt->i != 0) && (bt->i != imxb))
          || ((bt->j != 0) && (bt->j != jmxb))) {
         for (i = 0; i < bt->i; i++)
            free (bt->grid[i]);
         free (bt->grid);
         bt->i = 0;
         bt->j = 0;
      }
      /* Allocate memory if needed. */
      if ((bt->i == 0) && (bt->j == 0)) {
         bt->i = imxb;
         bt->j = jmxb;
         bt->val_i = imxb;
         bt->val_j = jmxb;
         bt->start_i = 0;
         bt->start_j = 0;
         bt->stop_i = imxb;
         bt->stop_j = jmxb;
         bt->grid =
            (basingrid_type **) malloc (bt->i * sizeof (basingrid_type *));
         for (i = 0; i < bt->i; i++) {
            bt->grid[i] =
               (basingrid_type *) malloc (bt->j * sizeof (basingrid_type));
         }
      }
      /* Reset variables (sometimes 2 basins may be same size so we don't
       * have to reallocate memory, but we do have to reset the variables. */
      if (bt->fileName != NULL) {
         free (bt->fileName);
         bt->fileName = NULL;
      }
      if (bt->fileName2 != NULL) {
         free (bt->fileName2);
         bt->fileName2 = NULL;
      }
      if (bt->num_cuts != 0) {
         bt->num_cuts = 0;
         free (bt->cuts);
      }
      for (i = 0; i < bt->i; i++) {
         for (j = 0; j < bt->j; j++) {
            /* init depth to land height */
            bt->grid[i][j].depth = 99.9;
            bt->grid[i][j].pen = -1;
            bt->grid[i][j].point.x = -1;

/* Init DD3 Data */
            bt->grid[i][j].ground = 999;
            bt->grid[i][j].momentum = 999;
            bt->grid[i][j].dry = 0;
            bt->grid[i][j].tree = 0;
            bt->grid[i][j].J_cut = -1;
            bt->grid[i][j].I_cut = -1;
            bt->grid[i][j].chnlWidth = 0;
            bt->grid[i][j].BankHeight = 0;

            /* i+1, j+1 to convert to [1..imxb] system */
            pq2ltln (i + 1, j + 1, &DX, &(ct->bsn), &(ct->slsh));
            bt->grid[i][j].pt.lat = DX.lat;
            bt->grid[i][j].pt.lon = -DX.lon;
            /* Convert basin to mercator, and find the bounding box of the
             * basin. */
            /* Also adjusts the longitude. */
            point = bt->grid[i][j].pt;
/*       llx2llm_orig (&point); */
/*        llx2llm_orig (&(bt->grid[i][j].pt)); */
            if ((i == imin) && (j == jmin)) {
               up = point;
               lw = point;
            } else if ((i >= imin) && (i <= imax) && (j >= jmin)
                       && (j <= jmax)) {
               if (point.lat > up.lat)
                  up.lat = point.lat;
               else if (point.lat < lw.lat)
                  lw.lat = point.lat;
               if (point.lon < up.lon)
                  up.lon = point.lon;
               else if (point.lon > lw.lon)
                  lw.lon = point.lon;
            }
         }
      }
      /* bt->up has no default value? */
      bt->scrTop.lat = -1;
      bt->scrTop.lon = -1;
      bt->scrBot.lat = -1;
      bt->scrBot.lon = -1;
/*
    bt->rat_x = -1;
    bt->rat_y = -1;
*/
/*    sprintf(interp->result, "%f %f %f %f", 360-up.lg, up.lt, 360-lw.lg, lw.lt); */
      sprintf (interp->result, "%f %f %f %f", -up.lon, up.lat, -lw.lon,
               lw.lat);
      return TCL_OK;
   }
}
#endif

/*****************************************************************************
 * <convertBoundCmd> :: Arthur Taylor TDL  (halo_conBound)
 * Purpose:
 *     Assuming the basin's mathematical formulas are initialized
 *   (see halo_convert_load), this returns the lat/lon bounds of the basin.
 *   returns clarke lat/lons.
 *****************************************************************************/
#ifdef INCLUDE_TCL
static int convertBoundCmd (ClientData clientData, Tcl_Interp * interp,
                            int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   LatLon DX, up, lw;
   int i, j;

   if ((argc != 1)) {
      sprintf (interp->result, "usage: %.50s", argv[0]);
      return TCL_ERROR;
   }
   for (i = 0; i < ct->bsn.ximax; i++) {
      for (j = 0; j < ct->bsn.yjmax; j += ct->bsn.yjmax - 1) {
         /* i+1, j+1 to convert to [1..imxb] system */
         pq2ltln (i + 1, j + 1, &DX, &(ct->bsn), &(ct->slsh));
         if ((i == 0) && (j == 0)) {
            up = DX;
            lw = DX;
         } else {
            if (DX.lat > up.lat)
               up.lat = DX.lat;
            else if (DX.lat < lw.lat)
               lw.lat = DX.lat;
            if (DX.lon > up.lon)
               up.lon = DX.lon;
            else if (DX.lon < lw.lon)
               lw.lon = DX.lon;
         }
      }
   }
   for (i = 0; i < ct->bsn.ximax; i += ct->bsn.ximax - 1) {
      for (j = 1; j < ct->bsn.yjmax - 1; j++) {
         /* i+1, j+1 to convert to [1..imxb] system */
         pq2ltln (i + 1, j + 1, &DX, &(ct->bsn), &(ct->slsh));
         if (DX.lat > up.lat)
            up.lat = DX.lat;
         else if (DX.lat < lw.lat)
            lw.lat = DX.lat;
         if (DX.lon > up.lon)
            up.lon = DX.lon;
         else if (DX.lon < lw.lon)
            lw.lon = DX.lon;
      }
   }
   sprintf (interp->result, "%f %f %f %f", up.lat, up.lon, lw.lat, lw.lon);
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <convertShadow> :: Arthur Taylor TDL
 * Purpose:
 *   This returns the "Shadow" of a basin on a line defined by a point on
 *   the x-y plane, and a direction from north.
 *
 * Variables: (I=input) (O=output) (G=global)
 *    Dir (I) Is in radians.
 *
 * Returns:
 * History: 9/4/99 AAT Created.
 * Notes:
 ****************************************************************************/
#define PI 3.1415926535
#ifdef INCLUDE_TCL
static void convertShadow (convert_type * ct, double dir, double x0,
                           double y0, double *min_dist, double *max_dist)
{
   double x1, y1, dist;
   int i, j;

   for (i = 0; i < ct->bsn.ximax; i++) {
      for (j = 0; j < ct->bsn.yjmax; j += ct->bsn.yjmax - 1) {
         /* i+1, j+1 to convert to [1..imxb] system */
         pq2xy (i + 1, j + 1, &x1, &y1, &(ct->bsn));
         dist = (x1 - x0) * sin (dir) + (y1 - y0) * cos (dir); /* assumes
                                                                * "compass"
                                                                * dir */
         if ((i == 0) && (j == 0)) {
            *min_dist = dist;
            *max_dist = dist;
         } else {
            if (*min_dist > dist)
               *min_dist = dist;
            if (*max_dist < dist)
               *max_dist = dist;
         }
      }
   }
   for (i = 0; i < ct->bsn.ximax; i += ct->bsn.ximax - 1) {
      for (j = 1; j < ct->bsn.yjmax - 1; j++) {
         /* i+1, j+1 to convert to [1..imxb] system */
         pq2xy (i + 1, j + 1, &x1, &y1, &(ct->bsn));
         dist = (x1 - x0) * sin (dir) + (y1 - y0) * cos (dir); /* assumes
                                                                * "compass"
                                                                * dir */
         if (*min_dist > dist)
            *min_dist = dist;
         if (*max_dist < dist)
            *max_dist = dist;
      }
   }
}
#endif

/* Increased the number of points checked on the outside ring (from 8 to 128)
 * due to Ivan and "small" epns. */
/* Original error:
 *   If basin is smaller than R*Pi/4 where R is Radius of 34 knot winds,
 *   then it can "slip inside" the perimeter.  For Ivan, R was 243 miles,
 *   so the basin had to be smaller than 243 * pi / 4 = 190 miles.  The
 *   epns subgrid was 163 on the diagonal.
 * New algorithm:
 *   Switching to 128 sample points gives: R*PI/64, so the subgrid now has
 *   to be smaller than 12 miles.  Using R=300 as largest 34 knot wind radii
 *   gives limit of 15 miles...
 * Basins are more than 15 miles on an edge. */
/* Would rather use a "compute speed on edge of basin algortithm". */
int InsideBasin (convert_type * ct, LatLon pt, double d)
{
   int i;
   double temp, Px, Py;
   double Gx, Gy, p, q;

   cl2cnf (pt.lat, &pt.lat, &temp);
   sll2xy (pt, &Px, &Py, &(ct->slsh));
   for (i = 0; i < 128; i++) {
      Gx = Px + d * cos (PI / 64. * i);
      Gy = Py + d * sin (PI / 64. * i);
      xy2pq (Gx, Gy, &p, &q, &(ct->bsn));
      if ((p >= 1) && (p <= ct->bsn.ximax) && (q >= 1)
          && (q <= ct->bsn.yjmax)) {
         return 1;
      }
   }
   return 0;
}

/* Increased the number of points checked on the outside ring (from 8 to 128)
 * due to Ivan and "small" epns. */
/* Original error:
 *   If basin is smaller than R*Pi/4 where R is Radius of 34 knot winds,
 *   then it can "slip inside" the perimeter.  For Ivan, R was 243 miles,
 *   so the basin had to be smaller than 243 * pi / 4 = 190 miles.  The
 *   epns subgrid was 163 on the diagonal.
 * New algorithm:
 *   Switching to 128 sample points gives: R*PI/64, so the subgrid now has
 *   to be smaller than 12 miles.  Using R=300 as largest 34 knot wind radii
 *   gives limit of 15 miles...
 * Basins are more than 15 miles on an edge. */
/* Would rather use a "compute speed on edge of basin algortithm". */
int InsideBasin2 (convert_type * ct, LatLon pt, double d, int pMin, int pMax,
                  int qMin, int qMax)
{
   int i;
   double temp, Px, Py;
   double Gx, Gy, p, q;

   cl2cnf (pt.lat, &pt.lat, &temp);
   sll2xy (pt, &Px, &Py, &(ct->slsh));
   for (i = 0; i < 128; i++) {
      Gx = Px + d * cos (PI / 64. * i);
      Gy = Py + d * sin (PI / 64. * i);
      xy2pq (Gx, Gy, &p, &q, &(ct->bsn));
      if ((p >= pMin) && (p <= pMax) && (q >= qMin) && (q <= qMax)) {
         return 1;
      }
   }
   return 0;
}

/*****************************************************************************
 * <conInsideCmd> :: Arthur Taylor TDL  (halo_conInside)
 * Purpose:
 *   Assumes basin already loaded.
 *     Determines if a point is inside or outside a basin
 *    IDEA: Take point, convert to xy plane, Represent neighborhood around
 *      point using an 8 point sample. Find out if any of those are inside
 *      the basin.  If not... point is outside.  Otherwise inside.
 *    Analysis.. not completely accurate.  It is simple.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     lat/lon    (I) Point in question.
 *     dist       (I) Allowed distance (mi)
 * Returns: TCL_ERROR or TCL_OK
 *          Gives Tcl/Tk: 1 if inside 0 if outside
 * History: 11/16/1999 AAT Created.
 * Notes:
 *   Seems a shame to keep having to recalculate the center of the basin.
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int conInsideCmd (ClientData clientData, Tcl_Interp * interp,
                         int argc, char **argv)
{
   int f_in;
   convert_type *ct = (convert_type *) clientData;
   LatLon pt;
   double d;

   if ((argc != 4)) {
      sprintf (interp->result, "usage: %.50s <lat> <lon> <dist>", argv[0]);
      return TCL_ERROR;
   }
   pt.lat = atof (argv[1]);
   pt.lon = atof (argv[2]);
   d = atof (argv[3]);
   f_in = InsideBasin (ct, pt, d);
   sprintf (interp->result, "%d", f_in);
   return TCL_OK;
}
#endif

#define PI_180 0.01745329251994
#define C180_PI 57.29577951308
/* See mercator.c DistComputeCmd for more details */
/* dist is in stat mi, bear is in deg from north */
#ifdef INCLUDE_TCL
static void Earth_Jump (LatLon pt1, double dist, double bear, LatLon * pt2)
{
   double R, A, B, C;

   pt1.lat *= PI_180;
   pt1.lon *= PI_180;
   /* Need dist in nautical miles. */
   dist = dist / 1.151;
   R = 60. * C180_PI;
   dist = dist / R;
   bear = PI - (bear * PI_180);

   pt2->lat =
      asin (cos (dist) * sin (pt1.lat) -
            sin (dist) * cos (bear) * cos (pt1.lat));
   A = cos (dist) * cos (pt1.lat);
   B = sin (dist) * cos (bear) * sin (pt1.lat);
   C = sin (dist) * sin (bear);
   pt2->lon =
      atan2 (A * sin (pt1.lon) + B * sin (pt1.lon) + C * cos (pt1.lon),
             A * cos (pt1.lon) + B * cos (pt1.lon) - C * sin (pt1.lon));
   pt2->lat = pt2->lat * C180_PI;
   pt2->lon = pt2->lon * C180_PI;
}
#endif

/*****************************************************************************
 * <conMeowGenCmd> :: Arthur Taylor TDL  (halo_conMeowGen)
 * Purpose:
 *   Assumes basin already loaded.
 *     To generate a set of tracks for a MEOW based on a lat/lon point, and
 *   a direction.  If num != -1 it creates that many tracks otherwise it
 *   creates enough tracks to cover the basin in the x-y plane, and a little
 *   more (uses the R34 value for how much beyond).
 *     Returns a list of the resulting Seed points (as lat/lon), which lie on
 *   the line (in the x-y plane) normal to the line (in the x-y plane) defined
 *   by the given seed point and direction.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     XYFlag     (I) 0=Tracks on XY plane. 1=Tracks on Spherical earth
 *     lat/lon    (I) Central Seed point.
 *     dir        (I) Direction (0 is north.) (in degrees)
 *     apart      (I) Distance between the tracks. (mi)
 *     num        (I) -1 means use as many as you need.
 *     r34        (I) The radius of 34 knot winds for this storm.
 * Returns: TCL_ERROR or TCL_OK
 *          Gives Tcl/Tk: list of seed points.
 * History: 11/4/1999 AAT Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int conMeowGenCmd (ClientData clientData, Tcl_Interp * interp,
                          int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   LatLon pt, seed;
   double dir, Dir, r34, min_dist, max_dist, temp, x0, y0, apart;
   double start, cur, x1, y1;
   int num, cnt, cnt2, XYFlag;
   char f_cont;
   char tempBuffer[100];

   if ((argc != 8)) {
      sprintf (interp->result,
               "usage: %.50s <XYFlag 0-XY plane, 1-Spherical Earth> "
               " <lat> <lon> <dir> <apart> <num or -1> <r34>", argv[0]);
      return TCL_ERROR;
   }
   XYFlag = atoi (argv[1]);
   pt.lat = atof (argv[2]);
   pt.lon = atof (argv[3]);
   dir = atof (argv[4]) + 90;
   if (dir > 360) {
      dir -= 360;
   }
   apart = atof (argv[5]);
   num = atoi (argv[6]);
   r34 = atof (argv[7]);

   Dir = dir;
   dir = dir * PI / 180.;
   cl2cnf (pt.lat, &pt.lat, &temp);
   sll2xy (pt, &x0, &y0, &(ct->slsh));
   convertShadow (ct, dir, x0, y0, &min_dist, &max_dist);
   min_dist -= r34;
   max_dist += r34;

   start = 0;
   while (start < min_dist)
      start += apart;
   while (start > max_dist)
      start -= apart;

   /* guaranteed that min <= start <= max */
   cnt = 0;
   cur = start;
   x1 = x0 + cur * sin (dir);
   y1 = y0 + cur * cos (dir);
   sxy2ll (x1, y1, &seed, &(ct->slsh));
   cnf2cl (seed.lat, &seed.lat, &temp);
   sprintf (tempBuffer, "\"%d %f %f\" ", cnt, seed.lat, seed.lon);
   Tcl_AppendResult (interp, tempBuffer, NULL);
   f_cont = 1;
   while (f_cont == 1) {
      cnt++;
      cur += apart;
      if ((cur > max_dist) || ((num != -1) && (cnt >= (num / 2) + 1))) {
         f_cont = 0;
      } else {
         if (XYFlag == 0) {
            /* jump forward and store */
            x1 = x0 + cur * sin (dir);
            y1 = y0 + cur * cos (dir);
            sxy2ll (x1, y1, &seed, &(ct->slsh));
            cnf2cl (seed.lat, &seed.lat, &temp);
            sprintf (tempBuffer, "\"%d %f %f\" ", cnt, seed.lat, seed.lon);
            Tcl_AppendResult (interp, tempBuffer, NULL);
         } else {
            Earth_Jump (seed, apart, Dir, &seed);
            sprintf (tempBuffer, "\"%d %f %f\" ", cnt, seed.lat, seed.lon);
            Tcl_AppendResult (interp, tempBuffer, NULL);
         }
      }
   }
   /* reset seed point. */
   cur = start;
   if (XYFlag == 1) {
      x1 = x0 + cur * sin (dir);
      y1 = y0 + cur * cos (dir);
      sxy2ll (x1, y1, &seed, &(ct->slsh));
      cnf2cl (seed.lat, &seed.lat, &temp);
   }
   f_cont = 1;
   cnt2 = 0;
   while (f_cont == 1) {
      cnt++;
      cnt2++;
      cur -= apart;
      if ((cur < min_dist) || ((num != -1) && (cnt >= num))) {
         f_cont = 0;
      } else {
         if (XYFlag == 0) {
            /* jump forward and store */
            x1 = x0 + cur * sin (dir);
            y1 = y0 + cur * cos (dir);
            sxy2ll (x1, y1, &seed, &(ct->slsh));
            cnf2cl (seed.lat, &seed.lat, &temp);
            sprintf (tempBuffer, "\"-%d %f %f\" ", cnt2, seed.lat, seed.lon);
            Tcl_AppendResult (interp, tempBuffer, NULL);
         } else {
            Earth_Jump (seed, -apart, Dir, &seed);
            sprintf (tempBuffer, "\"-%d %f %f\" ", cnt2, seed.lat, seed.lon);
            Tcl_AppendResult (interp, tempBuffer, NULL);
         }
      }
   }
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <conMeowTraceCmd> :: Arthur Taylor TDL  (halo_conMeowTrace)
 * Purpose:
 *   Assumes basin already loaded.
 *     To follow a line (defined by seed point and dir) on the x-y plane
 *     Also to allow one to find the nearest point to a given lat/lon point
 *   on a given x-y line.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     lat/lon    (I) Central Seed point.
 *     dir        (I) Direction (0 is north.) (in degrees)
 *     flag       (I) 0 we want to follow line a given distance.
 *                    1 we want nearest lat/lon on line.
 *     dist or Lat (I) See flag
 *     ---- or Lon (I) See flag
 *
 * Returns: TCL_ERROR or TCL_OK
 *          Gives Tcl/Tk: desired lat/lon on line.
 * History: 11/4/1999 AAT Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int conMeowTraceCmd (ClientData clientData, Tcl_Interp * interp,
                            int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   LatLon pt, pt2, ans;
   double dir, temp, x0, y0, x1, y1, dist;
   char flag;

   if ((argc != 7) && (argc != 6)) {
      sprintf (interp->result,
               "usage: %.50s <lat> <lon> <dir> <flag 0-dist along line,"
               "1-nearest lat/lon to lat/lon> <dist or lat> <lon>", argv[0]);
      return TCL_ERROR;
   }
   pt.lat = atof (argv[1]);
   pt.lon = atof (argv[2]);
   cl2cnf (pt.lat, &pt.lat, &temp);
   sll2xy (pt, &x0, &y0, &(ct->slsh));
   dir = atof (argv[3]);
   if (dir > 360) {
      dir -= 360;
   }
   dir = dir * PI / 180.;
   flag = (char) atoi (argv[4]);
   if (flag == 0) {
      dist = atof (argv[5]);
   } else if (flag == 1) {
      pt2.lat = atof (argv[5]);
      pt2.lon = atof (argv[6]);
      cl2cnf (pt2.lat, &pt2.lat, &temp);
      sll2xy (pt2, &x1, &y1, &(ct->slsh));
      dist = (x1 - x0) * sin (dir) + (y1 - y0) * cos (dir); /* assumes
                                                             * "compass" dir 
                                                             */
   } else {
      sprintf (interp->result,
               "usage: %.50s <lat> <lon> <dir> <flag 0-dist along line,"
               "1-nearest lat/lon to lat/lon> <dist or lat> <lon>", argv[0]);
      return TCL_ERROR;
   }
   x1 = x0 + dist * sin (dir);
   y1 = y0 + dist * cos (dir);
   sxy2ll (x1, y1, &ans, &(ct->slsh));
   cnf2cl (ans.lat, &ans.lat, &temp);
   sprintf (interp->result, "%f %f", ans.lat, ans.lon);
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <con_ltln2pqCmd> :: Arthur Taylor TDL  (halo_ltln2pq)
 * Purpose:
 *     To convert from lat/lon to p/q.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     abrev     (I) Basin abreviation  (For validation only)
 *     type      (I) Basin type (e/p/h) (For validation only)
 *     lat       (I) Latitude on clarke spheroid.  NOT! mercator.
 *     lon       (I) Longitude on clarke spheroid. NOT! mercator.
 *
 * Returns: p q in a list
 * History: 2/2/98 AAT Created.
 * Notes:
 *   May want to change so it returns "Outside Grid" if it is.
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int con_ltln2pqCmd (ClientData clientData, Tcl_Interp * interp,
                           int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   LatLon DX;
   double p, q;

   if (argc != 5) {
      sprintf (interp->result, "usage: %.50s <abrev> <type> <lat> <lon> {On "
               "clarke spheroid, not mercator}", argv[0]);
      return TCL_ERROR;
   }
   if ((argv[2][0] != ct->bsn.type)
       || (strncmp (argv[1], ct->bsn.code, 3) != 0)) {
      sprintf (interp->result,
               "We Loaded %s type %c, which is not %s type %s\n",
               ct->bsn.code, ct->bsn.type, argv[1], argv[2]);
      return TCL_ERROR;
   }
   DX.lat = atof (argv[3]);
   DX.lon = atof (argv[4]);
   ltln2pq (DX, &p, &q, &(ct->bsn), &(ct->slsh));
   sprintf (interp->result, "%15.10f %15.10f", p, q);
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <con_pq2ltlnCmd> :: Arthur Taylor TDL  (halo_pq2ltln)
 * Purpose:
 *     To convert from p/q to lat/lon.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the basin data.
 *   interp      (I) The Tcl/Tk interpreter.
 *   argv:
 *     abrev     (I) Basin abreviation  (For validation only)
 *     type      (I) Basin type (e/p/h) (For validation only)
 *     p         (I) p element on SLOSH grid.
 *     q         (I) q element on SLOSH grid.
 *
 * Returns: p q in a list
 * History: 2/2/98 AAT Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static int con_pq2ltlnCmd (ClientData clientData, Tcl_Interp * interp,
                           int argc, char **argv)
{
   convert_type *ct = (convert_type *) clientData;
   LatLon DX;
   double p, q;

   if (argc != 5) {
      sprintf (interp->result, "usage: %.50s <abrev> <type> <p> <q>",
               argv[0]);
      return TCL_ERROR;
   }
   if ((argv[2][0] != ct->bsn.type)
       || (strncmp (argv[1], ct->bsn.code, 3) != 0)) {
      sprintf (interp->result,
               "We Loaded %s type %c, which is not %s type %s\n",
               ct->bsn.code, ct->bsn.type, argv[1], argv[2]);
      return TCL_ERROR;
   }
   p = atof (argv[3]);
   q = atof (argv[4]);
   pq2ltln (p, q, &DX, &(ct->bsn), &(ct->slsh));
   sprintf (interp->result, "%15.10f %15.10f", DX.lat, DX.lon);
   return TCL_OK;
}
#endif

/*****************************************************************************
 * <convertLoadCmdDel> :: Arthur Taylor TDL
 * Purpose:
 *     To free the memory associated with the conversion commands.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   clientData  (I) A pointer to the conversion data.
 *
 * Returns: NULL
 * History: 2/2/98 AAT Created.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
static void convertLoadCmdDel (ClientData clientData)
{
   convert_type *ct = (convert_type *) clientData;
   free (ct);
}
#endif

/*****************************************************************************
 * <Convert_Init> :: Arthur Taylor TDL
 * Purpose:
 *     To Initialize the commands needed to convert a lat/lon to SLOSH p/q
 *   co-ordinates and vice-versa.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   interp      (I) The Tcl/Tk interpreter.
 *
 * Returns: TCL_OK;
 * History: 2/2/98 AAT Commented.
 * Notes:
 ****************************************************************************/
#ifdef INCLUDE_TCL
int Convert_Init (Tcl_Interp * interp)
{
   convert_type *ct;

   ct = (convert_type *) malloc (sizeof (convert_type));
   Tcl_CreateCommand (interp, "halo_convert_Load", convertLoadCmd,
                      (ClientData) ct,
                      (Tcl_CmdDeleteProc *) convertLoadCmdDel);
   Tcl_CreateCommand (interp, "halo_conBasin_Dim", convertBasinDimCmd,
                      (ClientData) ct,
                      (Tcl_CmdDeleteProc *) convertLoadCmdDel);
   Tcl_CreateCommand (interp, "halo_conLoad_llxSet", convertLoadSetCmd,
                      (ClientData) ct,
                      (Tcl_CmdDeleteProc *) convertLoadCmdDel);
   Tcl_CreateCommand (interp, "halo_conSave_llxSet", convertSaveLLXCmd,
                      (ClientData) ct,
                      (Tcl_CmdDeleteProc *) convertLoadCmdDel);
   Tcl_CreateCommand (interp, "halo_ltln2pq", con_ltln2pqCmd,
                      (ClientData) ct, (Tcl_CmdDeleteProc *) NULL);
   Tcl_CreateCommand (interp, "halo_pq2ltln", con_pq2ltlnCmd,
                      (ClientData) ct, (Tcl_CmdDeleteProc *) NULL);
   Tcl_CreateCommand (interp, "halo_conBound", convertBoundCmd,
                      (ClientData) ct, (Tcl_CmdDeleteProc *) NULL);
   Tcl_CreateCommand (interp, "halo_conMeowGen", conMeowGenCmd,
                      (ClientData) ct, (Tcl_CmdDeleteProc *) NULL);
   Tcl_CreateCommand (interp, "halo_conMeowTraceGen", conMeowTraceCmd,
                      (ClientData) ct, (Tcl_CmdDeleteProc *) NULL);
   Tcl_CreateCommand (interp, "halo_conInside", conInsideCmd,
                      (ClientData) ct, (Tcl_CmdDeleteProc *) NULL);
   return TCL_OK;
}
#endif
