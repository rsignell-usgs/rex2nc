/*****************************************************************************
 * wind.c
 *
 * DESCRIPTION:
 *    This module computes and displays the SLOSH wind field.
 *
 * HISTORY
 *   9/1998 Arthur Taylor RDC/TDL: Created
 *  10/2003 AAT (RSIS/MDL): Revisited
 *
 * NOTES
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "wind.h"

/* Conversion from 1 minute winds to 10 minute winds or vice versa. */
#define w1_w10 1.15
/* 1 nautical = 1.15... miles */
/* X in miles, then Y = mi_nautical * X */
#define MI_NAUTICAL 1.15077944802

/* #define PI_360 0.008726646259972 */
#ifndef PI_180
#define PI_180 0.01745329251994
#endif
/* #define C180_PI 57.29577951308 */

#define PI 3.14159265359
#define OMEGA 7.292116E-5
/* #define OMEGA 7.292115854322E-5  (assuming 86164.09054 sec) */
/* Used PI=3.1415926535 */
#define RHOA 2.298e-3
/* RHOA air density (MB hr^2/mi^2) */
#define RHOWG 29.89
/* RHOWG water density*gravity (MB/ft) */
#define RHOW 4905.
/* RHOW water density (MB s^2/(ft*mi) */
#define EPSLN 1.e-4

#define MIN(A,B) ((A)<(B))?(A):(B)
#define MAX(A,B) ((A)>(B))?(A):(B)
#define PYTHAG(X) (sqrt((1.-(X))*(1.+(X))))
#define IC12 800

/* MAX_RMAX = 138 because of legacy code. */
#define MAX_RMAX 137.
#define MIN_RMAX 1.

/* MAX_PDROP = 200 because I can't imagine a pressure of < 812 mb. */
/* MAX_PDROP = 158 because force1 isn't designed for more (need more
   iterations). */
#define MAX_PDROP 200.
#define MIN_PDROP 2.

/* MAX_VMAX = 318 because that is the limit of a tornado in MPH. */
#define MAX_VMAX 318.
#define MIN_VMAX 1.

/*****************************************************************************
 * intwnd() --
 *
 * Chester Jelesnianski / TDL
 *
 * PURPOSE
 *    Sets up an initial guess for the max wind of a stationary storm using
 * initial input storm parameters.  The initial guess is an entry value to
 * force1() for iteration to exact maximum winds.
 *
 * ARGUMENTS
 * pdrop = Storm's Delta (MBS) (Input)
 *   lat = Latitude of storm (deg) (Input)
 *     m = Slope of (vmax = m * rmax + b) (Output)
 *     b = Intercept of (vmax = m * rmax + b) (Output)
 *
 * RETURNS: vmax (Maximum wind mph)
 *
 * HISTORY
 *   8/1980 Chester Jelesnianski: Created
 *   9/1998 Arthur Taylor (RDC/TDL): Converted to C.
 *  10/2003 AAT (RSIS/MDL): Revisited.
 *   9/2004 AAT (MDL): Switched to returning m and b.
 *
 * NOTES
 *   For a given pressure drop the max wind varies almost linearly with storm
 * size.  (See a wind/radius plot)  The slope and intercept of the pressure
 * line varies weakly with latitude or coriolis.  The given slope (wm(kpres))
 * and intercept (wb(kpres)) are for latitude 30 deg.  The entry pressure at
 * time 'ibgnt' is modified to correct for lat other than 30 deg.
 ****************************************************************************/
static void intwnd (double pdrop, double lat, float *m, float *b)
{
   float WM[150] =
      { .10, .14, .17, .20, .23, .25, .28, .30, .32, .33, .35, .37, .38, .39,
      .41, .42, .43, .44, .45, .46, .46, .47, .48, .49, .49, .50, .51, .52,
      .52, .53,
      .53, .54, .54, .55, .56, .56, .57, .57, .58, .58, .59, .59, .60, .60,
      .61, .62,
      .62, .63, .63, .64, .65, .65, .66, .66, .67, .68, .69, .69, .70, .71,
      .71, .72,
      .72, .73, .73, .74, .74, .75, .75, .76, .76, .76, .77, .77, .77, .77,
      .78, .78,
      .78, .78, .79, .79, .79, .79, .80, .80, .80, .80, .80, .81, .81, .81,
      .81, .81,
      .82, .82, .82, .82, .82, .82, .83, .83, .83, .83, .84, .84, .84, .84,
      .84, .84,
      .84, .85, .85, .85, .85, .85, .85, .85, .86, .86, .86, .86, .86, .86,
      .86, .86,
      .86, .86, .87, .87, .87, .87, .87, .87, .87, .87, .87, .87, .87, .87,
      .87, .87,
      .87, .87, .87, .87, .87, .87, .87, .87
   };
   float WB[150] =
      { 12., 23., 29., 33., 37., 39., 41., 43., 45., 47., 49., 51., 53., 55.,
      57., 58., 60., 62., 63., 65., 67., 68., 70., 71., 73., 74., 76., 77.,
      78., 80.,
      81., 82., 84., 85., 86., 88., 89., 90., 92., 93., 94., 95., 96., 98.,
      99., 100.,
      101., 102., 104., 105., 106., 107., 108., 109., 110., 112., 113., 114.,
      115.,
      116., 117., 118., 119., 120., 121., 122., 123., 123., 124., 125., 126.,
      127.,
      128., 129., 130., 131., 132., 132., 133., 134., 135., 136., 137., 137.,
      138.,
      139., 140., 141., 141., 142., 143., 144., 145., 145., 146., 147., 148.,
      148.,
      149., 150., 151., 151., 152., 153., 154., 154., 155., 156., 157., 157.,
      158.,
      158., 159., 160., 160., 161., 162., 162., 163., 164., 164., 165., 166.,
      166.,
      167., 167., 168., 169., 169., 170., 171., 171., 172., 172., 173., 174.,
      174.,
      175., 175., 176., 176., 177., 177., 178., 179., 179., 180., 180., 181.,
      182.
   };
   double pf;
   int kpres;

   pf = pdrop / 67.7;
   kpres = pdrop + (pf * pf * .3 * (30 - lat)) + 0.5;
   if (kpres > 151) {
      /* 9/2004 fixed following from 140 to 150. */
      kpres = 150;
   } else if (kpres < 1) {
      kpres = 1;
   }
   /* -1 on next line because we are in C, not fortran */
   *m = -WM[kpres - 1];
   *b = WB[kpres - 1];
}

/*****************************************************************************
 * force1() --
 *
 * Jye Chen / TDL
 *
 * PURPOSE
 *    This procedure computes at mile intervals from the center, for 800 miles,
 * (1) sin of inflow angle, (2) cos of inflow angle, (3) surface pressure,
 * (4) static height.
 *    Refrence for the algorithm and theory is:
 * "Jelesnianski and Taylor, 1973, A preliminary view of storm surges before
 *  and after storm modifications, NOAA Technical Memorandum ERL WMPO-3,
 *  Weather Modification Program Office, ERL NOAA, U. S. Department of
 *  Commerce, Washington D.C., May, p. 33"
 *
 * ARGUMENTS
 *   pdrop = Storm's Delta (MBs) (Input)
 *    rmax = Storm's Radius of max winds (st. MI) (Input)
 *     lat = Latitude of storm (clarke deg) (Input)
 *    S, C = sin and cos of inflow angles (Output)
 * P, delp = static heights, and surface pressure. (MB/st mile) (ft) (Output)
 *    Vmax = Maximum wind (mph) (Input/Output)
 *
 * RETURNS: (int) 0 if ok, -1 if error occured
 *
 * HISTORY
 *   8/1985 Jye Chen (TDL): Created
 *   9/1998 Arthur Taylor (RDC/TDL): Converted to C.
 *  10/2003 AAT (RSIS/MDL): Revisited
 *
 * NOTES
 ****************************************************************************/
static int force1 (double pdrop, double rmax, double lat, double *S,
                   double *C, double *P, double *delp, double *Vmax)
{
   int jset;
   double corhr, cdrag = 1.;
   double vmax, rmax_2, VRmax;
   double xks, xkn;     /* Wind friction values */
   double r0, r1, delr; /* Radius of max winds */
   double veer0, veer1; /* Tangental velocity * r0, or r1 */
   double cotphi, sinphi, cosphi, osinfi, dosinf;
   double w0, w1, fact, slp0, slp1;
   double padd;         /* Pressure (MB) from a radius (jset) to infinity */
   double dgdw, delw, dpdr1, dpdr0;
   double press, prs, pinft;
   double term_r0;
   char i_stop, kk_stop, main_stop;
   int i, kk, k, j2;
   float intwnd_m, intwnd_b;

   if (lat < 0)
      lat = -1 * lat;

   jset = MIN (IC12 - 51, MAX (150, (int) (12 * rmax)));
   corhr = 2 * (OMEGA) * sin (lat * PI_180) * 3600;
   intwnd (pdrop, lat, &intwnd_m, &intwnd_b);
   vmax = intwnd_m * rmax + intwnd_b;
   rmax_2 = rmax * rmax;
   VRmax = vmax * rmax;
   /* if lake winds vmax = .85 * vmax, and cdrag = 4 * sqrt(22 / rmax) */
   /* Begin iteration */
   main_stop = 0;
   do {
      /* form 2 frictional values for a storm */
      xks = .01 * sqrt (rmax / (.3 * vmax + 60.)) * cdrag;
      xkn = .01 * sqrt (rmax / (.4 * vmax + 80.)) * cdrag;

      j2 = MIN (3. * rmax, 50.);
      r0 = jset + 50;
      delr = 1.;
      term_r0 = r0 / (rmax_2 + r0 * r0);
      veer0 = (2 * VRmax * term_r0) * r0;
      cotphi = (corhr * r0 / veer0 + xkn) / xks;
      sinphi = 1. / sqrt (1. + cotphi * cotphi);
      cosphi = PYTHAG (sinphi);
      osinfi = 1. / sinphi;
      w0 = veer0 * cosphi;
      press = 0.;
      dpdr0 = RHOA * (veer0 / r0) * (veer0 / r0) *
         (xks / sinphi - (1. / r0 - 2. * term_r0));

      /* Runge-Kutta-Heun-Taylor method to solve: dwdr = w * xks / sin(phi)
       * - corhr * r - xkn * vr = G(w,r) where w = v * r * cos(phi) */
      i_stop = 0;
      fact = 0;
      for (i = 1; ((i <= 1500) && (i_stop == 0)); i++) {
         r1 = r0 - delr;
         if (r1 < delr - .01) {
            /* break out of for loop...line 165 */
            i_stop = 1;
         } else {
            fact = 2. / delr;
            slp0 = w0 * xks / sinphi - corhr * r0 - xkn * veer0;
            veer1 = 2 * VRmax * r1 / (rmax_2 + r1 * r1) * r1;
            w1 = veer1 * cosphi;
            /* Newton-Raphson method to find w1 so that .5 * (G(w1,r1) +
             * G(w0,r0)) = (w0 - w1) / (r0 - r1) */
            kk_stop = 0;
            for (kk = 1; ((kk <= 100) && (kk_stop == 0)); kk++) {
               slp1 = w1 * xks / sinphi - corhr * r1 - xkn * veer1;
               padd = fact * (w0 - w1) - slp0;
               dgdw = xks / (sinphi * sinphi * sinphi);
               delw = (slp1 - padd) / (dgdw + fact);
               dosinf = -delw * cosphi / (veer1 * sinphi * sinphi * sinphi);
               osinfi = osinfi + dosinf;
               sinphi = 1. / osinfi;
               cosphi = PYTHAG (sinphi);
               w1 = veer1 * cosphi;
               if (w1 <= 0) {
                  /* break out of iterations...line 180 */
#ifdef DEBUG
                  fprintf (stderr, "force 1 : ocos(phi) < 0\n");
#endif
                  /* ERROR ocos(phi) < 0 */
                  return -1;
               } else {
                  if (fabs (dosinf / osinfi) <= EPSLN) {
                     /* break out of kk loop. */
                     kk_stop = 1;
                  }
               }
            }           /* end of kk loop line 120 */
            if (kk_stop == 0) {
#ifdef DEBUG
               fprintf (stderr, "force 1 : insufficient steps in "
                        "Newton-Raphson loop\n");
#endif
               /* ERROR insufficient steps in Newton-Raphson loop */
               return -1;
            }
            r0 = r1;
            /* Don't need to update term_r0 inside i loop */
            k = r1 + 1.01;
            veer0 = veer1;
            w0 = w1;
            S[k - 1] = sinphi;
            C[k - 1] = cosphi;
            dpdr1 = RHOA * (veer1 / r1) * (veer1 / r1) * (xks / sinphi -
                                                          (1. / r1 -
                                                           2. * r1 /
                                                           (rmax_2 +
                                                            r1 * r1)));
            P[k - 1] = dpdr1;
            press = press - (dpdr0 + dpdr1) / fact;
            delp[k - 1] = press;
            dpdr0 = dpdr1;
            if ((k - j2 - 1) <= 0) {
               delr = .1;
            }
         }
      }                 /* end of i loop line 160 */
      delp[0] = press - dpdr0 / fact;
      P[0] = 0;
      S[0] = 0;
      C[0] = 1;
      r0 = jset - 1;
      term_r0 = r0 / (rmax_2 + r0 * r0);
      padd = 2. * RHOA * VRmax * VRmax *
         (xks * (atan (rmax / r0) / rmax + term_r0) / S[jset - 1]) +
         .5 * RHOA * ((2. * VRmax * term_r0) * (2. * VRmax * term_r0));
      pinft = padd + delp[jset - 1] - delp[0];
      prs = pdrop;
      if (fabs (pinft - prs) > .1) {
         vmax = vmax * sqrt (prs / pinft);
         VRmax = vmax * rmax;
      } else {
         main_stop = 1;
      }
   } while (main_stop == 0);
   /* line 200 */
   sinphi = S[jset - 1 - 1];
   cosphi = C[jset - 1 - 1];
   for (k = jset + 1; k <= IC12; k++) {
      r0 = k + 1;
      S[k - 1] = sinphi;
      C[k - 1] = cosphi;
      term_r0 = r0 / (rmax_2 + r0 * r0);
      P[k - 1] = RHOA * (2. * VRmax * term_r0) * (2. * VRmax * term_r0) *
         (xks / sinphi - (1. / r0 - 2. * term_r0));
      delp[k - 1] = delp[k - 1 - 1] + .5 * (P[k - 1] + P[k - 1 - 1]);
   }
   press = padd + delp[jset - 1];
   for (i = 1; i <= IC12; i++) {
      delp[i - 1] = (press - delp[i - 1]) / RHOWG;
      P[i - 1] = P[i - 1] / RHOW;
   }
   *Vmax = vmax;
   return 0;
}

/*****************************************************************************
 * Quickforce1() --
 *
 * Jye Chen / TDL
 *
 * PURPOSE
 *    This procedure computes at mile intervals from the center, for 800 miles,
 * (1) sin of inflow angle, (2) cos of inflow angle, (3) surface pressure,
 * (4) static height.
 *    Refrence for the algorithm and theory is:
 * "Jelesnianski and Taylor, 1973, A preliminary view of storm surges before
 *  and after storm modifications, NOAA Technical Memorandum ERL WMPO-3,
 *  Weather Modification Program Office, ERL NOAA, U. S. Department of
 *  Commerce, Washington D.C., May, p. 33"
 *
 * ARGUMENTS
 *   pdrop = Storm's Delta (MBs) (Input)
 *    rmax = Storm's Radius of max winds (st. MI) (Input)
 *     lat = Latitude of storm (clarke deg) (Input)
 *    S, C = sin and cos of inflow angles (Output)
 * P, delp = static heights, and surface pressure. (MB/st mile) (ft) (Output)
 *    Vmax = Maximum wind (mph) (Input/Output)
 *
 * RETURNS: (int) 0 if ok, -1 if error occured
 *
 * HISTORY
 *   8/1985 Jye Chen (TDL): Created
 *   9/1998 Arthur Taylor (RDC/TDL): Converted to C.
 *  10/2003 AAT (RSIS/MDL): Revisited
 *   9/2004 AAT (MDL): Made it so it doesn't fill out S,C,P,delp except
 *          at rmax.
 *
 * NOTES
 ****************************************************************************/
static int QuickForce1 (double pdrop, double rmax, double lat, double *S,
                        double *C, double *Vmax)
{
   int jset;
   double corhr, cdrag = 1.;
   double vmax, rmax_2, VRmax;
   double xks, xkn;     /* Wind friction values */
   double r0, r1, delr; /* Radius of max winds */
   double veer0, veer1; /* Tangental velocity * r0, or r1 */
   double cotphi, sinphi, cosphi, osinfi, dosinf;
   double w0, w1, fact, slp0, slp1;
   double padd;         /* Pressure (MB) from a radius (jset) to infinity */
   double dgdw, delw, dpdr1, dpdr0;
   double press, pinft;
   double term_r0;
   char i_stop, kk_stop, main_stop;
   int i, kk, k, j2;
   float intwnd_m, intwnd_b;
   double delp_jset1 = pdrop; /* set to something to avoid "warning." */

   if (lat < 0) {
      lat = -1 * lat;
   }

   jset = MIN (IC12 - 51, MAX (150, (int) (12 * rmax)));
   corhr = 2 * (OMEGA) * sin (lat * PI_180) * 3600;
   intwnd (pdrop, lat, &intwnd_m, &intwnd_b);
   vmax = intwnd_m * rmax + intwnd_b;
   rmax_2 = rmax * rmax;
   VRmax = vmax * rmax;
   /* if lake winds vmax = .85 * vmax, and cdrag = 4 * sqrt(22 / rmax) */
   /* Begin iteration */
   main_stop = 0;
   do {
      /* form 2 frictional values for a storm */
      xks = .01 * sqrt (rmax / (.3 * vmax + 60.)) * cdrag;
      xkn = .01 * sqrt (rmax / (.4 * vmax + 80.)) * cdrag;

      j2 = MIN (3. * rmax, 50.);
      r0 = jset + 50;
      delr = 1.;
      term_r0 = r0 / (rmax_2 + r0 * r0);
      veer0 = (2 * VRmax * term_r0) * r0;
      cotphi = (corhr * r0 / veer0 + xkn) / xks;
      sinphi = 1. / sqrt (1. + cotphi * cotphi);
      cosphi = PYTHAG (sinphi);
      osinfi = 1. / sinphi;
      w0 = veer0 * cosphi;
      press = 0.;
      dpdr0 = RHOA * (veer0 / r0) * (veer0 / r0) *
         (xks / sinphi - (1. / r0 - 2. * term_r0));

      /* Runge-Kutta-Heun-Taylor method to solve: dwdr = w * xks / sin(phi)
       * - corhr * r - xkn * vr = G(w,r) where w = v * r * cos(phi) */
      i_stop = 0;
      fact = 0;
      for (i = 1; ((i <= 1500) && (i_stop == 0)); i++) {
         r1 = r0 - delr;
         if (r1 < delr - .01) {
            /* break out of for loop...line 165 */
            i_stop = 1;
         } else {
            fact = 2. / delr;
            slp0 = w0 * xks / sinphi - corhr * r0 - xkn * veer0;
            veer1 = 2 * VRmax * r1 / (rmax_2 + r1 * r1) * r1;
            w1 = veer1 * cosphi;
            /* Newton-Raphson method to find w1 so that .5 * (G(w1,r1) +
             * G(w0,r0)) = (w0 - w1) / (r0 - r1) */
            kk_stop = 0;
            for (kk = 1; ((kk <= 100) && (kk_stop == 0)); kk++) {
               slp1 = w1 * xks / sinphi - corhr * r1 - xkn * veer1;
               padd = fact * (w0 - w1) - slp0;
               dgdw = xks / (sinphi * sinphi * sinphi);
               delw = (slp1 - padd) / (dgdw + fact);
               dosinf = -delw * cosphi / (veer1 * sinphi * sinphi * sinphi);
               osinfi = osinfi + dosinf;
               sinphi = 1. / osinfi;
               cosphi = PYTHAG (sinphi);
               w1 = veer1 * cosphi;
               if (w1 <= 0) {
                  /* break out of iterations...line 180 */
#ifdef DEBUG
                  fprintf (stderr, "force 1 : ocos(phi) < 0\n");
#endif
                  /* ERROR ocos(phi) < 0 */
                  return -1;
               } else {
                  if (fabs (dosinf / osinfi) <= EPSLN) {
                     /* break out of kk loop. */
                     kk_stop = 1;
                  }
               }
            }           /* end of kk loop line 120 */
            if (kk_stop == 0) {
#ifdef DEBUG
               fprintf (stderr, "force 1 : insufficient steps in "
                        "Newton-Raphson loop\n");
#endif
               /* ERROR insufficient steps in Newton-Raphson loop */
               return -1;
            }
            r0 = r1;
            /* Don't need to update term_r0 inside i loop */
            k = r1 + 1.01;
            veer0 = veer1;
            w0 = w1;
            S[k - 1] = sinphi;
            C[k - 1] = cosphi;
            dpdr1 = RHOA * (veer1 / r1) * (veer1 / r1) * (xks / sinphi -
                                                          (1. / r1 -
                                                           2. * r1 /
                                                           (rmax_2 +
                                                            r1 * r1)));
            press = press - (dpdr0 + dpdr1) / fact;
            if (k - 1 == jset - 1) {
               delp_jset1 = press;
            }
            dpdr0 = dpdr1;
            if ((k - j2 - 1) <= 0) {
               delr = .1;
            }
         }
      }                 /* end of i loop line 160 */
      S[0] = 0;
      C[0] = 1;
      r0 = jset - 1;
      term_r0 = r0 / (rmax_2 + r0 * r0);
      padd = 2. * RHOA * VRmax * VRmax *
         (xks * (atan (rmax / r0) / rmax + term_r0) / S[jset - 1]) +
         .5 * RHOA * ((2. * VRmax * term_r0) * (2. * VRmax * term_r0));
      pinft = padd + delp_jset1 - (press - dpdr0 / fact);
      if (fabs (pinft - pdrop) > .1) {
         vmax = vmax * sqrt (pdrop / pinft);
         VRmax = vmax * rmax;
      } else {
         main_stop = 1;
      }
   } while (main_stop == 0);
   /* line 200 */
   *Vmax = vmax;
   return 0;
}

/*****************************************************************************
 * QuickWind_Vmax()
 *
 * PURPOSE
 *    To Compute the SLOSH Vmax, given forward motion... Returns in MPH.
 *
 * ARGUMENTS
 *   pdrop = Storm's Delta Pressure (MBs) (Input)
 *    rmax = Storm's Radius of max winds (st. MI) (Input)
 *     lat = Latitude of storm (clarke deg) (Input)
 *     dir = Storm's direction 0 for north. (Input)
 *  wspeed = Storm's speed (MPH) (Input)
 *    Vmax = The Vmax given forward motion. (Output)
 *
 * Returns: (int) 0 if ok, -1 if error occured in force1
 *
 * History:
 *   9/2004 Arthur Taylor (MDL) Revised Wind_Vmax.  Doesn't care about S, C,
 *          P, delp
 *
 * NOTES
 ****************************************************************************/
static int QuickWind_Vmax (double pdrop, double rmax, double lat, double dir,
                           double wspeed, double *Vmax)
{
   double S[IC12], C[IC12]; /* Sin and Cos of the inflow angles. */
   double vmax, xp, yp, dr, ck1, sk1;
   double a, b, wspeed_cos, wspeed_sin;
   double phi, theta;
   int k;

   if (QuickForce1 (pdrop, rmax, lat, S, C, &vmax) != 0) {
#ifdef DEBUG
      fprintf (stderr, "Error in force 1 : %f %f %f %f \n", pdrop, rmax,
               lat, vmax);
#endif
      return -1;
   }

   dr = rmax - (int) (rmax);
   k = MIN ((int) (rmax), IC12 - 2);
   ck1 = C[k] + dr * (C[k + 1] - C[k]); /* inflow angle... */
   sk1 = S[k] + dr * (S[k + 1] - S[k]);

   phi = tan (sk1 / ck1);
   dir = 90. - dir;
   theta = (3 * PI) / 2. + (dir * PI_180 - phi);
   if (theta > 2 * PI) {
      theta -= 2 * PI;
   }
   xp = rmax * cos (theta);
   yp = rmax * sin (theta);

   /* Relative north is 90 deg mathematical */
   wspeed_cos = wspeed * cos (dir * PI_180);
   wspeed_sin = wspeed * sin (dir * PI_180);
   a = rmax * wspeed_cos - 2 * vmax * (yp * ck1 + xp * sk1);
   b = rmax * wspeed_sin + 2 * vmax * (xp * ck1 - yp * sk1);

   /* *Vmax = (rmax/(rmax*rmax + rs*rs))*sqrt((a*a + b*b)); rs == rmax */
   *Vmax = sqrt (a * a + b * b) / (2. * rmax); /* vmax is in MPH */
   return 0;
}

/*****************************************************************************
 * Wind_Vmax()
 *
 * PURPOSE
 *    To Compute the SLOSH Vmax, given forward motion... Returns in MPH.
 *
 * ARGUMENTS
 *   pdrop = Storm's Delta Pressure (MBs) (Input)
 *    rmax = Storm's Radius of max winds (st. MI) (Input)
 *     lat = Latitude of storm (clarke deg) (Input)
 *     dir = Storm's direction 0 for north. (Input)
 *  wspeed = Storm's speed (MPH) (Input)
 *    S, C = Sin and Cos of inflow angles. (Output)
 * P, delp = Pressure fields. (Output)
 *    Vmax = The Vmax given forward motion. (Output)
 *
 * Returns: (int) 0 if ok, -1 if error occured in force1
 *
 * History:
 *  10/1999 Arthur Taylor (RSIS/TDL) Created.
 *  10/2003 AAT (MDL/RSIS) Revisited. (Renamed from SLOSH_Vmax)
 *
 * NOTES
 ****************************************************************************/
int Wind_Vmax (double pdrop, double rmax, double lat, double dir,
               double wspeed, double *S, double *C, double *P,
               double *delp, double *Vmax)
{
   double vmax, xp, yp, dr, ck1, sk1;
   double a, b, wspeed_cos, wspeed_sin;
   double phi, theta;
   int k;

   if (force1 (pdrop, rmax, lat, S, C, P, delp, &vmax) != 0) {
      return -1;
   }

   dr = rmax - (int) (rmax);
   k = MIN ((int) (rmax), IC12 - 2);
   ck1 = C[k] + dr * (C[k + 1] - C[k]); /* inflow angle... */
   sk1 = S[k] + dr * (S[k + 1] - S[k]);

   phi = tan (sk1 / ck1);
   dir = 90. - dir;
   theta = (3 * PI) / 2. + (dir * PI_180 - phi);
   if (theta > 2 * PI) {
      theta -= 2 * PI;
   }
   xp = rmax * cos (theta);
   yp = rmax * sin (theta);

   /* Relative north is 90 deg mathematical */
   wspeed_cos = wspeed * cos (dir * PI_180);
   wspeed_sin = wspeed * sin (dir * PI_180);
   a = rmax * wspeed_cos - 2 * vmax * (yp * ck1 + xp * sk1);
   b = rmax * wspeed_sin + 2 * vmax * (xp * ck1 - yp * sk1);

   /* *Vmax = (rmax/(rmax*rmax + rs*rs))*sqrt((a*a + b*b)); rs == rmax */
   *Vmax = sqrt (a * a + b * b) / (2. * rmax); /* vmax is in MPH */
   return 0;
}

/*****************************************************************************
 * SLOSH_Vmax()
 *
 * PURPOSE
 *    To compute the SLOSH Vmax, given forward motion... This is the maximum
 * 1 minute average winds in MPH.
 *
 * ARGUMENTS
 * tenMinFact = Number to convert from 10 min avg winds to 1 min avg winds
 *      pdrop = Storm's Delta Pressure (MBs) (Input)
 *       rmax = Storm's Radius of max winds (st. MI) (Input)
 *        lat = Latitude of storm (clarke deg) (Input)
 *        dir = Storm's direction 0 for north. (Input)
 *      speed = Storm's forward speed (MPH) (Input)
 *       Vmax = The max 1 minute avg sustained windspeed. (Output)
 *
 * Returns: 0 if ok, -1 if error occured in force1
 *
 * History:
 *  10/2003 Arthur Taylor (MDL/RSIS) Created.
 *
 * NOTES
 ****************************************************************************/
int SLOSH_Vmax (double tenMinFact, double pdrop, double rmax, double lat,
                double dir, double speed, double *Vmax)
{
   if (QuickWind_Vmax (pdrop, rmax, lat, dir, speed, Vmax) != 0) {
      return -1;
   }
   *Vmax = *Vmax * tenMinFact;
   return 0;
}

/*****************************************************************************
 * SLOSH_Rmax()
 *
 * PURPOSE
 *    Given Vmax and Pressure, compute the SLOSH Rmax (in st. miles). It does
 * this through an iterative process of guessing an Rmax, computing the
 * associated Vmax, and checking if it is "close enough" to the given Vmax.
 *
 * ARGUMENTS
 * tenMinFact = Number to convert from 10 min avg winds to 1 min avg winds
 *      pdrop = Storm's Delta Pressure (MBs) (Input)
 *       Vmax = The max 1 minute avg sustained windspeed. (Input)
 *        lat = Latitude of storm (clarke deg) (Input)
 *        dir = Storm's direction 0 for north. (Input)
 *      speed = Storm's forward speed (MPH) (Input)
 *       Rmax = Storm's Radius of max winds (st. MI) (Input/Output)
 *
 * Returns: 0 if ok, -1 if error occured in force1
 *
 * History:
 *  10/2003 Arthur Taylor (MDL/RSIS) Created.
 *
 * NOTES
 *   After 17 steps the range of 200 reduces to .003, so the iteration should
 * stop after 16 steps.
 ****************************************************************************/
int SLOSH_Rmax (double tenMinFact, double pdrop, double Vmax, double lat,
                double dir, double speed, double *Rmax)
{
   double rmax;         /* Current guess for rmax. */
   double top = MAX_RMAX; /* Max value for rmax. */
   double bot = MIN_RMAX; /* Min value for rmax. */
   double vmax;         /* vmax assuming "rmax". */
   int i;               /* Make's sure the search isn't infinite. */
   float m, b;          /* describes linear relation between vmax and rmax. */

   /* Check input parameters... */
   if (pdrop < MIN_PDROP) {
#ifdef DEBUG
      fprintf (stderr, "delta p < %f is unrealistic\n", MIN_PDROP);
#endif
      pdrop = MIN_PDROP;
   }
   if (pdrop > MAX_PDROP) {
#ifdef DEBUG
      fprintf (stderr, "delta p > %f is unrealistic\n", MAX_PDROP);
#endif
      pdrop = MAX_PDROP;
   }
   if (Vmax > MAX_VMAX) {
#ifdef DEBUG
      fprintf (stderr, "vmax > %f is unrealistic\n", MAX_VMAX);
#endif
      Vmax = MAX_VMAX;
   }

   Vmax = Vmax / tenMinFact;
   if (Vmax < MIN_VMAX + speed / 2.) {
#ifdef DEBUG
      fprintf (stderr,
               "rmax: speed = %f, vmax (%f) < %f + %f = %f is unrealistic\n",
               speed, Vmax, MIN_VMAX, speed / 2., MIN_VMAX + speed / 2.);
#endif
      Vmax = MIN_VMAX + speed / 2.;
   }
   if (lat < 0) {
      intwnd (pdrop, -1 * lat, &m, &b);
   } else {
      intwnd (pdrop, lat, &m, &b);
   }
   rmax = (Vmax - (speed / 2.) - b) / m;
   if ((rmax < bot) || (rmax > top)) {
      rmax = (top + bot) / 2.0;
   }

   /* 138 / 2^16 < .01 */
   for (i = 0; i < 16; i++) {
      if (QuickWind_Vmax (pdrop, rmax, lat, dir, speed, &vmax) != 0) {
         *Rmax = rmax;
         return -1;
      }
      if ((Vmax - vmax) > .01) {
         top = rmax;
         rmax = (rmax + bot) / 2;
      } else if ((Vmax - vmax) < -.01) {
         bot = rmax;
         rmax = (rmax + top) / 2;
      } else {
         break;
      }
      if ((top - bot) < .01) {
         break;
      }
   }
   *Rmax = rmax;
   if ((MAX_RMAX - rmax) < .01) {
      /* fprintf (stderr, "caution: rmax max reached.\n"); */
      return 1;
   }
   if ((rmax - MIN_RMAX) < .01) {
      /* fprintf (stderr, "caution: rmax min reached\n"); */
      return 2;
   }
   return 0;
}

/*****************************************************************************
 * SLOSH_pDrop()
 *
 * PURPOSE
 *    Given Vmax and Rmax, compute the SLOSH Delta Pressure. It does this
 * through an iterative process of guessing a delta pressure, computing the
 * associated Vmax, and checking if it is "close enough" to the given Vmax.
 *
 * ARGUMENTS
 * tenMinFact = Number to convert from 10 min avg winds to 1 min avg winds
 *       rmax = Storm's Radius of max winds (st. MI) (Output)
 *       Vmax = The max 1 minute avg sustained windspeed. (Input)
 *        lat = Latitude of storm (clarke deg) (Input)
 *        dir = Storm's direction 0 for north. (Input)
 *      speed = Storm's forward speed (MPH) (Input)
 *      PDrop = Storm's Delta Pressure (MBs) (Input/Output)
 *
 * Returns: 0 if ok, -1 if error occured in force1
 *
 * History:
 *  10/2003 Arthur Taylor (MDL/RSIS) Created.
 *
 * NOTES
 *   After 17 steps the range of 200 reduces to .003, so the iteration should
 * stop after 16 steps.
 ****************************************************************************/
int SLOSH_pDrop (double tenMinFact, double rmax, double Vmax, double lat,
                 double dir, double speed, double *PDrop)
{
   double pDrop;        /* Current guess for delta pressure. */
   double top = MAX_PDROP; /* Max value for delta pressure. */
   double bot = MIN_PDROP; /* Min value for delta pressure. */
   double vmax;         /* vmax assuming "pDrop". */
   size_t i;            /* Make's sure the search isn't infinite. */

   /* Check input parameters... */
   if (rmax < MIN_RMAX) {
#ifdef DEBUG
      fprintf (stderr, "rmax < %f is unrealistic\n", MIN_RMAX);
#endif
      rmax = MIN_RMAX;
   }
   if (rmax > MAX_RMAX) {
#ifdef DEBUG
      fprintf (stderr, "rmax > %f is unrealistic\n", MAX_RMAX);
#endif
      rmax = MAX_RMAX;
   }
   if (Vmax > MAX_VMAX) {
#ifdef DEBUG
      fprintf (stderr, "vmax > %f is unrealistic\n", MAX_VMAX);
#endif
      Vmax = MAX_VMAX;
   }

   pDrop = *PDrop;
   if ((pDrop < bot) || (pDrop > top)) {
      pDrop = (top + bot) / 2.0;
   }
   Vmax = Vmax / tenMinFact;
   if (Vmax < MIN_VMAX + speed / 2.) {
#ifdef DEBUG
      fprintf (stderr,
               "pdrop: speed = %f, vmax (%f) < %f + %f = %f is unrealistic\n",
               speed, Vmax, MIN_VMAX, speed / 2., MIN_VMAX + speed / 2.);
#endif
      Vmax = MIN_VMAX + speed / 2.;
   }
   for (i = 0; i < 16; i++) {
      if (QuickWind_Vmax (pDrop, rmax, lat, dir, speed, &vmax)
          != 0) {
         *PDrop = pDrop;
         return -1;
      }
      if ((Vmax - vmax) > .01) {
         bot = pDrop;
         pDrop = (pDrop + top) / 2;
      } else if ((Vmax - vmax) < -.01) {
         top = pDrop;
         pDrop = (pDrop + bot) / 2;
      } else {
         break;
      }
      if ((top - bot) < .01) {
         break;
      }
   }
   *PDrop = pDrop;
   return 0;
}

/* pdrop in MBs, lat/lon in negative west, rmax in stat mile,
   wspeed is storm forward motion in MPH, dir is storm direction from north,
   U, V are in MPH (magnitude is srqt (sum square) in 1 min avg MPH)  
 */
int windProbe (LatLon track, double pdrop, double rmax, double wspeed,
               double dir, LatLon station, double *U, double *V) {
   static double S[800],C[800],P[800],delp[800];
   static double vmax = 0;
   static double lclPdrop = 0;
   static double lclRmax = 0;
   static double lclLat = 0;
   double xp, yp, rmax_2, rsq, rs, dr, ck1, sk1;
   double a, b, cc, wspeed_cos, wspeed_sin;
   int k;

/*
   if (argc != 9) {
    sprintf (interp->result, "usage: %.50s <storm lon> <storm lat clarke deg> "
          "<storm pdrop MBs> <storm rmax st. miles> <storm wspeed MPH> "
          "<storm direction deg> <station lon> <station lat clarke deg>", argv[0]);
    return TCL_ERROR;
  }
  track.lon = -atof(argv[1]);
  track.lat = atof(argv[2]);
  pdrop = atof(argv[3]);
  rmax = atof(argv[4]);
  wspeed = atof(argv[5]);
  dir = 90.-atof(argv[6]);
  station.lon = -atof(argv[7]);
  station.lat = atof(argv[8]);
*/
  while (track.lon > 360) track.lon -= 360;
  while (track.lon < 0) track.lon += 360;
  while (station.lon > 360) station.lon -= 360;
  while (station.lon < 0) station.lon += 360;

  if ((pdrop != lclPdrop) || (rmax != lclRmax) || (track.lat != lclLat)) {
     lclPdrop = pdrop;
     lclRmax = rmax;
     lclLat = track.lat;
     if (force1 (lclPdrop, lclRmax, lclLat, S, C, P, delp, &vmax) != 0) {
        printf ("force1 returned -1 to (windProbe)");
        return -1;
     }
  }

  xp = 60 * MI_NAUTICAL * (station.lon - track.lon) * cos((station.lat + track.lat)/2. * PI/180.);
  yp = 60 * MI_NAUTICAL * (station.lat - track.lat);
  rmax_2 = rmax*rmax;

  /* Relative north is 90 deg mathematical */
  wspeed_cos = wspeed * cos(dir*PI_180);
  wspeed_sin = wspeed * sin(dir*PI_180);
/*
  Needed only for lake winds.
  X12_23 = wspeed * cos(dir*PI_180+7.*PI/6.-.5*PI);
  X12_24 = wspeed * sin(dir*PI_180+7.*PI/6.-.5*PI);
*/
  rsq = xp*xp+yp*yp;
  rs = sqrt(rsq);
  k = (int)(rs)+1;
  dr = rs+1. - k;
  k = MIN(k,790);
  ck1 = C[k-1]+dr*(C[k+1-1]-C[k-1]);
  sk1 = S[k-1]+dr*(S[k+1-1]-S[k-1]);
  a = rs*wspeed_cos -2*vmax*(yp*ck1 + xp*sk1);
  b = rs*wspeed_sin +2*vmax*(xp*ck1 - yp*sk1);
      /* Skipping ... because we don't have lake winds..
       *  rhol=fabs(xp*cos(dir*PI_180)+yp*sin(dir*PI_180)
       *  ccc=(rmax_2 + rsq)*rhol/(rmax_2+rhol*rhol)
       *  a=a+ccc*X12_23
       *  b=b+ccc*X12_24
       *******************/
  cc = rmax/(rmax_2 + rsq);
     /*
      U = (a*cc*.44704);  convert to meters/sec from MPH.
      V = (b*cc*.44704);
      */
/*  U = (a*cc/1.15077944802)*w1_w10;  * convert to knots *
    V = (b*cc/1.15077944802)*w1_w10;  * w1_w10 is for 10min->1min average winds */
  *U = (a*cc)*w1_w10;
  *V = (b*cc)*w1_w10;
/*
  mag = sqrt(U*U+V*V);
  ansList = Tcl_NewListObj(0, NULL);
  Tcl_ListObjAppendElement(interp, ansList, Tcl_NewDoubleObj(mag));
  Tcl_ListObjAppendElement(interp, ansList, Tcl_NewDoubleObj(U));
  Tcl_ListObjAppendElement(interp, ansList, Tcl_NewDoubleObj(V));
  Tcl_SetObjResult(interp, ansList);
*/

/*  sprintf (interp->result, "%8.2f ", mag); */
      /* I Don't compute Pressure field...
       * pk1 = delp[k-1] + dr*(delp[k+1-1]-delp[k-1]);
       * P[i-1][j-1] = 120-pk1*RHOWG*10+.5;
       **************************/
  /* end computation of vector winds */

  return 0;
}
