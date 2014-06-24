#ifndef WIND_H
#define WIND_H

#include "libaat.h"

int Wind_Vmax (double pdrop, double rmax, double lat, double dir,
                      double wspeed, double *S, double *C, double *P,
                      double *delp, double *Vmax);

/* tenMinFact = Number to convert from 10 min avg winds to 1 min avg winds
 *      pdrop = Storm's Delta Pressure (MBs) (Input)
 *       rmax = Storm's Radius of max winds (st. MI) (Input)
 *        lat = Latitude of storm (clarke deg) (Input)
 *        dir = Storm's direction 0 for north. (Input)
 *      speed = Storm's forward speed (MPH) (Input)
 *       Vmax = The max 1 minute avg sustained windspeed. (Output) */
int SLOSH_Vmax (double tenMinFact, double pdrop, double rmax, double lat,
                double dir, double speed, double *Vmax);

/* tenMinFact = Number to convert from 10 min avg winds to 1 min avg winds
 *      pdrop = Storm's Delta Pressure (MBs) (Input)
 *       Vmax = The max 1 minute avg sustained windspeed. (Input)
 *        lat = Latitude of storm (clarke deg) (Input)
 *        dir = Storm's direction 0 for north. (Input)
 *      speed = Storm's forward speed (MPH) (Input)
 *       Rmax = Storm's Radius of max winds (st. MI) (Output) */
int SLOSH_Rmax (double tenMinFact, double pdrop, double Vmax, double lat,
                double dir, double speed, double *Rmax);

/* tenMinFact = Number to convert from 10 min avg winds to 1 min avg winds
 *       rmax = Storm's Radius of max winds (st. MI) (Output)
 *       Vmax = The max 1 minute avg sustained windspeed. (Input)
 *        lat = Latitude of storm (clarke deg) (Input)
 *        dir = Storm's direction 0 for north. (Input)
 *      speed = Storm's forward speed (MPH) (Input)
 *      PDrop = Storm's Delta Pressure (MBs) (Output) */
int SLOSH_pDrop (double tenMinFact, double rmax, double Vmax, double lat,
                 double dir, double speed, double *PDrop);

/* pdrop in MBs, lat/lon in negative west, rmax in stat mile,
   wspeed is storm forward motion in MPH, dir is storm direction from north,
   U, V are in MPH (magnitude is srqt (sum square) in MPH)  
 */
int windProbe (LatLon track, double pdrop, double rmax, double wspeed,
               double dir, LatLon station, double *U, double *V);

#endif
