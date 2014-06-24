/* Define the Datum Flags here */
#define VERTICAL_DATUM_NGVD29_VAL 0
#define VERTICAL_DATUM_NAVD88_VAL 9999

#ifdef HAVE_NETCDF_SUPPORT

/* Main NetCDF Headers */
#include <netcdf.h>

/* Handle errors by printing an error message and exiting with a non-zero status */
#define NETCDF_ERRCODE 2
#define NETCDF_ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(NETCDF_ERRCODE);}

/* NetCDF Constants */
#ifndef NETCDF_NAME_STRLEN
#define NETCDF_NAME_STRLEN 64 /* Max string length of station names and attributes */
#endif

/* Max/min value in meters - should be slightly less than 99.9*0.3048 */
#ifndef NETCDF_MAXVALUE
#define NETCDF_MAXVALUE 25.0
#endif

#ifndef NETCDF_FILLVALUE
#define NETCDF_FILLVALUE -999.9
#endif

#ifndef NETCDF_SHUFFLE_FLAG
#define NETCDF_SHUFFLE_FLAG 0
#endif

#ifndef NETCDF_DEFLATE_FLAG
#define NETCDF_DEFLATE_FLAG 1
#endif

#ifndef NETCDF_DEFLATE_LEVEL
#define NETCDF_DEFLATE_LEVEL 4
#endif

#endif
