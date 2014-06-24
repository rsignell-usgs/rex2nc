/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

#ifndef HAVE_NETCDF_SUPPORT

/* Fake prototypes if we do not have NetCDF support */

static void netcdf_DatatypeInitTS() { printf ("NetCDF output not supported\n");}
static void netcdf_OpenTS(char *Filename) {}
static void netcdf_WriteHeaderTS(int num_pts) {}
static void netcdf_UpdateTimeUnitsTS(int year, int month, int day,  int hour, int min, int sec) {}
static void netcdf_SetVerticalDatumTS(int VerticalDatumFlag) {}
static void netcdf_OutputDataTS(int num_pts, point_type *pts) {}
static void netcdf_CloseTS() {}
static void netcdf_DatatypeFreeTS() {}
static void netcdf_ConvertUnitsTS() {}

#else

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Structures to hold NetCDF Data */
typedef struct {

  /* file handle */
  int ncid;

  /* dimensions */
  int        time_dim;
  int     num_pts_dim;
  int name_strlen_dim;

  /* ids */
  int station_name_id;
  int         time_id;
  int          lon_id;
  int          lat_id;
  int          eta_id;
  int            u_id;
  int            v_id;

  /* Units of Time */
  char *time_attribute_units;

  /* Vertical Datum */
  char *vertical_datum;

} NetcdfDataTypeTS;

/* File Local Variable to hold NetCDF Data */
NetcdfDataTypeTS      NetcdfDataTS;

/* Init the data type */
static void netcdf_DatatypeInitTS()
{
  NetcdfDataTS.                ncid = -1;
  
  NetcdfDataTS.            time_dim = -1;
  NetcdfDataTS.         num_pts_dim = -1;
  NetcdfDataTS.     name_strlen_dim = -1;
  
  NetcdfDataTS.     station_name_id = -1;
  NetcdfDataTS.             time_id = -1;
  NetcdfDataTS.              lon_id = -1;
  NetcdfDataTS.              lat_id = -1;
  NetcdfDataTS.              eta_id = -1;
  NetcdfDataTS.                u_id = -1;
  NetcdfDataTS.                v_id = -1;

  NetcdfDataTS.time_attribute_units = NULL;
  NetcdfDataTS.vertical_datum       = NULL;

}
  
/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Free the data type */
static void netcdf_DatatypeFreeTS()
{
   free(NetcdfDataTS.time_attribute_units);
   free(NetcdfDataTS.vertical_datum);
}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

static void netcdf_OpenTS(char *Filename)
{
  int status; /* return status variable */

  /* Open the NetCDF File */
  if (( status = nc_create(Filename, NC_CLOBBER | NC_NETCDF4, &NetcdfDataTS.ncid) ))
    NETCDF_ERR(status);

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Close the NetCDF File */
static void netcdf_CloseTS()
{
  int status; /* return status variable */

  /* Close the NetCDF File */
  if (( status = nc_close(NetcdfDataTS.ncid) ))
    NETCDF_ERR(status);

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Write the NetCDF header */
static void netcdf_WriteHeaderTS(int num_pts)
{
  /* NetCDF Variables */
  int   dimids_2d[2];
  char *attrib_value = NULL;
  int   status;
  float FillValue = NETCDF_FILLVALUE;

  /* Allocate memory for attribute value variable */
  attrib_value = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );
  
  /* NetCDF Global Attributes */
  strcpy (attrib_value, "timeSeries");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NC_GLOBAL,"CF:featureType", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "CF-1.3");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NC_GLOBAL,"Conventions",    strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
    
  /* Define Dimensions */
  if (( status = nc_def_dim(NetcdfDataTS.ncid, "time",               NC_UNLIMITED, &NetcdfDataTS.time_dim       ) ))
    NETCDF_ERR(status);
  if (( status = nc_def_dim(NetcdfDataTS.ncid, "number_of_stations", num_pts,      &NetcdfDataTS.num_pts_dim    ) ))
    NETCDF_ERR(status);
  if (( status = nc_def_dim(NetcdfDataTS.ncid, "name_strlen", NETCDF_NAME_STRLEN,  &NetcdfDataTS.name_strlen_dim) ))
    NETCDF_ERR(status);
  
  /* Define Variables */
  dimids_2d[0]=NetcdfDataTS.num_pts_dim;
  dimids_2d[1]=NetcdfDataTS.name_strlen_dim;
  if (( status=nc_def_var(NetcdfDataTS.ncid, "station_name", NC_CHAR, 2, dimids_2d, &NetcdfDataTS.station_name_id ) ))
    NETCDF_ERR(status);
  /* there is no standard_name for this variable */
  /* strcpy (attrib_value, "station_id"); */
  /* if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.station_name_id,"standard_name", strlen(attrib_value),  attrib_value  ) )) */
  /*    NETCDF_ERR(status);*/
  strcpy (attrib_value, "Station Name");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.station_name_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  if (( status=nc_def_var(NetcdfDataTS.ncid, "time", NC_DOUBLE, 1, &NetcdfDataTS.time_dim, &NetcdfDataTS.time_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.time_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Time");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.time_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "gregorian");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.time_id,"calendar",      strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.time_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  if (( status=nc_def_var(NetcdfDataTS.ncid, "lon", NC_DOUBLE, 1, &NetcdfDataTS.num_pts_dim, &NetcdfDataTS.lon_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "longitude");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.lon_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Longitude");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.lon_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "degrees_east");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.lon_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  if (( status=nc_def_var(NetcdfDataTS.ncid, "lat", NC_DOUBLE, 1, &NetcdfDataTS.num_pts_dim, &NetcdfDataTS.lat_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "latitude");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.lat_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Latitude");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.lat_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "degrees_north");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.lat_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  dimids_2d[0]=NetcdfDataTS.num_pts_dim;
  dimids_2d[1]=NetcdfDataTS.time_dim;
  if (( status=nc_def_var(NetcdfDataTS.ncid, "eta", NC_FLOAT, 2, dimids_2d, &NetcdfDataTS.eta_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time lat lon");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "water_surface_height_above_reference_datum");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
     NETCDF_ERR(status);
  strcpy (attrib_value, "Water Surface Height Above Reference Datum");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  dimids_2d[0]=NetcdfDataTS.num_pts_dim;
  dimids_2d[1]=NetcdfDataTS.time_dim;
  if (( status=nc_def_var(NetcdfDataTS.ncid, "u", NC_FLOAT, 2, dimids_2d, &NetcdfDataTS.u_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "eastward_wind");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.u_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Eastward Wind");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.u_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m/s");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.u_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time lat lon");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.u_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataTS.ncid,NetcdfDataTS.u_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);
  
  dimids_2d[0]=NetcdfDataTS.num_pts_dim;
  dimids_2d[1]=NetcdfDataTS.time_dim;
  if (( status=nc_def_var(NetcdfDataTS.ncid, "v", NC_FLOAT, 2, dimids_2d, &NetcdfDataTS.v_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "northward_wind");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.v_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Northward Wind");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.v_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m/s");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.v_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time lat lon");
  if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.v_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataTS.ncid,NetcdfDataTS.v_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);
  
#ifdef HAVE_HDF_SUPPORT
  /* Setup deflate */
  if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.lon_id,          NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.lat_id,          NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.station_name_id, NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  
  /* For some reason if time is deflated, the CF compliance checkers die
     maybe because it is the unlimited variable.
     if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.time_id,         NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
     NETCDF_ERR(status);
  */

  if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,          NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.u_id,            NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataTS.ncid,NetcdfDataTS.v_id,            NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
#endif

  /* Final attribute to be define so can free memory */
  free (attrib_value);
  
  /* Leave "define" mode */
  if (( status = nc_enddef(NetcdfDataTS.ncid) ))
    NETCDF_ERR(status);
  
}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Update the units of Time */
static void netcdf_SetVerticalDatumTS(int VerticalDatumFlag)
{

  int status;

  /* Define the units of time */
  if ( NetcdfDataTS.vertical_datum == NULL ) {
    
    /* Leave "define" mode */
    if (( status = nc_redef(NetcdfDataTS.ncid) ))
      NETCDF_ERR(status);
    
    /* Allocate memory for the units variable */
    NetcdfDataTS.vertical_datum = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );

    /* Set the correct vertical datum */
    switch (VerticalDatumFlag) {

    case VERTICAL_DATUM_NGVD29_VAL:
      /* Define the NGVD29 vertical Datum */
      sprintf (NetcdfDataTS.vertical_datum, "urn:ogc:def:datum:epsg::5102");
      if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"VerticalDatum", strlen(NetcdfDataTS.vertical_datum), NetcdfDataTS.vertical_datum ) ))
	NETCDF_ERR(status);
      break;

    case VERTICAL_DATUM_NAVD88_VAL:
      /* Define the NAVD88 vertical Datum */
      sprintf (NetcdfDataTS.vertical_datum, "urn:ogc:def:datum:epsg::5103");
      if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.eta_id,"VerticalDatum", strlen(NetcdfDataTS.vertical_datum), NetcdfDataTS.vertical_datum ) ))
	NETCDF_ERR(status);
      break;
      
    default:
      printf ("Unknown vertical datum flag value: %i\n",VerticalDatumFlag);
      printf ("  Not setting the VerticalDatum attribute\n") ;
      break;
      
    }
    
    /* Leave "define" mode */
    if (( status = nc_enddef(NetcdfDataTS.ncid) ))
      NETCDF_ERR(status);
  
  }  

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Update the units of Time */
static void netcdf_UpdateTimeUnitsTS(int year, int month, int day,  int hour, int min, int sec)
{

  /* NetCDF Variables */
  int   status;

  /* Define the units of time */
  if ( NetcdfDataTS.time_attribute_units == NULL ) {
    
    /* Leave "define" mode */
    if (( status = nc_redef(NetcdfDataTS.ncid) ))
      NETCDF_ERR(status);
    
    /* Allocate memory for the units variable */
    NetcdfDataTS.time_attribute_units = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );

    /* Define the time */
    sprintf (NetcdfDataTS.time_attribute_units, "minutes since %04d-%02d-%02d %02d:%02d:%02d", year, month, day,  hour, min, sec);
    if (( status =  nc_put_att_text(NetcdfDataTS.ncid,NetcdfDataTS.time_id,"units", strlen(NetcdfDataTS.time_attribute_units), NetcdfDataTS.time_attribute_units ) ))
      NETCDF_ERR(status);
    
    /* Leave "define" mode */
    if (( status = nc_enddef(NetcdfDataTS.ncid) ))
      NETCDF_ERR(status);
  
  }  

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Output Data */
static void netcdf_OutputDataTS(int num_pts, point_type *pts)
{
  
  /* NetCDF Variables */
  char    *station_name    = NULL;
  float   *array_float_1d  = NULL;
  double  *array_double_1d = NULL;
  size_t time_start[1], station_name_start[2], etauv_start[2];
  size_t time_count[1], station_name_count[2], etauv_count[2];
  
  /* Generic variables */
  int n, t, status, num_times, base_index;

  /* Output the times */

  /* Assume all stations have the same time array */
  base_index=0;
  num_times = pts[base_index].numSInfo;

  /* Allocate memory for the 1D double arrays - Time */
  array_double_1d = (double *) realloc (array_double_1d, num_times * sizeof (double) );
  
  for (t = 0; t < num_times; t++) {
    array_double_1d[t] = (pts[base_index].sInfo[t].clock - pts[base_index].sInfo[0].clock) / 60.;
  }
  
  time_start[0]=0;
  time_count[0]=num_times;
  if (( status =  nc_put_vara_double (NetcdfDataTS.ncid,NetcdfDataTS.time_id,time_start, time_count, array_double_1d) ))
    NETCDF_ERR(status);
  
  /* Allocate memory for the 1D double arrays - Lat-Lon*/
  array_double_1d = (double *) realloc (array_double_1d, num_pts * sizeof (double) );
  
  for (n = 0; n < num_pts; n++) {
    array_double_1d[n]=pts[n].station.lon;
  }
  if (( status =  nc_put_var_double(NetcdfDataTS.ncid,NetcdfDataTS.lon_id,array_double_1d) ))
    NETCDF_ERR(status);
  
  for (n = 0; n < num_pts; n++) {
    array_double_1d[n]=pts[n].station.lat;
  }
  if (( status =  nc_put_var_double(NetcdfDataTS.ncid,NetcdfDataTS.lat_id,array_double_1d) ))
    NETCDF_ERR(status);
  
  /* Free up memory */
  free (array_double_1d);
  
  /* Allocate memory for the station names */
  station_name = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );
  
  /* Create/Output the station names (i,j) */
  for (n = 0; n < num_pts; n++) {
    sprintf (station_name, "(%d,%d)", pts[n].x, pts[n].y);
    
    station_name_start[0]=n;
    station_name_start[1]=0;
    station_name_count[0]=1;
    station_name_count[1]=strlen(station_name)+1;
    
    if (( status =  nc_put_vara_text(NetcdfDataTS.ncid,NetcdfDataTS.station_name_id, station_name_start, station_name_count, station_name) ))
      NETCDF_ERR(status);
    
  }
  
  /* Free up memory */
  free (station_name);
  
  /* Output water level */
  
  /* Allocate memory for the surge */
  array_float_1d = (float *) malloc ( num_times * sizeof (float) );

  /* Store and output surge */
  for (n = 0; n < num_pts; n++) {

    etauv_start[0]=n;
    etauv_start[1]=0;
    etauv_count[0]=1;
    etauv_count[1]=num_times;
    
    for (t = 0; t < num_times; t++) {
      array_float_1d[t]=pts[n].sInfo[t].surge;
    }
    if (( status =  nc_put_vara_float (NetcdfDataTS.ncid,NetcdfDataTS.eta_id,etauv_start, etauv_count, array_float_1d) ))
      NETCDF_ERR(status);
  }
  
  /* free up memory */
  free (array_float_1d);
  
  /* Allocate memory for U/V */
  array_double_1d = (double *) malloc ( num_times * sizeof (double) );
  
  /* Store and output surge */
  for (n = 0; n < num_pts; n++) {
    
    etauv_start[0]=n;
    etauv_start[1]=0;
    etauv_count[0]=1;
    etauv_count[1]=num_times;
    
    for (t = 0; t < pts[n].numSInfo; t++) {
      array_double_1d[t]=pts[n].sInfo[t].U;
    }
    if (( status =  nc_put_vara_double (NetcdfDataTS.ncid,NetcdfDataTS.u_id,etauv_start, etauv_count, array_double_1d) ))
      NETCDF_ERR(status);
    
    for (t = 0; t < pts[n].numSInfo; t++) {
      array_double_1d[t]=pts[n].sInfo[t].V;
    }
    if (( status =  nc_put_vara_double (NetcdfDataTS.ncid,NetcdfDataTS.v_id,etauv_start, etauv_count, array_double_1d) ))
      NETCDF_ERR(status);
    
  }
  
  /* free up memory */
  free (array_double_1d);

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

static void netcdf_ConvertUnitsTS(int num_pts, point_type *pts)
{
  /* 
     Convert U/V from MPH to m/s  and surge from ft to m
  */
  int n, t;
  
  for (n = 0; n < num_pts; n++) {
    for (t = 0; t < pts[n].numSInfo; t++) {
      pts[n].sInfo[t].U     = pts[n].sInfo[t].U     * 0.44704;
      pts[n].sInfo[t].V     = pts[n].sInfo[t].V     * 0.44704;
      pts[n].sInfo[t].surge = pts[n].sInfo[t].surge * 0.3048;
    }
  }
   
}
 
/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

#endif
