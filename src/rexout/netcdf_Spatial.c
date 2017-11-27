/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

#ifndef HAVE_NETCDF_SUPPORT

/* Fake prototypes if we do not have NetCDF support */

static void netcdf_DatatypeInitSpatial() { printf ("NetCDF output not supported\n");}
static void netcdf_OpenSpatial(char *Filename) {}
static void netcdf_UpdateTimeUnitsSpatial(int year, int month, int day,  int hour, int min, int sec) {}
static void netcdf_SetVerticalDatumSpatial(int VerticalDatumFlag) {}
static void netcdf_ConvertUnitsSpatial(short int imxb, short int jmxb, point_type *pts) {}
static void netcdf_OutputDataSpatial(short int imxb, short int jmxb, point_type *pts) {}
static void netcdf_CloseSpatial() {}
static void netcdf_DatatypeFreeSpatial() {}
static void netcdf_WriteHeaderSpatial(short int imxb, short int jmxb) {}

#else

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Define the max function */
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

/* Structures to hold NetCDF Data */
typedef struct {

  /* file handle */
  int ncid;

  /* dimensions */
  int        time_dim;
  int        imxb_dim;
  int        jmxb_dim;

  /* ids */
  int station_name_id;
  int         time_id;
  int          lon_id;
  int          lat_id;
  int          eta_id;
  int       etamax_id;
  int            u_id;
  int            v_id;

  /* Units of Time */
  char *time_attribute_units;

  /* Vertical Datum */
  char *vertical_datum;

} NetcdfDataTypeSpatial;

/* File Local Variable to hold NetCDF Data */
NetcdfDataTypeSpatial      NetcdfDataSpatial;

/* Init the data type */
static void netcdf_DatatypeInitSpatial()
{
  NetcdfDataSpatial.                ncid = -1;
  
  NetcdfDataSpatial.            time_dim = -1;
  NetcdfDataSpatial.            imxb_dim = -1;
  NetcdfDataSpatial.            jmxb_dim = -1;
  
  NetcdfDataSpatial.             time_id = -1;
  NetcdfDataSpatial.              lon_id = -1;
  NetcdfDataSpatial.              lat_id = -1;
  NetcdfDataSpatial.              eta_id = -1;
  NetcdfDataSpatial.           etamax_id = -1;
  NetcdfDataSpatial.                u_id = -1;
  NetcdfDataSpatial.                v_id = -1;

  NetcdfDataSpatial.time_attribute_units = NULL;
  NetcdfDataSpatial.vertical_datum       = NULL;

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Free the data type */
static void netcdf_DatatypeFreeSpatial()
{
   free(NetcdfDataSpatial.time_attribute_units);
   free(NetcdfDataSpatial.vertical_datum);
}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */
static void netcdf_OpenSpatial(char *Filename)
{
  int status; /* return status variable */

  /* Open the NetCDF File */
  if (( status = nc_create(Filename, NC_CLOBBER | NC_NETCDF4, &NetcdfDataSpatial.ncid) ))
    NETCDF_ERR(status);

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Write the NetCDF header */
static void netcdf_WriteHeaderSpatial(short int imxb, short int jmxb)
{
  /* NetCDF Variables */
  int   dimids_2d[2];
  int   dimids_3d[3];
  char *attrib_value = NULL;
  int   status;
  float FillValue = NETCDF_FILLVALUE;

  /* Allocate memory for attribute value variable */
  attrib_value = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );
  
  /* NetCDF Global Attributes */
  strcpy (attrib_value, "CF-1.3");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NC_GLOBAL,"Conventions",    strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
    
  /* Define Dimensions */
  if (( status = nc_def_dim(NetcdfDataSpatial.ncid, "time", NC_UNLIMITED, &NetcdfDataSpatial.time_dim    ) ))
    NETCDF_ERR(status);
  if (( status = nc_def_dim(NetcdfDataSpatial.ncid, "number_of_i_points", imxb,         &NetcdfDataSpatial.imxb_dim    ) ))
    NETCDF_ERR(status);
  if (( status = nc_def_dim(NetcdfDataSpatial.ncid, "number_of_j_points", jmxb,         &NetcdfDataSpatial.jmxb_dim    ) ))
    NETCDF_ERR(status);
  
  /* Define Variables */
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "time", NC_DOUBLE, 1, &NetcdfDataSpatial.time_dim, &NetcdfDataSpatial.time_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Time");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "gregorian");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,"calendar",      strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  dimids_2d[0]=NetcdfDataSpatial.jmxb_dim;
  dimids_2d[1]=NetcdfDataSpatial.imxb_dim;
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "lon", NC_DOUBLE, 2, dimids_2d, &NetcdfDataSpatial.lon_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "longitude");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lon_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Longitude");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lon_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "degrees_east");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lon_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  dimids_2d[0]=NetcdfDataSpatial.jmxb_dim;
  dimids_2d[1]=NetcdfDataSpatial.imxb_dim;
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "lat", NC_DOUBLE, 2, dimids_2d, &NetcdfDataSpatial.lat_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "latitude");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lat_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Latitude");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lat_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "degrees_north");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lat_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);

  dimids_2d[0]=NetcdfDataSpatial.jmxb_dim;
  dimids_2d[1]=NetcdfDataSpatial.imxb_dim;
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "etamax", NC_FLOAT, 2, dimids_2d, &NetcdfDataSpatial.etamax_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "water_surface_height_above_reference_datum");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Water Surface Height Above Reference Datum");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time: maximum");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"cell_methods",  strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "lat lon");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);

  dimids_3d[0]=NetcdfDataSpatial.time_dim;
  dimids_3d[1]=NetcdfDataSpatial.jmxb_dim;
  dimids_3d[2]=NetcdfDataSpatial.imxb_dim;
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "eta", NC_FLOAT, 3, dimids_3d, &NetcdfDataSpatial.eta_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time: point");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"cell_methods",  strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time lat lon");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "water_surface_height_above_reference_datum");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
     NETCDF_ERR(status);
  strcpy (attrib_value, "Water Surface Height Above Reference Datum");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  
  dimids_3d[0]=NetcdfDataSpatial.time_dim;
  dimids_3d[1]=NetcdfDataSpatial.jmxb_dim;
  dimids_3d[2]=NetcdfDataSpatial.imxb_dim;
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "u", NC_FLOAT, 3, dimids_3d, &NetcdfDataSpatial.u_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "eastward_wind");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.u_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Eastward Wind");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.u_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m/s");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.u_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time lat lon");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.u_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataSpatial.ncid,NetcdfDataSpatial.u_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);
  
  dimids_3d[0]=NetcdfDataSpatial.time_dim;
  dimids_3d[1]=NetcdfDataSpatial.jmxb_dim;
  dimids_3d[2]=NetcdfDataSpatial.imxb_dim;
  if (( status=nc_def_var(NetcdfDataSpatial.ncid, "v", NC_FLOAT, 3, dimids_3d, &NetcdfDataSpatial.v_id ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "northward_wind");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.v_id,"standard_name", strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "Northward Wind");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.v_id,"long_name",     strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "m/s");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.v_id,"units",         strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  strcpy (attrib_value, "time lat lon");
  if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.v_id,"coordinates",   strlen(attrib_value),  attrib_value  ) ))
    NETCDF_ERR(status);
  if (( status =  nc_put_att_float(NetcdfDataSpatial.ncid,NetcdfDataSpatial.v_id,"_FillValue", NC_FLOAT, 1, &FillValue ) ))
    NETCDF_ERR(status);
  
#ifdef HAVE_HDF_SUPPORT
  /* Setup deflate */
  if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lon_id,          NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.lat_id,          NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  
  /* For some reason if time is deflated, the CF compliance checkers die
     maybe because it is the unlimited variable.
     if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,         NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
     NETCDF_ERR(status);
  */

  if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,       NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,          NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.u_id,            NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
  if (( status =  nc_def_var_deflate(NetcdfDataSpatial.ncid,NetcdfDataSpatial.v_id,            NETCDF_SHUFFLE_FLAG, NETCDF_DEFLATE_FLAG, NETCDF_DEFLATE_LEVEL ) ))
    NETCDF_ERR(status);
#endif

  /* Final attribute to be define so can free memory */
  free (attrib_value);
  
  /* Leave "define" mode */
  if (( status = nc_enddef(NetcdfDataSpatial.ncid) ))
    NETCDF_ERR(status);
  
}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Update the units of Time */
static void netcdf_SetVerticalDatumSpatial(int VerticalDatumFlag)
{

  int status;

  /* Define the units of time */
  if ( NetcdfDataSpatial.vertical_datum == NULL ) {
    
    /* Leave "define" mode */
    if (( status = nc_redef(NetcdfDataSpatial.ncid) ))
      NETCDF_ERR(status);
    
    /* Allocate memory for the units variable */
    NetcdfDataSpatial.vertical_datum = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );

    /* Set the correct vertical datum */
    switch (VerticalDatumFlag) {

    case VERTICAL_DATUM_NGVD29_VAL:
      /* Define the NGVD29 vertical Datum */
      sprintf (NetcdfDataSpatial.vertical_datum, "urn:ogc:def:datum:epsg::5102");     
      if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"VerticalDatum", strlen(NetcdfDataSpatial.vertical_datum), NetcdfDataSpatial.vertical_datum ) ))
	NETCDF_ERR(status);
      if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"VerticalDatum", strlen(NetcdfDataSpatial.vertical_datum), NetcdfDataSpatial.vertical_datum ) ))
	NETCDF_ERR(status);
      break;

    case VERTICAL_DATUM_NAVD88_VAL:
      /* Define the NAVD88 vertical Datum */
      sprintf (NetcdfDataSpatial.vertical_datum, "urn:ogc:def:datum:epsg::5103");
      if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.eta_id,"VerticalDatum", strlen(NetcdfDataSpatial.vertical_datum), NetcdfDataSpatial.vertical_datum ) ))
	NETCDF_ERR(status);
      if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.etamax_id,"VerticalDatum", strlen(NetcdfDataSpatial.vertical_datum), NetcdfDataSpatial.vertical_datum ) ))
	NETCDF_ERR(status);
      break;
     
    default:
      printf ("Unknown vertical datum flag value: %i\n",VerticalDatumFlag);
      printf ("  Not setting the VerticalDatum attribute\n") ;
      break;

 
    }

    /* Leave "define" mode */
    if (( status = nc_enddef(NetcdfDataSpatial.ncid) ))
      NETCDF_ERR(status);
  
  }  

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Update the units of Time */
static void netcdf_UpdateTimeUnitsSpatial(int year, int month, int day,  int hour, int min, int sec)
{

  /* NetCDF Variables */
  int   status;

  /* Define the units of time */
  if ( NetcdfDataSpatial.time_attribute_units == NULL ) {
    
    /* Leave "define" mode */
    if (( status = nc_redef(NetcdfDataSpatial.ncid) ))
      NETCDF_ERR(status);
    
    /* Allocate memory for the units variable */
    NetcdfDataSpatial.time_attribute_units = (char *) malloc (NETCDF_NAME_STRLEN * sizeof (char) );

    /* Define the time */
    sprintf (NetcdfDataSpatial.time_attribute_units, "minutes since %04d-%02d-%02d %02d:%02d:%02d", year, month, day,  hour, min, sec);
    if (( status =  nc_put_att_text(NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,"units", strlen(NetcdfDataSpatial.time_attribute_units), NetcdfDataSpatial.time_attribute_units ) ))
      NETCDF_ERR(status);
    
    /* Leave "define" mode */
    if (( status = nc_enddef(NetcdfDataSpatial.ncid) ))
      NETCDF_ERR(status);
  
  }  

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

static void netcdf_ConvertUnitsSpatial(short int imxb, short int jmxb, point_type *pts)
{
  /* 
     Convert U/V from MPH to m/s  and surge from ft to m
  */
  int t, num_pts, n;
  
  /* Total number of points */
  num_pts = imxb*jmxb;

  /* Loop through all the points */
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

/* Output Data */
static void netcdf_OutputDataSpatial(short int imxb, short int jmxb, point_type *pts)
{
  
  /* NetCDF Variables */
  float   *array_float_1d  = NULL;
  double  *array_double_1d = NULL;
  size_t time_start[1], etauv_start[3], lonlat_start[2];
  size_t time_count[1], etauv_count[3], lonlat_count[2];
  
  /* Generic variables */
  int n, t, status, num_pts, num_times, base_index;

  /* Calculate the total number of points */
  num_pts=imxb*jmxb;

  /* Output the times */

  /* Assume all stations have the same time array */
  base_index=0;
  num_times = pts[base_index].numSInfo;

  /* ********************************************************************* */

  /* Allocate memory for the 1D double arrays - Time */
  array_double_1d = (double *) realloc (array_double_1d, num_times * sizeof (double) );
  
  for (t = 0; t < num_times; t++) {
    array_double_1d[t] = (pts[base_index].sInfo[t].clock - pts[base_index].sInfo[0].clock) / 60.;
  }
  
  time_start[0]=0;
  time_count[0]=num_times;
  if (( status =  nc_put_vara_double (NetcdfDataSpatial.ncid,NetcdfDataSpatial.time_id,time_start, time_count, array_double_1d) ))
    NETCDF_ERR(status);
  
  /* ********************************************************************* */

  /* Allocate memory for the 2D double arrays - Lat-Lon*/
  array_double_1d = (double *) realloc (array_double_1d, imxb*jmxb * sizeof (double) );
  
  /* Setup the lon/lat/etamax start and count values */
  lonlat_start[0]=0;
  lonlat_start[1]=0;
  lonlat_count[0]=jmxb;
  lonlat_count[1]=imxb;
  
  /* Same the lon values*/
  for (n = 0; n < num_pts; n++) {
    array_double_1d[n]=pts[n].station.lon;
  }
  if (( status =  nc_put_vara_double(NetcdfDataSpatial.ncid, 
                                     NetcdfDataSpatial.lon_id,
                                     lonlat_start,
                                     lonlat_count,
                                     array_double_1d) ))
    NETCDF_ERR(status);

  /* Same the lat values*/
  for (n = 0; n < num_pts; n++) {
    array_double_1d[n]=pts[n].station.lat;
  }
  if (( status =  nc_put_vara_double(NetcdfDataSpatial.ncid, 
                                     NetcdfDataSpatial.lat_id,
                                     lonlat_start,
                                     lonlat_count,
                                     array_double_1d) ))
    NETCDF_ERR(status);

  /* free up memory */
  free (array_double_1d);

  /* ********************************************************************* */

  /* Output water level, U and V */
  
  /* Allocate memory for the 2D float arrays - eta/u/v */
  array_float_1d = (float *) realloc (array_float_1d, imxb*jmxb * sizeof (float) );

  /* Store and output surge */
  for (t = 0; t < num_times; t++) {

    etauv_start[0]=t;
    etauv_start[1]=0;
    etauv_start[2]=0;

    etauv_count[0]=1;
    etauv_count[1]=jmxb;
    etauv_count[2]=imxb;
    
    for (n = 0; n < num_pts; n++) {
      array_float_1d[n]=pts[n].sInfo[t].surge;

      /* Set to fill value if necessary */
      if ( (array_float_1d[n] >  NETCDF_MAXVALUE) ||
	   (array_float_1d[n] < -NETCDF_MAXVALUE) )
	{   
	  array_float_1d[n] = NETCDF_FILLVALUE;
	    }

    }
    if (( status =  nc_put_vara_float (NetcdfDataSpatial.ncid,
				       NetcdfDataSpatial.eta_id,
				       etauv_start,
				       etauv_count,
				       array_float_1d) ))
      NETCDF_ERR(status);
  }

  /* Store and output U */
  for (t = 0; t < num_times; t++) {

    etauv_start[0]=t;
    etauv_start[1]=0;
    etauv_start[2]=0;

    etauv_count[0]=1;
    etauv_count[1]=jmxb;
    etauv_count[2]=imxb;
    
    for (n = 0; n < num_pts; n++) {
      array_float_1d[n]=pts[n].sInfo[t].U;
    }
    if (( status =  nc_put_vara_float (NetcdfDataSpatial.ncid,
				       NetcdfDataSpatial.u_id,
				       etauv_start,
				       etauv_count,
				       array_float_1d) ))
      NETCDF_ERR(status);
  }

  /* Store and output V */
  for (t = 0; t < num_times; t++) {

    etauv_start[0]=t;
    etauv_start[1]=0;
    etauv_start[2]=0;

    etauv_count[0]=1;
    etauv_count[1]=jmxb;
    etauv_count[2]=imxb;
    
    for (n = 0; n < num_pts; n++) {
      array_float_1d[n]=pts[n].sInfo[t].V;
    }
    if (( status =  nc_put_vara_float (NetcdfDataSpatial.ncid,
				       NetcdfDataSpatial.v_id,
				       etauv_start,
				       etauv_count,
				       array_float_1d) ))
      NETCDF_ERR(status);
  }

  
  /* ********************************************************************* */

  /* Compute etamax */
  for (n = 0; n < num_pts; n++) {
    array_float_1d[n]=NETCDF_FILLVALUE;
  }

  for (t = 0; t < num_times; t++) {
    for (n = 0; n < num_pts; n++) {

      if ( (pts[n].sInfo[t].surge > -NETCDF_MAXVALUE) &&
	   (pts[n].sInfo[t].surge <  NETCDF_MAXVALUE) )
	{   
	  array_float_1d[n]=max(array_float_1d[n],pts[n].sInfo[t].surge);
	}
    }
  }	    


  /* output it */
  if (( status =  nc_put_vara_float (NetcdfDataSpatial.ncid, 
                                     NetcdfDataSpatial.etamax_id,
                                     lonlat_start,
                                     lonlat_count,
                                     array_float_1d) ))
    NETCDF_ERR(status);
  
  /* ********************************************************************* */

  /* free up memory */
  free (array_float_1d);
  
}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

/* Close the NetCDF File */
static void netcdf_CloseSpatial()
{
  int status; /* return status variable */

  /* Close the NetCDF File */
  if (( status = nc_close(NetcdfDataSpatial.ncid) ))
    NETCDF_ERR(status);

}

/* ****************************************************************************************** */
/* ****************************************************************************************** */
/* ****************************************************************************************** */

#endif
