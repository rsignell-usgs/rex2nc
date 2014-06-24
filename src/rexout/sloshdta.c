/*****************************************************************************
 * sloshdta.c
 *
 * DESCRIPTION
 *    This file contains the code to read geographic data from the SLOSH model
 * bathy/topo files.
 *
 * HISTORY
 *  7/2007 Arthur Taylor (MDL): Created.
 *
 * NOTES
 ****************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "libaat.h"
#include "sloshdta.h"

#ifdef MEMWATCH
#include "memwatch.h"
#endif

/*****************************************************************************
 * read15f4() (Private) -- Arthur Taylor / MDL
 *
 * PURPOSE
 *    To convert a FORTRAN ASCII string into an array of 15 floats (F4 fomrat)
 *
 * ARGUMENTS
 * buffer = The ASCII string to read from. (Input)
 *    ray = Where to store the results. (Output)
 *
 * RETURNS: int
 * -1 = String is too short.
 *  0 = Ok
 *
 * HISTORY
 *  7/2007 Arthur Taylor (MDL): Created.
 *
 * NOTES
 ****************************************************************************/
static int read15f4(char *buffer, float ray[15])
{
   size_t i;            /* Loop counter over the array of floats. */
   char c;              /* Temporary char used to convert ASCII to float. */

   if (strlen(buffer) < 15 * 4 + 1) {
      myWarn_Err2Arg("Buffer '%s' is too short\n", buffer);
      return -1;
   }
   for (i=0; i < 15; i++) {
      c = buffer[4];
      buffer[4] = '\0';
      ray[i] = atof(buffer);
      buffer[4] = c;
      buffer += 4;
   }
   return 0;
}

/*****************************************************************************
 * readDepths() -- Arthur Taylor / MDL
 *
 * PURPOSE
 *    Read the depths from a SLOSH model bathy/topo file.
 *
 * ARGUMENTS
 *   file = The bathy/topo file to read from. (Input)
 * depths = The 2d array to store the results in. (Output)
 *   imxb = The I dimmension of the array. (Input)
 *   jmxb = The J dimmension of the array. (Input)
 *
 * RETURNS: int
 * -1 = Couldn't open the file for read.
 *  0 = Ok
 *
 * HISTORY
 *  7/2007 Arthur Taylor (MDL): Created.
 *
 * NOTES
 ****************************************************************************/
int readDepths(const char *file, float **depths, int imxb, int jmxb)
{
   FILE *fp;            /* Opened file pointer. */
   char *buffer = NULL; /* Current line from the file. */
   size_t buffLen = 0;  /* Allocated length of buffer. */
   char ac;             
   float fltRay[15];   
   int i;             
   int curI, curJ;

   if ((fp = fopen (file, "rt")) == NULL) {
      myWarn_Err2Arg("Couldn't open file %s for read", file);
      return -1;
   }

   /* Perform a grep on the file looking for the first instance of DEPTHS */
   while (reallocFGets(&buffer, &buffLen, fp) > 0) {
      ac = buffer[66];
      buffer[66] = '\0';
      if (strcmp ((buffer + 60), "DEPTHS") == 0) {
         buffer[66] = ac;
         break;
      }
      buffer[66] = ac;
   }

   /* If imxb < 15 we have a problem with the initial ray. */
   myAssert (imxb >= 15);
   if (read15f4(buffer, fltRay) != 0) {
      fclose (fp);
      free (buffer);
      return -2;
   }
   curJ = 0;
   curI = 0;
   for (i = 0; i < 15; i++) {
      depths[curI + i][curJ] = -1 * fltRay[i];
   }
   curJ++;
   while (reallocFGets(&buffer, &buffLen, fp) > 0) {
      if (read15f4(buffer, fltRay) != 0) {
         fclose (fp);
         free (buffer);
         return -2;
      }
      for (i = 0; (i < 15) && (curI + i < imxb); i++) {
         depths[curI + i][curJ] = -1 * fltRay[i];
      }
      curJ++;
      if (curJ == jmxb) {
         curI += 15;
         curJ = 0;
         if (curI >= imxb) {
            break;
         }
      }
   }
   reallocFGets(&buffer, &buffLen, fp);
   fclose (fp);
   free (buffer);
   return 0;
}

