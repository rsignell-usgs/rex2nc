#ifndef _TIO3_H
#define _TIO3_H

/*
 * tio3.h --
 *
 * Implements binary file I/O in an Endian independent manner (ie specify if
 * the file was created on an intel machine or opposite to intel, when you
 * open the file, and then ignore), and in such a way that it can be called
 * from FORTRAN 77.
 *
 *   Note: If using FORTRAN generated files, they sometimes stick extra bytes
 * at the beginning and end of a record.  These are "hidden" from the FORTRAN
 * user, but not from this library.  They usually contain the record length.
 * Example: in UNIX it usually stores an int*4, and in Microsoft Fortran, an
 * int *1 with 129 being a continuation flag.
 *
 * Copyright (c) 1997 RDC / US Dept of Commerce-NOAA-NWS
 * Copyright (c) 2001 RSIS / US Dept of Commerce-NOAA-NWS
 *
 * Author: Arthur Taylor (arthur.taylor@noaa.gov),
 *    Marine Group, Evaluations Branch, Meteorological Development Lab,
 *    U.S. National Weather Service.
 *
 * Procedure list in order of which they are defined:
 *   1) tOpen
 *   2) tFind
 *   3) tClose
 *   4) tEndian_Switch
 *   5) tGet_Endian
 *   6) tFlush
 *   7) tSeek
 *   8) tGetBit
 *   9) tRead
 *  10) tReadStruct
 *  11) tPutBit
 *  12) tWrite
 */

#include <stdio.h>

typedef enum { TFLAG_READ, TFLAG_WRITE, TFLAG_APPEND, TFLAG_RW } TFLAG_ACCESS;
typedef enum { TFLAG_MadeOnIntel, TFLAG_NotMadeOnIntel } TFLAG_SYSTEM;
/* MadeOnIntel    ==> LittleEndian
 * NotMadeOnIntel ==> BigEndian
 */

#ifdef WORDS_BIGENDIAN
  #undef _sysINTEL
#else
  #define _sysINTEL
#endif
/* Don't need to worry about _sysDOS since that handles fortran wrappers,
 * and we don't have fortran wrappers in rexout.  We may want to investigate
 * AC_F77_WRAPPERS
 */
/*
#ifdef _UNIX_
#else
  #ifdef _LINUX_
    #define _sysINTEL
  #else
    #define _sysDOS
    #define _sysINTEL
  #endif
#endif
*/

#ifndef BYTE
  typedef unsigned char BYTE;
/*  #define BYTE unsigned char*/
#endif

#ifndef BOOL
  typedef int BOOL;
/*  #define BOOL unsigned char */
#endif

/* Also have an _ERROR_CHECK_ flag? */

/* Idea behind stream is that we can load a buffer of data from the file,
   store it in stream, and then parse it out using tRead / tGetBit
   operations... */
typedef struct {
  FILE *fp;
  char *name;
  unsigned short int fid;
  TFLAG_ACCESS access;
  BOOL f_os; /* true if file is opposite endian from operating system. */
  BYTE gbuf, pbuf;
  unsigned char gbuf_loc, pbuf_loc; /* These are which bit from 1..8 we
           are looking at.  If 0, then gbuf is empty, but pbuf is full. */

  BOOL f_stream; /* false means use fp, true means use stream */
  BYTE *stream;
  unsigned short int stream_len, stream_loc;
} TIO_type;

/* -------------C routines, and calling prototypes...-------------- */
/* NULL on error, valid structure on success. */
TIO_type *tOpen (unsigned short int fid, char *name, TFLAG_ACCESS access,
                 TFLAG_SYSTEM sys);

/* returns NULL if it couldn't find this fid. */
TIO_type *tFind (unsigned short int fid);

/* Returns 1 "true" on error. */
BOOL tClose (TIO_type *tio);

/* following 3 could be done directly by adjusting tio.*/
void tEndian_Switch (TIO_type *tio);

/* true if file is opposite endian from operating system. */
BOOL tGet_Endian (TIO_type *tio);

int tFlush (BOOL f_put, TIO_type *tio);

/* Simply a call to fseek.  (returns non-zero on error.)
   origin == SEEK_SET, SEEK_CUR, SEEK_END */
int tSeek (TIO_type *tio, long int offset, int origin);

/* -----Read Routines-----*/

/* If EOF, returns EOF (-1).  On error 1 "true".  Normally: 0->"false". */
/* dst is usually NOT an array of elements. */
signed char tGetBit (void *dst, size_t dst_len, unsigned short int num_bits,
                     TIO_type *tio);

/* Returns number of elements read.  See feof and ferror, for EOF/errors. */
/* dst can be an array of elements each of "elem_size" length */
size_t tRead (void *Dst, size_t elem_size, size_t num_elem, TIO_type *tio);

/* Note: Because of the organization of tRead, one can do the following...
  If one knows that their file was written say on intel...
  #ifdef _sysINTEL
    FILE *tio;
    #define tRead fread
  #else
    TIO_type *tio;
  #endif
  Which would speed it up a little on the "same" endian machine.
*/

/* To read a 3 byte int, try:
  long int li_temp = 0;
  #ifdef _sysINTEL
    tRead (& li_temp, 3, 1, tio);
  #else
    tRead ( (((char *)&li_temp)+1), 3, 1, tio);
  #endif
  ... alternatively ...
  union {long int li_val, char s_val[4]} val;
  val.li_val = 0;
  #ifdef _sysINTEL
    tRead (& val.li_temp, 3, 1, tio);
  #else
    tRead (& (val.s_val[1]), 3, 1, tio);
  #endif
*/

/* ---tReadStruct---*/
/* dst_len is passed so that on "same" endian machine this is a simple call
   to fread.  On "other" endian machine, it uses elem_sizes to figure out
   how to swap bytes.  An example of elem_sizes is "1,4,2,8" => 1 char,
   1 long, 1 short, 1 double  Types larger than 9 should be A..Z
   (ascii 65..90). Allows "4*20" => 20 long ints. So you might have
   "2,G,1,1*20,2".  There is no "0*20" as the "same" endian machine wouldn't
   see it, so if I skipped either in dst or the file, it would only work on
   the "other" endian machine. */
/* If you have pointers, put them at the end, and decrease dst_len, or put
   them at the beginning and adjust both Dst and dst_len. */
/* For 3 byte ints stored as longs, put them at front or back (same as
   pointers), then use tRead to fill them. */
/* Returns number of elements read. See feof and ferror, for EOF/errors. */
BOOL tReadStruct (void *Dst, size_t dst_len, const char *elem_sizes,
                  TIO_type *tio);

/* -----Write Routines-----*/

/* On error 1 "true".  Normally: 0->"false". */
BOOL tPutBit (void *Src, size_t src_len, unsigned short int num_bits,
              TIO_type *tio);

/* Returns number of elements written.  See feof and ferror. */
/* src can be an array of elements each of "elem_size" length */
size_t tWrite (void *src, size_t elem_size, size_t num_elem, TIO_type *tio);
/* see notes about tRead... they should apply here as well. */


/* -------------C-called FORTRAN helper routine...------ */
int tSet (TFLAG_SYSTEM sys);

/* -------------FORTRAN called routines...-------------- */
/*
#ifdef _sysDOS
  #define F2C
#else
  #define HP
#endif

#ifdef F2C
  void topen_ (int *fid, char *name, int *flag, int *sys);
  void tread_ (int *fid, long int *ptr, int *num);
  void treadc_ (int *fid, char *ptr, int *num);
  void treadf_ (int *fid, float *ptr, int *num);
  void twrite_ (int *fid, long int *ptr, int *num);
  void twritec_ (int *fid, char *ptr, int *num);
  void twrites_ (int *fid, short int *ptr, int *num);
  void tclose_ (int *fid);
#endif

* HP-Fortran Calls *
#ifdef HP
  void topen (int *fid, char *name, int *flag, int *sys);
  void tread (int *fid, long int *ptr, int *num);
  void treadc (int *fid, char *ptr, int *num);
  void treads (int *fid, short int *ptr, int *num);
  void treadf (int *fid, float *ptr, int *num);
  void twrite (int *fid, long int *ptr, int *num);
  void twritec (int *fid, char *ptr, int *num);
  void twrites (int *fid, short int *ptr, int *num);
  void tclose (int *fid);
#endif
*/

#endif
