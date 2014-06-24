/*
 * tio3.c --
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
#include "tio3.h"

#include <string.h>
#include <stdlib.h>

typedef struct TIO_List {
  TIO_type *tio;
  struct TIO_List *next;
} TIO_List;

static TIO_List *head = NULL;

/* don't need SYSTEM variable... have _sysINTEL instead. */

/*---------------------------------------------------------------------------
 * tOpen -- Exported
 *
 *   Opens a new file.
 *
 * Variables:
 *   fid         The file ID to associate with this file.
 *   name        The name of the file to open.
 *   access      What kind of file access to open it with.
 *   sys         Same endian 0 or opposite 1 as system.
 *               Alternatively... 2 file made on intel,
 *                                3 file not made on intel.
 *
 * Returns:
 *   NULL on error, pointer to valid structure otherwise..
 *
 * Notes:
 *-------------------------------------------------------------------------*/
TIO_type *tOpen (unsigned short int fid, char *name, TFLAG_ACCESS access,
                 TFLAG_SYSTEM sys) {
  TIO_List *ptr;
  FILE *fp;

  /* Check if input is valid. */
  if (name == NULL)
    return NULL;

  /* Check if we already have fid, or we have already opened name. */
  ptr = head;
  while (ptr != NULL) {
    if ((ptr->tio->fid == fid) || (strcmp (name, ptr->tio->name) == 0))
      return NULL;
    ptr = ptr->next;
  }

  /* Attempt to Open the file. */
  switch (access) {
    case TFLAG_READ:
      if ((fp = fopen (name, "rb")) == NULL)
        return NULL;
      break;
    case TFLAG_WRITE:
      fp = fopen (name, "wb");
      break;
    case TFLAG_APPEND:
      fp = fopen (name, "ab");
      break;
    case TFLAG_RW:
      fp = fopen (name, "r+b");
      break;
    default:
      return NULL;
  }

  /* We managed to open the file, so we will attempt to set asside space for
     the TIO data structure... */
  if ((ptr = (TIO_List *) malloc (sizeof(TIO_List))) == NULL) {
    fclose (fp);
    return NULL;
  }
  if ((ptr->tio = (TIO_type *) malloc (sizeof(TIO_type))) == NULL) {
    free (ptr);
    fclose (fp);
    return NULL;
  }
  if ((ptr->tio->name = (char *) malloc ((strlen(name) +1) * sizeof(char)))
         == NULL) {
    free (ptr->tio);
    free (ptr);
    fclose (fp);
    return NULL;
  }

  ptr->tio->fp = fp;
  strcpy (ptr->tio->name, name);
  ptr->tio->fid = fid;
  ptr->tio->access = access;
#ifdef _sysINTEL
  if (sys == TFLAG_MadeOnIntel) {
    ptr->tio->f_os = 0;
  } else {
    ptr->tio->f_os = 1;
  }
#else
  if (sys == TFLAG_MadeOnIntel) {
    ptr->tio->f_os = 1;
  } else {
    ptr->tio->f_os = 0;
  }
#endif
  ptr->tio->gbuf = 0;
  ptr->tio->gbuf_loc = 0;
  ptr->tio->pbuf = 0;
  ptr->tio->pbuf_loc = 8;
  ptr->tio->f_stream = 0;
  ptr->tio->stream = NULL;
  ptr->tio->stream_len = 0;
  ptr->tio->stream_loc = 0;

  /* link in our new data structure. */
  ptr->next = head;
  head = ptr;
  return ptr->tio;
}

/*---------------------------------------------------------------------------
 * tFind -- Exported
 *
 *   Returns a pointer to an already opened TIO_type given a file ID.
 *
 * Variables:
 *   fid         The file ID in question.
 *
 * Returns:
 *   NULL on error, pointer to valid structure otherwise.
 *
 * Notes:
 *-------------------------------------------------------------------------*/
TIO_type *tFind (unsigned short int fid) {
  TIO_List *ptr = head;
  while (ptr != NULL) {
    if (ptr->tio->fid == fid)
      return ptr->tio;
    ptr = ptr->next;
  }
  return NULL;
}

/*---------------------------------------------------------------------------
 * tClose -- Exported
 *
 *   Close a TIO file, and free the related data structures.
 *
 * Variables:
 *   tio         The Opened file ptr to close.
 *
 * Returns: 1 "true" on error, 0 "false" otherwise.
 *
 * Notes:
 *   May want to flush the bit streams just in case?
 *-------------------------------------------------------------------------*/
BOOL tClose (TIO_type *tio) {
  TIO_List *ptr = head, *ptr2 = NULL;

  /* find tio in the list maintained by head... */
  while (ptr != NULL) {
    if (ptr->tio == tio)
      break;
    ptr2 = ptr;
    ptr = ptr->next;
  }
  if (ptr == NULL) {  /* we didn't find it... */
    return 1;
  }

  /* unlink ptr from head list... */
  if (ptr2 == NULL) {
    head = ptr->next;
  } else {
    ptr2->next = ptr->next;
  }

  /* close / free data pointed to by ptr. */
  fclose (ptr->tio->fp);
  free (ptr->tio->name);
  free (ptr->tio->stream);
  free (ptr->tio);
  free (ptr);
  return 0;
}

/*---------------------------------------------------------------------------
 * tEndian_Switch -- Exported
 *
 *   Toggles the endian'ness of a file.
 *
 * Variables:
 *   tio         The Opened file ptr to adjust.
 *
 * Notes:
 *-------------------------------------------------------------------------*/
void tEndian_Switch (TIO_type *tio) {
  tio->f_os = (BOOL) ((tio->f_os) ? 0 : 1);
}

/*---------------------------------------------------------------------------
 * tGet_Endian -- Exported
 *
 *   Returns a boolean as to whether to treat the file as opposite endian to
 *   the operating system.
 *
 * Variables:
 *   tio         The opened file ptr to inquire from.
 *
 * Returns: "true" if opposite endian, 0 "false" otherwise.
 *
 * Notes:
 *-------------------------------------------------------------------------*/
BOOL tGet_Endian (TIO_type *tio) {
  return (tio->f_os);
}

/*---------------------------------------------------------------------------
 * tFlush -- Exported
 *
 *   Clears the bit buffer for tGetBit
 *
 * Variables:
 *   f_put       True, flush the put byte, false flush the get byte.
 *   tio         The opened file ptr to flush.
 *
 * Returns: number of bits it put onto the stream.
 *
 * Notes:
 *-------------------------------------------------------------------------*/
int tFlush (BOOL f_put, TIO_type *tio) {
  if (f_put) {
    if (tio->pbuf_loc != 8) {
      fputc ((int) tio->pbuf, tio->fp);
      tio->pbuf = 0;
      tio->pbuf_loc = 8;
      return 8;
    } else {
      tio->pbuf = 0;
      tio->pbuf_loc = 8;
      return 0;
    }
  } else {
    tio->gbuf = 0;
    tio->gbuf_loc = 0;
    return 0;
  }
}

/*----------------------------------------------------------------------
 * tSeek -- Exported
 *
 *     To seek to specific locations in the file.  (See fseek)
 *
 * Variables:
 *   tio         The opened file ptr to adjust.
 *   offset      distance from the origin.
 *   origin      SEEK_SET (0) from beginning,
 *               SEEK_CUR (1) from current,
 *               SEEK_END (2) from end.
 *
 * Returns: 0 if successful, non-zero on error.
 *
 * Notes:
 *   If fseek, seeks past end of file, the next call to fgetc is an EOF.
 *----------------------------------------------------------------------*/
int tSeek (TIO_type *tio, long int offset, int origin) {
  return fseek (tio->fp, offset, origin);
}


/*---------------------------------------------------------------------------
 * tGetBit -- Exported
 *
 *   To get bits (instead of bytes) from the file.  Stores the current byte,
 *   and passes the bits that were requested to the user.  Leftover bits, are
 *   stored in a char associated with the file, for future reads.
 *
 * Variables:
 *   tio         The opened file ptr to read from.
 *   dst         The storage place for the data read from file.
 *   dst_len     The size of dst (in bytes)
 *   num_bits    The number of bits to read from the file.
 *
 * Returns: EOF if EOF, 1 if error, 0 otherwise.
 *
 * Notes:
 *   1) Don't check tio->access, because user should be smart enough not
 *     to call this when inappropriate.  If they do, fgetc will return
 *     EOF which denotes an error in this case.
 *   2) Don't ask for more bits than you have room for.
 *   3) dst_len <= 127.
 *   4) dst is usually NOT an array of elements.
 *   5) Call tFlush (tio, 0) to flush the bit buffer.
 *-------------------------------------------------------------------------*/
signed char tGetBit (void *Dst, size_t dst_len, unsigned short int num_bits,
                     TIO_type *tio) {
  static unsigned char BitRay[] = {0, 1, 3, 7, 15, 31, 63, 127, 255};
  register BYTE buf_loc, buf, *ptr;
  BYTE *dst = Dst;           /* do I need this? Yes! */
  size_t num_bytes;
  BYTE dst_loc;
  int c;

  memset (Dst, 0, dst_len);

  if (num_bits == 0) {
    return 0;
  }

  /* Since num_bits is always used with -1, I might as well do --num_bits
     here. */
  num_bytes = ((--num_bits) / 8) +1;  /* 1..8 bits = 1 byte, ... */
  /* Check if dst has enough room for num_bits. */
  if (dst_len < num_bytes) {
    return 1;
  }

  /* num_bits was modified earlier. */
  dst_loc = (BYTE) ((num_bits % 8) +1);
  buf_loc = tio->gbuf_loc;
  buf = tio->gbuf;

#ifdef _sysINTEL           /* most significant byte is last... */
  ptr = dst + (num_bytes -1);
#else                      /* most significant byte is first... */
  ptr = dst + (dst_len - num_bytes);
#endif

  /* Deal with initial "remainder" part (most significant byte) in dst. */
  if (buf_loc >= dst_loc) {
    /* can now deal with entire "remainder". */
#ifdef _sysINTEL
    *(ptr--) |= (BYTE) ((buf & BitRay[buf_loc]) >> (buf_loc - dst_loc));
#else
    *(ptr++) |= (BYTE) ((buf & BitRay[buf_loc]) >> (buf_loc - dst_loc));
#endif
    buf_loc -= dst_loc;
  } else {
    /* need to do 2 calls to deal with entire "remainder". */
    if (buf_loc != 0) {
      *ptr |= (BYTE) ((buf & BitRay[buf_loc]) << (dst_loc - buf_loc));
    }
    /* buf_loc is now 0. so we need more data. */
    /* dst_loc is now dst_loc - buf_loc. */
    if ((c = fgetc(tio->fp)) == EOF) {
      tio->gbuf_loc = buf_loc;
      tio->gbuf = buf;
      return EOF;
    }
    buf = (BYTE) c;                        /* buf_loc should be 8 */
    buf_loc += (BYTE) (8 - dst_loc);       /* 8 - (dst_loc - buf_loc) */
    /* need mask in case right shift with sign extension?
       should be ok since buf is an unsigned char, so it fills with 0s. */
#ifdef _sysINTEL
    *(ptr--) |= (BYTE) (buf >> buf_loc);
#else
    *(ptr++) |= (BYTE) (buf >> buf_loc);
#endif
                     /* buf_loc should now be 8 - (dst_loc - buf_loc) */
  }

  /* note buf_loc < dst_loc from here on... either it is 0 or < 8. */
  /* also dst_loc is always 8 from here out... */
#ifdef _sysINTEL
  while (ptr >= dst) {
#else
  while (ptr < dst+dst_len) {
#endif
    if (buf_loc != 0) {
      *ptr |= (BYTE) ((buf & BitRay[buf_loc]) << (8 - buf_loc));
    }
    /* buf_loc is now 0. so we need more data. */
    if ((c = fgetc(tio->fp)) == EOF) {
      tio->gbuf_loc = buf_loc;
      tio->gbuf = buf;
      return EOF;
    }
    buf = (BYTE) c;
    /* need mask in case right shift with sign extension?
     should be ok since buf is an unsigned char, so it fills with 0s. */
#ifdef _sysINTEL
    *(ptr--) |= (BYTE) (buf >> buf_loc);
#else
    *(ptr++) |= (BYTE) (buf >> buf_loc);
#endif
  }

  tio->gbuf_loc = buf_loc;
  tio->gbuf = buf;
  return 0;
}

/*---------------------------------------------------------------------------
 * tRead -- Exported
 *
 *   To read into dst, which is an array of variables each of which is of
 *   size elem_size.  It does appropriate endian swapping as it reads.
 *
 * Variables:
 *   dst         The storage place for the data read from file.
 *   elem_size   The size of a given element of dst (in bytes)
 *   num_elem    The number of elements to read.
 *   tio         The opened file ptr to read from.
 *
 * Returns: Number of elements read.
 *   Slows it down to check for EOF
 *
 * Notes:
 *  1) Arguments are ordered identically to fread's arguments.
 *  2) Example of a call:
 *       int buff[100];
 *       tRead (&buff, sizeof(int), sizeof(buff)/sizeof(int), tio);
 *  3) It was felt that 1 call to fread with n/2 swaps is faster than
 *     n calls to fgetc.
 *  4) Don't check tio->access, because user should be smart enough not
 *     to call this when inappropriate, and if they do, fread will put
 *     an error they can see by calling ferror.
 *-------------------------------------------------------------------------*/
size_t tRead (void *Dst, size_t elem_size, size_t num_elem, TIO_type *tio) {
/*int tRead (void *Dst, size_t elem_size, size_t num_elem, TIO_type *tio)*/
  size_t ans;
  size_t j;
  BYTE temp, *dst, *ptr, *ptr2;

  ans = fread (Dst, elem_size, num_elem, tio->fp);

  /* File was not made on this system.  Endian fix needed. */
  if (tio->f_os) {
    if ((elem_size != 1) && (ans == num_elem)) {
      dst = Dst;
      for (j=0; j < elem_size * num_elem; j += elem_size) {
        ptr = dst + j;
        ptr2 = dst + (j + elem_size - 1);
        while (ptr2 > ptr) {
          temp = *ptr;
          *(ptr++) = *ptr2;
          *(ptr2--) = temp;
        }
      }
    }
  }
  return ans;
}

/*---------------------------------------------------------------------------
 * tReadStruct -- Exported
 *
 *   To read from file into dst which is a structure of some sort.
 *   It also does appropriate endian swapping.
 *
 * Variables:
 *   dst         The storage place for the data read from file.
 *   dst_len     The length of the dst in bytes.
 *   elem_sizes  A char string describing the struct.  See notes for specs.
 *   tio         The opened file ptr to read from.
 *
 * Returns: 1 "true" if an error occured. 0 otherwise.
 *   see feof and ferror for EOF/errors.
 *
 * Notes:
 *  Specs of elem_sizes:
 *    example "1,4,2,8" => 1 char, 1 long, 1 short, 1 double
 *    Types larger than 9 (ascii 48..57) should be A..Z (ascii 65..90).
 *    (This is to reduce the calls to "atoi")
 *    "4*20," => 20 long ints. So you might have:
 *    "2,G,1,1*20,2".
 *    There is no "0*20" as the "same" endian machine wouldn't see it, so
 *    if I skipped either in dst or the file, it would only work on the
 *    other" endian machine.
 *
 *  1) It was felt that 1 call to fread with n/2 swaps is faster than
 *     n calls to fgetc.
 *  2) Don't check tio->access, because user should be smart enough not
 *     to call this when inappropriate, and if they do, fread will put
 *     an error they can see by calling ferror.
 *  3) If you have pointers, put them at the end, and decrease dst_len, or
 *     put them at the beginning and adjust both Dst and dst_len.
 *  4) For 3 byte ints stored as longs, put them at front or back (same as
 *     pointers), then use tRead to fill them. *
 *-------------------------------------------------------------------------*/
BOOL tReadStruct (void *Dst, size_t dst_len, const char *elem_sizes,
                    TIO_type *tio) {
  const char *ptr1;
  BYTE *dst, temp, *dst1, *dst2;
  size_t cur_size, cur_len;
  unsigned short int array_size;
  unsigned int j;

  if (fread (Dst, dst_len, 1, tio->fp) != 1) {
    return 1;
  }
   /* error check overflow of dst... */
  if (tio->f_os) {
    cur_len = 0;
    ptr1 = elem_sizes;
    dst = Dst;
    while ((*ptr1 != '\0') && (cur_len < dst_len)) {
      if (*ptr1 == '1') {
        if (*(++ptr1) == '*') {
          /* strtoul stops at ',' or '*' or '\0' putting remainder in NULL */
          cur_len += (unsigned short int) strtoul ((++ptr1), (char **)NULL, 10);
          if ((ptr1 = strchr (ptr1, ',')) == NULL)
            break; /* simpler just to break than go back to '\0' space. */
          ptr1++;
        } else {
          cur_len ++;
          if (*ptr1 == ',') { /* ptr1 should be on ',' or '\0' now... */
            ptr1 ++;
          }
        }
      } else if (*ptr1 != '0') {
/* Swapping needed. */
        cur_size = (*ptr1 - '0');
        if (cur_size > 9) {
          cur_size = (*ptr1 - 'A' + 10);
        }
        if (*(++ptr1) == '*') {
          /* strtoul stops at ',' or '*' or '\0' putting remainder in NULL */
          array_size = (unsigned short int) strtoul((++ptr1),
                                                    (char **) NULL, 10);
          for (j=0; j < cur_size * array_size; j += cur_size) {
            dst1 = dst + cur_len + j;
            dst2 = dst1 + cur_size - 1;
            while (dst2 > dst1) {
              temp = *dst1;
              *(dst1++) = *dst2;
              *(dst2--) = temp;
            }
          }
          cur_len += (array_size * cur_size);
          if ((ptr1 = strchr (ptr1, ',')) == NULL)
            break; /* simpler just to break than go back to '\0' space. */
          ptr1++;
        } else {
          dst1 = dst + cur_len;
          dst2 = dst1 + cur_size -1;
          while (dst2 > dst1) {
            temp = *dst1;
            *(dst1++) = *dst2;
            *(dst2--) = temp;
          }
          cur_len += cur_size;
          if (*ptr1 == ',') {  /* ptr1 should be on ',' or '\0' now... */
            ptr1++;
          }
        }
      } else {
        return 1;
      }
    }
    if (cur_len > dst_len) { /* not enough room set asside. */
      return 1;
    }
  }
  return 0;
}

/*---------------------------------------------------------------------------
 * tPutBit -- Exported
 *
 *   To write bits from src out to file.  Leftover bits that aren't on a
 *   full byte boundary, are stored in buffer.  Make sure you remember to
 *   flush the buffer when you have finished your last PutBit call, by
 *   calling tFlush.
 *
 * Variables:
 *   tio         The opened file ptr to read from.
 *   src         The data to put out to file.
 *   src_len     The size of a src (in bytes)
 *   num_bits    The number of bits to write to the file.
 *
 * Returns: 1 if error, 0 otherwise
 *
 * Notes:
 *   1) Does not check tio->access, because user should be smart enough not
 *     to call this when inappropriate.  If they do, fgetc will return
 *     EOF which denotes an error in this case.
 *   2) Don't try to write more bits than you have.
 *   3) src is usually NOT an array of elements.
 *   4) Call tFlush (tio, 1) to flush the bit buffer.
 *-------------------------------------------------------------------------*/
BOOL tPutBit (void *Src, size_t src_len, unsigned short int num_bits,
              TIO_type *tio) {
  register BYTE buf_loc, buf, *ptr;
  BYTE *src = Src;
  size_t num_bytes;
  BYTE src_loc;

  if (num_bits == 0) {
    return 0;
  }
  /* Since num_bits is always used with -1, I might as well do --num_bits
     here. */
  num_bytes = ((--num_bits) / 8) +1;  /* 1..8 bits = 1 byte, ... */
  /* Check if src has enough bits for us to put out. */
  if (src_len < num_bytes) {
    return 1;
  }

  /* num_bits was modified earlier. */
  src_loc = (BYTE) ((num_bits % 8) +1);
  buf_loc = tio->pbuf_loc;
  buf = tio->pbuf;

  /* Get to start of interesting part of src. */
#ifdef _sysINTEL
  ptr = src + (num_bytes -1);
#else
  ptr = src + (src_len - num_bytes);
#endif

  /* Deal with most significant byte in src. */
  if (buf_loc >= src_loc) {
    /* can store entire MSB in buf. */
    /* Mask? ... Safer to do so... Particularly if user has a number where
       she wants us to start saving half way through. */
#ifdef _sysINTEL
    buf |= (BYTE) ((*(ptr--) & ((1 << src_loc) -1)) << (buf_loc - src_loc));
#else
    buf |= (BYTE) ((*(ptr++) & ((1 << src_loc) -1)) << (buf_loc - src_loc));
#endif
    buf_loc -= src_loc;
  } else {
    /* need to do 2 calls to store the MSB. */
    if (buf_loc != 0) {
      buf |= (BYTE) ((*ptr & ((1 << src_loc) -1)) >> (src_loc - buf_loc));
    }
    /* buf_loc is now 0, so we write it out. */
    if (fputc ((int) buf, tio->fp) == EOF) {
      tio->pbuf_loc = buf_loc;
      tio->pbuf = buf;
      return 1;
    }
    buf = (BYTE) 0;
    /* src_loc is now src_loc - buf_loc */
    /* store rest of ptr in buf...
       So left shift by 8 - (src_loc -buf_loc)
       and set buf_loc to 8 - (src_loc - buf_loc) */
    buf_loc += (BYTE) (8 - src_loc);
#ifdef _sysINTEL
    buf |= (BYTE) (*(ptr--) << buf_loc);
#else
    buf |= (BYTE) (*(ptr++) << buf_loc);
#endif
  }
  /* src_loc should always be considered 8 from now on.. */

#ifdef _sysINTEL
  while (ptr >= src) {
#else
  while (ptr < src+src_len) {
#endif
    if (buf_loc == 0) {  /* simple case where buf and src line up..*/
      if (fputc ((int) buf, tio->fp) == EOF) {
        tio->pbuf_loc = buf_loc;
        tio->pbuf = buf;
        return 1;
      }
#ifdef _sysINTEL
      buf = (BYTE) *(ptr--);
#else
      buf = (BYTE) *(ptr++);
#endif
    } else {
      /* no mask since src_loc is considered 8. */
      /* need mask in case right shift with sign extension?
         should be ok since *ptr is an unsigned char so it fills with 0s. */
      buf |= (BYTE) ((*ptr) >> (8 - buf_loc));
      /* buf_loc is now 0, so we write it out. */
      if (fputc ((int) buf, tio->fp) == EOF) {
        tio->pbuf_loc = buf_loc;
        tio->pbuf = buf;
        return 1;
      }
      buf = (BYTE) 0;
      /* src_loc is 8-buf_loc... */
      /* need to left shift by 8 - (8-buf_loc) */
#ifdef _sysINTEL
      buf |= (BYTE) (*(ptr--) << buf_loc);
#else
      buf |= (BYTE) (*(ptr++) << buf_loc);
#endif
    }
  }
  /* We would rather not keep a full bit buffer. */
  if (buf_loc == 0) {
    if (fputc ((int) buf, tio->fp) == EOF) {
      tio->pbuf_loc = buf_loc;
      tio->pbuf = buf;
      return 1;
    }
    buf_loc = 8;
    buf = (BYTE) 0;
  }
  tio->pbuf_loc = buf_loc;
  tio->pbuf = buf;
  return 0;
}

/*---------------------------------------------------------------------------
 * tWrite -- Exported
 *
 *   To write from src to file.  src may be an array of variables each of
 *   which is of size elem_size.  It does appropriate endian swapping as it
 *   writes.
 *
 * Variables:
 *   src         The source for the data to write to the file.
 *   elem_size   The size of a given element of src (in bytes)
 *   num_elem    The number of elements to write.
 *   tio         The opened file ptr to write to.
 *
 * Returns: Number of elements written.
 *   see feof and ferror for EOF/errors.
 *
 * Notes:
 *  1) Arguments are ordered identically to fwrite's arguments.
 *  2) Example of a call:
 *       int buff[100];
 *       tWrite (&buff, sizeof(int), sizeof(buff)/sizeof(int), tio);
 *  3) Couldn't do n/2 swaps and 1 call to fread because we would either
 *     corrupt the user's src, or need to call malloc.  Instead use n
 *     calls to fgetc.
 *  4) Don't check tio->access, because user should be smart enough not
 *     to call this when inappropriate, and if they do, fwrite will put
 *     an error they can see by calling ferror.
 *-------------------------------------------------------------------------*/
size_t tWrite (void *Src, size_t elem_size, size_t num_elem, TIO_type *tio) {
  register BYTE *ptr;
  size_t i, j;
  BYTE *src;

  if (tio->f_os) {
    if (elem_size == 1) {
      return fwrite (Src, elem_size, num_elem, tio->fp);
    } else {
      src = Src;
      /* can't do n/2 swaps because we would corrupt src, which the
         user might need. */
      ptr = src - elem_size -1;
      for (j=0; j < num_elem; j++) {
       /* can do this because ptr loops through all (not n/2) elements */
        ptr += 2 * elem_size;
        for (i=0; i < elem_size; i++) {
          if (fputc ((int) *(ptr--), tio->fp) == EOF)
            return j;
        }
      }
      return num_elem;
    }
  } else {
    return fwrite (Src, elem_size, num_elem, tio->fp);
  }
}

TFLAG_SYSTEM _system = 0;

int tSet (TFLAG_SYSTEM sys) {
  if ((sys!=0) && (sys!=1))
    return -1;
  _system = sys;
  return 0;
}

/*****************************************************************************
 * f2c Calls :: Arthur Taylor TDL
 * Purpose:
 *     To enable FORTRAN code that has been sent through f2c to be able to
 *   call the exported procedures of this module.
 * History: 10/3/97 AAT Commented.
 * Notes:
 ****************************************************************************/
#ifdef F2C
void topen_ (int *fid, char *name, int *flag, int *sys) {
  /* if sys == 2, then FORTRAN is using the endian'ness of the
     last call by 'C' to tSet. */
  if (*sys == 2) {
    tOpen ((unsigned short int) *fid, name, *flag, _system);
  } else {
    tOpen ((unsigned short int) *fid, name, *flag, *sys);
  }
}
void tread_ (int *fid, long int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (long int), *num, tp);
}
void treadc_ (int *fid, char *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (char), *num, tp);
}
void treadf_ (int *fid, float *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (float), *num, tp);
}
void twrite_ (int *fid, long int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tWrite (ptr, sizeof(long int), *num, tp);
}
void twritec_ (int *fid, char *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tWrite (ptr, sizeof(char), *num, tp);
}
void twrites_ (int *fid, short int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tWrite (ptr, sizeof(short int), *num, tp);
}
void tclose_ (int *fid) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tClose (tp);
}
#endif
/*****************************************************************************
 * Hp-FORTRAN Calls :: Arthur Taylor TDL
 * Purpose:
 *     To enable FORTRAN code that has been compiled on the HP or in UNIX to
 *   be able to call the exported procedures of this module.
 * History: 10/3/97 AAT Commented.
 * Notes:
 ****************************************************************************/
#ifdef HP
void topen (int *fid, char *name, int *flag, int *sys) {
  /* if sys == 2, then FORTRAN is using the endian'ness of the
     last call by 'C' to tSet. */
  if (*sys == 2) {
    tOpen ((unsigned short int) *fid, name, *flag, _system);
  } else {
    tOpen ((unsigned short int) *fid, name, *flag, *sys);
  }
}
void tread (int *fid, long int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (long int), *num, tp);
}
void treadc (int *fid, char *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (char), *num, tp);
}
void treads (int *fid, short int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (short int), *num, tp);
}
void treadf (int *fid, float *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tRead (ptr, sizeof (float), *num, tp);
}
void twrite (int *fid, long int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tWrite (ptr, sizeof(long int), *num, tp);
}
void twritec (int *fid, char *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tWrite (ptr, sizeof(char), *num, tp);
}
void twrites (int *fid, short int *ptr, int *num) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tWrite (ptr, sizeof(short int), *num, tp);
}
void tclose (int *fid) {
  TIO_type *tp = tFind ((unsigned short int) *fid);
  tClose (tp);
}
#endif
