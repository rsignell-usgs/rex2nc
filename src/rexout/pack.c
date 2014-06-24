#include "pack.h"

/* Undefine if debugging. */
#define NO_INPUT_CHECK

/*****************************************************************************
 * <Stuff_xxx> :: Arthur Taylor TDL
 * Purpose:
 *     To Stuff basin data, using my current compression scheme for basin
 *   data.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   FID         (I) File Id of opened binary file. (using TIO.c)
 *   val         (I) The slosh value to store (10*height or 999)
 *   f_flag      (I) 0 => normal 1-> flush the stream and bit_buffer
 *   endian      (I) The endian of the machine 0 for DOS, 1 for HP
 *
 * Returns: -1 on error, 0..n number of bits saved this time.
 * History:
 *   4/20/98 AAT Started
 * Notes:
 *   Uses a pcx like encoding scheme.
 *   >  huffman codes are only good if we have a delta?
 *   >  A better method might use delta's for non-land flags.
 *   >  Also use huff codes for run lengths.
 *   Simple scheme....
 *   Always print 8 bit fixed run [1...256]
 *   followed by fixed 9 bit literal [-15.0..36.0]
 *   literal 999 = 36.1
 *   >  for speed sake may want entire array passed in.
 ****************************************************************************/
int Stuff_xxx (TIO_type *tp, int val, char f_flag) {
  static unsigned int prev = 0, run = 0;
  int count = 0;

  if (f_flag == 1) {
/*  output prev/run */
    if (run != 0) {
      run = run - 1;
      tPutBit (&run, sizeof (run), 8, tp);
      tPutBit (&prev, sizeof (prev), 9, tp);
      count += 17;
      run = 0;
    }
    /* Clear the put byte. Only time count += TPUTBIT */
    count += tFlush (1, tp);
    return count;
  }
  if (((val > 360) && (val != 999)) || (val < -150))
    return -1;
  if (val == 999)
    val = 361;
  val += 150;

  if (val == prev) {
    run++;
    if (run == 256) {
      run = run - 1;
      tPutBit (&run, sizeof (run), 8, tp);
      tPutBit (&prev, sizeof (prev), 9, tp);
      count += 17;
      run = 0;
    }
  } else {
/*    output prev/run (if there was one)*/
    if (run != 0) {
      run = run - 1;
      tPutBit (&run, sizeof (run), 8, tp);
      tPutBit (&prev, sizeof (prev), 9, tp);
      count += 17;
    }
    prev = val;
    run = 1;
  }
  return count;
}

/*****************************************************************************
 * <Stuff2_xxx> :: Arthur Taylor TDL
 * Purpose:
 *     To Stuff basin data, using my current compression scheme for basin
 *   data.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   FID         (I) File Id of opened binary file. (using TIO.c)
 *   val         (I) The slosh value to store (10*height or 999)
 *   f_flag      (I) 0 => normal 1-> flush the stream and bit_buffer
 *   endian      (I) The endian of the machine 0 for DOS, 1 for HP
 *
 * Returns: -1 on error, 0..n number of bits saved this time.
 * History:
 *   4/20/98 AAT Started
 * Notes:
 *   Uses a pcx like encoding scheme.
 *   >  huffman codes are only good if we have a delta?
 *   >  A better method might use delta's for non-land flags.
 *   >  Also use huff codes for run lengths.
 *   Simple scheme....
 *   Always print 8 bit fixed run [1...256]
 *   followed by fixed 10 bit literal [-32.0..70.0]
 *   literal 999 = 70.1
 *   >  for speed sake may want entire array passed in.
 ****************************************************************************/
int Stuff2_xxx (TIO_type *tp, int val, char f_flag) {
  static unsigned int prev = 0, run = 0;
  int count = 0;

  if (f_flag == 1) {
/*  output prev/run */
    if (run != 0) {
      run = run - 1;
      tPutBit (&run, sizeof (run), 8, tp);
      tPutBit (&prev, sizeof (prev), 10, tp);
      count += 18;
      run = 0;
    }
    /* Clear the put byte. Only time count += TPUTBIT */
    count += tFlush (1, tp);
    return count;
  }
  if (((val > 700) && (val != 999)) || (val < -320))
    return -1;
  if (val == 999)
    val = 701;
  val += 320;

  if (val == prev) {
    run++;
    if (run == 256) {
      run = run - 1;
      tPutBit (&run, sizeof (run), 8, tp);
      tPutBit (&prev, sizeof (prev), 10, tp);
      count += 18;
      run = 0;
    }
  } else {
/*    output prev/run (if there was one)*/
    if (run != 0) {
      run = run - 1;
      tPutBit (&run, sizeof (run), 8, tp);
      tPutBit (&prev, sizeof (prev), 10, tp);
      count += 18;
    }
    prev = val;
    run = 1;
  }
  return count;
}

/*****************************************************************************
 * <UnStuff_xxx> :: Arthur Taylor TDL
 * Purpose:
 *     To Un-Stuff basin data, using my current compression scheme for basin
 *   data.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   FID         (I) File Id of opened binary file. (using TIO.c)
 *   val         (0) The slosh value to store (10*height or 999)
 *   f_flag      (I) 0 => normal 1-> flush the stream and bit_buffer
 *   endian      (I) The endian of the machine 0 for DOS, 1 for HP
 *
 * Returns: -1 on error, 0 otherwise
 * History:
 *   4/20/98 AAT Started
 * Notes:
 *   Uses a pcx like encoding scheme.
 *   >  huffman codes are only good if we have a delta?
 *   >  A better method might use delta's for non-land flags.
 *   >  Also use huff codes for run lengths.
 *   Simple scheme....
 *   Always print 8 bit fixed run [1...256]
 *   followed by fixed 9 bit literal [-15.0..36.0]
 *   literal 999 = 36.1
 *   > for speed sake, may want entire array passed in?
 ****************************************************************************/
int UnStuff_xxx (TIO_type *tp, int *val, char f_flag) {
  static int prev = 0; /*prev = [-150..360] or 999 */
  static unsigned int run = 0;
  unsigned int temp;

  if (f_flag == 1) {
    /* Clear the get byte. */
    tFlush (0, tp);
    run = 0;
    *val = 0;
    return 0;
  }
  if (run != 0) {
    *val = prev;
    run--;
    return 0;
  } else {
    tGetBit (&run, sizeof (run), 8, tp);
    run += 1;
    tGetBit (&temp, sizeof (temp), 9, tp);
    prev = temp -150;
    if (prev == 361)
      prev = 999;
    *val = prev;
    run--;
    return 0;
  }
}

/*****************************************************************************
 * <UnStuff2_xxx> :: Arthur Taylor TDL
 * Purpose:
 *     To Un-Stuff basin data, using my current compression scheme for basin
 *   data.
 *
 * Variables: (I=input) (O=output) (G=global)
 *   FID         (I) File Id of opened binary file. (using TIO.c)
 *   val         (0) The slosh value to store (10*height or 999)
 *   f_flag      (I) 0 => normal 1-> flush the stream and bit_buffer
 *
 * Returns: -1 on error, 0 otherwise
 * History:
 *   4/20/98 AAT Started
 * Notes:
 *   Uses a pcx like encoding scheme.
 *   >  huffman codes are only good if we have a delta?
 *   >  A better method might use delta's for non-land flags.
 *   >  Also use huff codes for run lengths.
 *   Simple scheme....
 *   Always print 8 bit fixed run [1...256]
 *   followed by fixed 10 bit literal [-32.0..70.0]
 *   literal 999 = 70.1
 *   > for speed sake, may want entire array passed in?
 ****************************************************************************/
int UnStuff2_xxx (TIO_type *tp, int *val, char f_flag) {
  static int prev = 0; /*prev = [-320..700] or 999 */
  static unsigned int run = 0;
  unsigned int temp;

  if (f_flag == 1) {
    /* Clear the get byte. */
    tFlush (0, tp);
    run = 0;
    *val = 0;
    return 0;
  }
  if (run != 0) {
    *val = prev;
    run--;
    return 0;
  } else {
    tGetBit (&run, sizeof (run), 8, tp);
    run += 1;
    tGetBit (&temp, sizeof (temp), 10, tp);
    prev = temp -320;
    if (prev == 701)
      prev = 999;
    *val = prev;
    run--;
    return 0;
  }
}

/* Stuff_latlon routine uses the huffman code choosen for lat/lon files,
   and methods such as a [10000] for change in sign, [0][5bits+extra] for
   literal [1][4bits+extra][4bits+extra] for length, distance pair.
*/
/* Stuff_slosh would use a huffman code choosen for SLOSH data, and might
   combine the lit alphabet with the length alphabet.
*/

/* f_flag == 0 normal f_flag == 1 flush */
/* f_flag == 2 if end of stream. */
/* len == -1 if dist is a literal. */

unsigned char ltln_ray[32] = {0,0,1,2,2,2,3,3,
                         4,4,5,5,6,6,7,7,
                         8,8,9,9,10,10,11,11,
                         11,13,14,14,16,17,19,19};
unsigned long int ltln_val[32] = {0,1,2,4,8,12,16,24,
                      32,48,64,96,128,192,256,384,
                      512,768,1024,1536,2048,3072,4096,6144,
                      8192,10240,18432,34816L,51200L,116736L,247808L,772096L};

/* Idea... look at up to 6 bit at a time (Size of a end of stream char)...
   buffer extra.  */
int UnStuff_LatLon (TIO_type *tp, long int *x, char *sgn, char f_flag) {
  unsigned char lit, index, lit2, first, second;
  unsigned long int enc_x;
  static unsigned char last = 0, len = 0;

  if (f_flag == 0) {
/* need 4,5, or 6 bits depending on last, and len. */
    if (len == 2) {
      tGetBit (&lit2, sizeof (lit2), 4, tp);
      lit = (unsigned char) (((last & 3) << 4) + lit2);
    } else if (len == 1) {
      tGetBit (&lit2, sizeof (lit2), 5, tp);
      lit = (unsigned char) (((last & 1) << 5) + lit2);
    } else {
      tGetBit (&lit, sizeof (lit), 6, tp);
    }
    first = (unsigned char) (lit&48);
    if (first == 32) {
      *x = (*sgn) * ((lit & 12) >> 2);
      last = lit;
      len = 2;
      return 0;
    } else if (first == 16) {
      *x = (*sgn) * ((lit & 28) >> 2);
      last = lit;
      len = 2;
      return 0;
    } else if (first == 0) {
      if (lit >= 8) {
        *x = (*sgn) * ((lit >> 1) + 4);
        last = lit;
        len = 1;
        return 0;
      } else if (lit >= 4) {
        if (lit == 7) {
        /* End of Stream.. return -2 */
          len = 0;
          return -2;
        } else {
          *x = (*sgn) * (lit + 8);
          len = 0;
          return 0;
        }
      } else {
        len = 0;
        tGetBit (&lit2, sizeof (lit2), 3, tp);
        index = (unsigned char) (((lit & 3) << 3) + lit2);
        enc_x = 0;
        if (ltln_ray[index] != 0) {
          tGetBit (&enc_x, sizeof (enc_x), ltln_ray[index], tp);
        }
/* NOTE: dist - 15 before.  We add 15 at this end.*/
        *x = (*sgn) * (ltln_val[index] + enc_x + 15);
        return 0;
      }
    } else {
      *sgn = (char) (-1 *(*sgn));
      first = (unsigned char) (lit&12);
      second = (unsigned char) (lit&3);
      len = 0;
      if (first == 8) {
        *x = (*sgn) * second;
        return 0;
      } else if (first == 4) {
        *x = (*sgn) * (second + 4);
        return 0;
/* First has to == 0 (ie else if (first == 0) ) */
      } else if (first == 0) {
        if (second >= 2) {
          tGetBit (&lit, sizeof (lit), 1, tp);
          second = (unsigned char) ((second << 1) + lit);
          *x = (*sgn) * (second + 4);
          return 0;
        } else if (second == 1) {
          tGetBit (&lit, sizeof (lit), 2, tp);
          if (lit == 3) {
            /* End of Stream.. return -2 */
            return -2;
          } else {
            *x = (*sgn) * (lit + 12);
            return 0;
          }
        } else {
          tGetBit (&index, sizeof (index), 5, tp);
          enc_x = 0;
          if (ltln_ray[index] != 0) {
            tGetBit (&enc_x, sizeof (enc_x), ltln_ray[index], tp);
          }
/* NOTE: dist - 15 before.  We add 15 at this end.*/
          *x = (*sgn) * (ltln_val[index] + enc_x + 15);
          return 0;
        }
#ifdef NO_INPUT_CHECK
      }
#else
      } else {
        FILE *fp;
        fp = fopen ("error.txt", "at");
        fprintf (fp, "shouldn't have 2 sign changed in a row.\n");
        fclose (fp);
        return -1;
      }
#endif
    }
    return 0;
  } else if (f_flag == 1) {
    tFlush (0, tp);
    last = 0;
    len = 0;
    return 0;
  } else {
    return 0;
  }
}
