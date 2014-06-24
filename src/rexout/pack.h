#ifndef PACK_H
#define PACK_H

#include "tio3.h"

/* Type 1 rex compression... */
int Stuff_xxx (TIO_type *tp, int val, char f_flag);
int UnStuff_xxx (TIO_type *tp, int *val, char f_flag);

/* Type 2 rex compression... */
int Stuff2_xxx (TIO_type *tp, int val, char f_flag);
int UnStuff2_xxx (TIO_type *tp, int *val, char f_flag);

int UnStuff_LatLon (TIO_type *tp, long int *x, char *sgn, char f_flag);

#endif
