/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** myuni.h - header for RNG functions in myuni.c ********** */

#include <float.h>
#include <limits.h>

/* set max. random number seed value suitable for rinit() */
#if 900000001L > INT_MAX
#error LVB WARNING: type int not suitable, try with a 32-bit or larger system
#else
#define MAX_SEED 900000000
#endif  /* if 900000001L > INT_MAX */

/* external uni functions */
double uni(void);
void rinit(int ijkl);
