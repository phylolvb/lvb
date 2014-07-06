/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** mymaths.h - interface for mymaths.c ********** */

#ifndef LVB_MYMATHS_H
#define LVB_MYMATHS_H

double exp_wrapper(double) /*@globals errno@*/ /*@modifies nothing@*/ ;
double log_wrapper(double) /*@globals errno@*/ /*@modifies nothing@*/ ;
double pow_wrapper(double, double) /*@globals errno@*/ /*@modifies nothing@*/ ;

#endif /* LVB_MYMATHS_H */
