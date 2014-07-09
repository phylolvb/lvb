/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */
    
/* *** getstartt.c - function to determine starting temperature *** */

#include "lvb.h"

double get_initial_t(Dataptr matrix, const Branch *const inittree, long root, long m, long n,
		const long *weights, Lvb_bool log_progress)
/* Determine the starting temperature for the annealing search 
 * by finding the temperature T at which 65% of proposed 
 * positive transitions (changes in the tree structure which increase
 * the tree length) are accepted. Starting at t = LVB_EPS, the
 * algorithm will gradually increase the temperature, estimating the 
 * ratio of accepted to proposed postive transitions at each step
 * using a sample of sample_size transitions. When the ratio reaches the 
 * desired value the search stops and the current temperature is 
 * returned as starting temperature.
 * Note: The procedures for creating mutations and deciding on 
 * whether to accept them have been adopted from the anneal()
 * function.
*/ 
{
	/* Variables for the generation of transitions (adopted from anneal()) */
    double deltah;		/* change in energy (1 - C.I.) */
    long deltalen;		/* change in length with new tree */
    long iter;		/* iteration of mutate/evaluate loop */
    long len;			/* length of current tree */
    long prev_len = UNSET;	/* length of previous tree */
    long lendash;		/* length of proposed new tree */
    long lenmin;		/* minimum length for any tree */
    double pacc;		/* prob. of accepting new config. */
    double r_lenmin;		/* minimum length for any tree */
    long rootdash;		/* root of new configuration */
    double t = LVB_EPS;		/* current temperature */
    Branch *x;			/* current configuration */
    Branch *xdash;		/* proposed new configuration */

    /* Variables specific to the get_initial_temperature() procedure*/
    int acc_pos_trans = 0;        /* Number of accepted positve transitions */
    double increment_size = 0.00001; /* Step size by which the temperature is increased */
    int prop_pos_trans = 0;       /* Number of proposed positve transitions */
    double r_acc_to_prop = 0;   /* Ratio of accepted to proposed positve transitions */
    int sample_size = 100;                /* Sample size used to estimate the ratio */
    
    /* Create "local" dynamic heap memory and initialise tree 
     * structures like in anneal() */
    x = treealloc(matrix);
    xdash = treealloc(matrix);

    treecopy(matrix, x, inittree);	/* current configuration */
    len = getplen(x, root, m, n, weights);
    
    lenmin = getminlen(matrix);
    r_lenmin = (double) lenmin;
    
    /* Log progress to standard output if chosen*/
    if (log_progress) printf("\nDetermining the Starting Temperature ...\n");

    while (r_acc_to_prop <= 0.65)
    {
		/* Collect a sample of sample_size permutations at the current temperature 
		* and compute the ratio of proposed vs accepted worse changes*/
		for (iter = 0; iter <= sample_size; iter++)
		{
			/* Create an alternative tree topology (adopted from anneal()) */
			prev_len = len;

			/* occasionally re-root, to prevent influence from root position */
			if ((iter % REROOT_INTERVAL) == 0) root = arbreroot(matrix, x, root);

			lvb_assert(t > DBL_EPSILON);

			/* mutation: alternate between the two mutation functions */
			rootdash = root;
			if (iter & 0x01) mutate_nni(matrix, xdash, x, root);	/* local change */
			else mutate_spr(matrix, xdash, x, root);	/* global change */

			lendash = getplen(xdash, rootdash, m, n, weights);
			lvb_assert (lendash >= 1L);
			deltalen = lendash - len;
			deltah = (r_lenmin / (double) len) - (r_lenmin / (double) lendash);
			
			if (deltah > 1.0)	/* getminlen() problem with ambiguous sites */
				deltah = 1.0;

			/* Check whether the change is accepted (Again adopted from anneal()*/
			if (deltalen <= 0)	/* accept the change */
			{
				/* update current tree and its stats */
				prev_len = len;
				len = lendash;
				treeswap(&x, &root, &xdash, &rootdash);
			}	
			else {
				prop_pos_trans++; /* Another positive transition has been generated*/

				if (-deltah < t * log_wrapper(LVB_EPS)) {
					pacc = 0.0;
					/* Call uni() even though its not required. It
					* would have been called in LVB 1.0A, so this
					* helps make results identical to results with
					* that version. */
					(void) uni();
				}
				else	/* possibly accept the change */
				{
					pacc = exp_wrapper(-deltah/t);
					if (uni() < pacc)	/* do accept the change */
					{
						prev_len = len;
						len = lendash;
						treeswap(&x, &root, &xdash, &rootdash);
						acc_pos_trans++;  /* The change has been accepted */
					}
				}
			}
		}

		/* Calculate the ratio of accepted to proposed positve transitions 
		 * at the current temperature*/
		r_acc_to_prop = (double) acc_pos_trans / prop_pos_trans;

		/* Increase t and make sure it stays within range*/      
		t += increment_size;
		if (t >= 1 || t <= 0)
			return 1;

		/* Reset variables for next temperature */
		prop_pos_trans = 0;
		acc_pos_trans = 0;
    }
    
    /* free "local" dynamic heap memory */
    free(x);
    free(xdash);
    
    /* Log progress if chosen*/
    if (log_progress)
        printf("Starting Temperature is: %-.8f\n", (t - increment_size));
    
    /* Return the temperature last used */
    return (t - increment_size);

} /* end get_initial_t() */
