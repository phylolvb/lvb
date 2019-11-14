/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2019 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/* ========== getstartt.c - function to determine starting temperature ========== */

#include "Lvb.h"

#ifndef NP_Implementation
double get_initial_t(Dataptr matrix, const Branch *const inittree, Params rcstruct, long root, int myMPIid, Lvb_bool log_progress)
#else
double get_initial_t(Dataptr matrix, const Branch *const inittree, Params rcstruct, long root, const long *weights, Lvb_bool log_progress)
#endif

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
    long lendash;		/* length of proposed new tree */
    double pacc;		/* prob. of accepting new config. */
    double r_lenmin;		/* minimum length for any tree */
    long rootdash;		/* root of new configuration */
    double t = LVB_EPS;		/* current temperature */
    Branch *x;			/* current configuration */
    Branch *xdash;		/* proposed new configuration */

    /* Variables specific to the get_initial_temperature() procedure*/
    int acc_pos_trans = 0;        /* Number of accepted positve transitions */
	#ifndef NP_Implementation
    double increment_size = INITIAL_INCREMENT; /* Step size by which the temperature is increased */
	#else
	long lenmin; // minimum length of any tree
	double increment_size = 0.00001; // Step size by which temperature is increased
	#endif
    int prop_pos_trans = 0;       /* Number of proposed positve transitions */
    double r_acc_to_prop = 0;   /* Ratio of accepted to proposed positve transitions */
    int sample_size = 100;                /* Sample size used to estimate the ratio */
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */
    const double log_wrapper_LVB_EPS = log_wrapper(LVB_EPS);

    /* Create "local" dynamic heap memory and initialise tree 
     * structures like in anneal() */
    x = treealloc(matrix, LVB_TRUE);
    xdash = treealloc(matrix, LVB_TRUE);

    treecopy(matrix, x, inittree, LVB_TRUE);	/* current configuration */
    alloc_memory_to_getplen(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	
	#ifndef NP_Implementation
    len = getplen(matrix, x, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
    r_lenmin = (double) matrix->min_len_tree;
	#else
	len = getplen(matrix, x, rcstruct, root, weights, p_todo_arr, p_todo_arr_sum_changes, p_runs);
    
    lenmin = getminlen(matrix);
    r_lenmin = (double) lenmin;
	#endif

    /* Log progress to standard output if chosen*/
	if (log_progress) printf(" Starting temperature: ");
    while (r_acc_to_prop <= 0.65)
    {

		/* Collect a sample of sample_size permutations at the current temperature 
		* and compute the ratio of proposed vs accepted worse changes*/
		for (iter = 0; iter <= sample_size; iter++)
		{
			/* Create an alternative tree topology (adopted from anneal()) */

			/* occasionally re-root, to prevent influence from root position */
			if ((iter % REROOT_INTERVAL) == 0) root = arbreroot(matrix, x, root);

			lvb_assert(t > DBL_EPSILON);

			/* mutation: alternate between the two mutation functions */
			rootdash = root;
			if (iter & 0x01) mutate_spr(matrix, xdash, x, root);	/* global change */
			else mutate_nni(matrix, xdash, x, root);	/* local change */

			#ifndef NP_Implementation
			lendash = getplen(matrix, xdash, rcstruct, rootdash, p_todo_arr, p_todo_arr_sum_changes, p_runs);
			#else
			lendash = getplen(matrix, xdash, rcstruct, rootdash, weights, p_todo_arr, p_todo_arr_sum_changes, p_runs);
			#endif

			lvb_assert (lendash >= 1L);
			deltalen = lendash - len;
			deltah = (r_lenmin / (double) len) - (r_lenmin / (double) lendash);
			
			if (deltah > 1.0)	/* getminlen() problem with ambiguous sites */
				deltah = 1.0;

			/* Check whether the change is accepted (Again adopted from anneal()*/
			if (deltalen <= 0)	/* accept the change */
			{
				/* update current tree and its stats */
				len = lendash;
				treeswap(&x, &root, &xdash, &rootdash);
			}	
			else {
				prop_pos_trans++; /* Another positive transition has been generated*/

				if (-deltah < t * log_wrapper_LVB_EPS) {
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
		if (t >= 1 || t <= 0) return 1;

		/* Reset variables for next temperature */
		prop_pos_trans = 0;
		acc_pos_trans = 0;
    }
    /* free "local" dynamic heap memory */
    free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
    free(x);
    free(xdash);
    
    /* Log progress if chosen*/
#ifndef NP_Implementation
#ifdef MAP_REDUCE_SINGLE
    if (log_progress)
    	// printf("%-.8g   Process: %d\n", (t - increment_size), myMPIid);
		printf(" %-.8g\n", (t - increment_size));
#else
    if (log_progress)
        printf("Starting Temperature is:%-.8g   Process:%d   Seed:%d\n", (t - increment_size), myMPIid, rcstruct.seed);
#endif
#else
	if (log_progress)
		printf(" %-.8f\n", (t - increment_size));
#endif
    
    /* Return the temperature last used */
    return (t - increment_size);

} /* end get_initial_t() */