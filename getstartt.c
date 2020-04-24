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

#include "lvb.h"

double AnnealStartingTemperature(Dataptr matrix, const Branch *const inittree, Params rcstruct, long root,
		Lvb_bool log_progress)
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
 * whether to accept them have been adopted from the Anneal()
 * function.
*/ 
{
	/* Variables for the generation of transitions (adopted from Anneal()) */
    double deltah;		/* change in energy (1 - C.I.) */
    long deltalen;		/* change in length with new tree */
    long iter;		/* iteration of mutate/evaluate loop */
    long len;			/* length of current tree */
    long lendash;		/* length of proposed new tree */
    long lenmin;		/* minimum length for any tree */
    double pacc;		/* prob. of accepting new config. */
    double r_lenmin;		/* minimum length for any tree */
    long rootdash;		/* root of new configuration */
    double t = LVB_EPS;		/* current temperature */
    Branch *x;			/* current configuration */
    Branch *xdash;		/* proposed new configuration */

    /* Variables specific to the AnnealStartingTemperatureemperature() procedure*/
    int acc_pos_trans = 0;        /* Number of accepted positve transitions */
    double increment_size = 0.00001; /* Step size by which the temperature is increased */
    int prop_pos_trans = 0;       /* Number of proposed positve transitions */
    double r_acc_to_prop = 0;   /* Ratio of accepted to proposed positve transitions */
    int sample_size = 100;                /* Sample size used to estimate the ratio */
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */
    const double log_wrapper_LVB_EPS = log_wrapper(LVB_EPS);

    /* Create "local" dynamic heap memory and initialise tree 
     * structures like in Anneal() */
    x = AllocBlankTreeArray(matrix, LVB_TRUE);
    xdash = AllocBlankTreeArray(matrix, LVB_TRUE);

    CopyTreeArray(matrix, x, inittree, LVB_TRUE);	/* current configuration */
    AllocCurrentTreeLength(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
    len = CurrentTreeLength(matrix, x, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
    
	lenmin = MinimumTreeLength(matrix);
    r_lenmin = (double) lenmin;
    
    /* Log progress to standard output if chosen*/
    // if (log_progress) printf("\nDetermining the Starting Temperature ...\n");

    while (r_acc_to_prop <= 0.65)
    {

		/* Collect a sample of sample_size permutations at the current temperature 
		* and compute the ratio of proposed vs accepted worse changes*/
		for (iter = 0; iter <= sample_size; iter++)
		{
			/* Create an alternative tree topology (adopted from Anneal()) */

			/* occasionally re-root, to prevent influence from root position */
			if ((iter % REROOT_INTERVAL) == 0) root = RandomTreeRoot(matrix, x, root);

			lvb_assert(t > DBL_EPSILON);

			/* mutation: alternate between the two mutation functions */
			rootdash = root;
			if (iter & 0x01) HeuristicSPR(matrix, xdash, x, root);	/* global change */
			else HeuristicNNI(matrix, xdash, x, root);	/* local change */

			lendash = CurrentTreeLength(matrix, xdash, rcstruct, rootdash, p_todo_arr, p_todo_arr_sum_changes, p_runs);
			lvb_assert (lendash >= 1L);
			deltalen = lendash - len;
			deltah = (r_lenmin / (double) len) - (r_lenmin / (double) lendash);
			
			if (deltah > 1.0)	/* MinimumTreeLength() problem with ambiguous sites */
				deltah = 1.0;

			/* Check whether the change is accepted (Again adopted from Anneal()*/
			if (deltalen <= 0)	/* accept the change */
			{
				/* update current tree and its stats */
				len = lendash;
				SwapTrees(&x, &root, &xdash, &rootdash);
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
						SwapTrees(&x, &root, &xdash, &rootdash);
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
    FreeCurrentTreeLength(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
    free(x);
    free(xdash);
    
    /* Log progress if chosen*/
    if (log_progress)
        printf("  SA Starting Temperature: %-.8f\n", (t - increment_size));
    
    /* Return the temperature last used */
    return (t - increment_size);

} /* end AnnealStartingTemperature() */
