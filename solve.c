/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** solve.c - solving functions ********** */

#include "lvb.h"

static void lenlog(FILE *lengthfp, long iteration, long length, double temperature)
/* write a message to file pointer lengthfp; iteration gives current iteration;
 * crash verbosely on write error */
{
    fprintf(lengthfp, "%-15.8f%-15ld%-15ld\n", temperature, iteration, length); 
    if (ferror(lengthfp))
    {
	crash("file error when logging search progress");
    }

} /* end lenlog() */

long deterministic_hillclimb(Dataptr matrix, Treestack *bstackp, const Branch *const inittree,
		long root, FILE * const lenfp, const long *weights,
		long *current_iter, Lvb_bool log_progress)
/* perform a deterministic hill-climbing optimization on the tree in inittree,
 * using NNI on all internal branches until no changes are accepted; return the
 * length of the best tree found; current_iter should give the iteration number
 * at the start of this call and will be used in any statistics sent to lenfp,
 * and will be updated on return */
{
    long i;				/* loop counter */
    long j;				/* loop counter */
    long todo_cnt = 0;			/* count of internal branches */
    long len;				/* current length */
    long prev_len;			/* previous length */
    long lendash;			/* length of proposed new config */
    long rootdash = root;		/* root of proposed new config */
    long deltalen;			/* change in length */
    Lvb_bool newtree;			/* accepted a new configuration */
    Branch *p_current_tree;				/* current configuration */
    Branch *p_proposed_tree;			/* proposed new configuration */
    static long todo[MAX_BRANCHES];	/* array of internal branch numbers */
    static Lvb_bool leftright[] = {	/* to loop through left and right */
    				LVB_FALSE, LVB_TRUE };

    /* "local" dynamic heap memory */
    p_current_tree = treealloc(matrix);
    p_proposed_tree = treealloc(matrix);

    treecopy(matrix, p_current_tree, inittree);      /* current configuration */
    len = getplen(matrix, p_current_tree, root, weights);
    prev_len = len;

    /* identify internal branches */
    for (i = matrix->n; i < matrix->nbranches; i++) todo[todo_cnt++] = i;
    lvb_assert(todo_cnt == matrix->nbranches - matrix->n);

    do {
		newtree = LVB_FALSE;
		for (i = 0; i < todo_cnt; i++) {
			for (j = 0; j < 2; j++) {
				mutate_deterministic(matrix, p_proposed_tree, p_current_tree, root, todo[i], leftright[j]);
				lendash = getplen(matrix, p_proposed_tree, rootdash, weights);
				lvb_assert (lendash >= 1L);
				deltalen = lendash - len;
				if (deltalen <= 0) {
					if (deltalen < 0)  /* very best so far */
					{
						treestack_clear(bstackp);
						len = lendash;
					}
					if (treestack_push(matrix, bstackp, p_proposed_tree, rootdash) == 1) {
						newtree = LVB_TRUE;
						treeswap(&p_current_tree, &root, &p_proposed_tree, &rootdash);
					}
				}
				if ((log_progress == LVB_TRUE) && ((len != prev_len) || ((*current_iter % STAT_LOG_INTERVAL) == 0))) {
					lenlog(lenfp, *current_iter, len, 0);
				}
				prev_len = len;
				*current_iter += 1;
			}
		}
    } while (newtree == LVB_TRUE);

    /* free "local" dynamic heap memory */
    free(p_current_tree);
    free(p_proposed_tree);

    return len;
}

long anneal(Dataptr matrix, Treestack *bstackp, const Branch *const inittree, long root,
		const double t0, const long maxaccept, const long maxpropose,
		const long maxfail, FILE *const lenfp, const long *weights, long *current_iter,
		const int cooling_schedule, Lvb_bool log_progress)
/* seek parsimonious tree from initial tree in inittree (of root root)
 * with initial temperature t0, and subsequent temperatures obtained by
 * multiplying the current temperature by (t1 / t0) ** n * t0 where n is
 * the ordinal number of this temperature, after at least maxaccept changes
 * have been accepted or maxpropose changes have been proposed, whichever is
 * sooner;
 * return the length of the best tree(s) found after maxfail consecutive
 * temperatures have led to no new accepted solution;
 * lenfp is for output of current tree length and associated details;
 * *current_iter should give the iteration number at the start of this call and
 * will be used in any statistics sent to lenfp, and will be updated on
 * return */
{
    long accepted = 0;		/* changes accespted */
    Lvb_bool dect;		/* should decrease temperature */
    double deltah;		/* change in energy (1 - C.I.) */
    long deltalen;		/* change in length with new tree */
    long failedcnt = 0; 	/* "failed count" for temperatures */
    long iter = 0;		/* iteration of mutate/evaluate loop */
    long len;			/* length of current tree */
    long prev_len = UNSET;	/* length of previous tree */
    long lenbest;		/* bet length found so far */
    long lendash;		/* length of proposed new tree */
    long lenmin;		/* minimum length for any tree */
    double ln_t;		/* ln(current temperature) */
    long t_n = 0;		/* ordinal number of current temperature */
    Lvb_bool newtree;		/* accepted a new configuration */
    double pacc;		/* prob. of accepting new config. */
    Lvb_bool probaccd;		/* have accepted based on Pacc */
    long proposed = 0;		/* trees proposed */
    double r_lenmin;		/* minimum length for any tree */
    long rootdash;		/* root of new configuration */
    double t = t0;		/* current temperature */
    double grad_geom = 0.99;		/* "gradient" of the geometric schedule */
    double grad_linear = 3.64 * LVB_EPS; /* gradient of the linear schedule */
    Branch *p_current_tree;			/* current configuration */
    Branch *p_proposed_tree;		/* proposed new configuration */

    /* "local" dynamic heap memory */
    p_current_tree = treealloc(matrix);
    p_proposed_tree = treealloc(matrix);

    treecopy(matrix, p_current_tree, inittree);	/* current configuration */
    len = getplen(matrix, p_current_tree, root, weights);
    dect = LVB_FALSE;		/* made LVB_TRUE as necessary at end of loop */

    lvb_assert( ((float) t >= (float) LVB_EPS) && (t <= 1.0) && (grad_geom >= LVB_EPS) && (grad_linear >= LVB_EPS));

    lenbest = len;
    treestack_push(matrix, bstackp, inittree, root);	/* init. tree initially best */
    if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {
        fprintf(lenfp, "\nTemperature:   Rearrangement: Length:\n");
    }

    lenmin = getminlen(matrix);
    r_lenmin = (double) lenmin;

    while (1) {
        if ((log_progress == LVB_TRUE) && ((*current_iter % STAT_LOG_INTERVAL) == 0)) {
        	lenlog(lenfp, *current_iter, len, t);
        }
        /*lenlog(lenfp, *current_iter, len, t);*/

		prev_len = len;
		*current_iter += 1;

		/* occasionally re-root, to prevent influence from root position */
		if ((*current_iter % REROOT_INTERVAL) == 0)
			root = arbreroot(matrix, p_current_tree, root);

		lvb_assert(t > DBL_EPSILON);
		newtree = LVB_FALSE;
		probaccd = LVB_FALSE;

		/* mutation: alternate between the two mutation functions */
		rootdash = root;
		if (iter & 0x01) mutate_nni(matrix, p_proposed_tree, p_current_tree, root);	/* local change */
		else mutate_spr(matrix, p_proposed_tree, p_current_tree, root);	/* global change */

		lendash = getplen(matrix, p_proposed_tree, rootdash, weights);
		lvb_assert (lendash >= 1L);
		deltalen = lendash - len;
		deltah = (r_lenmin / (double) len) - (r_lenmin / (double) lendash);
		if (deltah > 1.0) deltah = 1.0; /* getminlen() problem with ambiguous sites */

		if (deltalen <= 0)	/* accept the change */
		{
			if (lendash <= lenbest)	/* store tree if new */
			{
				/*printf("%ld\n", *current_iter);*/
				if (lendash < lenbest) treestack_clear(bstackp);	/* discard old bests */
				if (treestack_push(matrix, bstackp, p_proposed_tree, rootdash) == 1)
					newtree = LVB_TRUE;	/* new tree accepted*/
			}
			/* update current tree and its stats */
			prev_len = len;
			len = lendash;
			treeswap(&p_current_tree, &root, &p_proposed_tree, &rootdash);

			/* very best so far */
			if (lendash < lenbest) lenbest = lendash;
		}
		else	/* poss. accept change for the worse */
		{
			/* Mathematically,
			 *     Pacc = e ** (-1/T * deltaH)
			 *     therefore ln Pacc = -1/T * deltaH
			 *
			 * Computationally, if Pacc is going to be small, we
			 * can assume Pacc is 0 without actually working it
			 * out.
			 * i.e.,
			 *     if ln Pacc < ln eps, let Pacc = 0
			 * substituting,
			 *     if -deltaH / T < ln eps, let Pacc = 0
			 * rearranging,
			 *     if -deltaH < T * ln eps, let Pacc = 0
			 * This lets us work out whether Pacc will be very
			 * close to zero without dividing anything by T. This
			 * should prevent overflow. Since T is no less
			 * than eps and ln eps is going to have greater
			 * magnitude than eps, underflow when calculating
			 * T * ln eps is not possible. */
			if (-deltah < t * log_wrapper(LVB_EPS))
			{
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
					probaccd = LVB_TRUE;
					treeswap(&p_current_tree, &root, &p_proposed_tree, &rootdash);
				}
			}
			if (probaccd == LVB_TRUE){
				prev_len = len;
				len = lendash;
			}
		}
		proposed++;
		if (newtree == LVB_TRUE)
			accepted++;

		/* decide whether to reduce temperature */
		if (accepted >= maxaccept)	/* enough new trees */
		{
			failedcnt = 0;  /* this temperature a 'success' */
			dect = LVB_TRUE;
		}
		else if (proposed >= maxpropose)	/* enough proposals */
		{
			failedcnt++;
			if (failedcnt >= maxfail && t < FROZEN_T)	/* system frozen */
			{
				/* Preliminary experiments yielded that the freezing
				 * criterion used in previous versions of LVB is not
				 * suitable for the new version and regularly results
				 * in premature termination of the search. An easy fix
				 * is to only apply the freezing criterion in the
				 * temperature ranges of previous versions of LVB
				 * (LVB_EPS < t < 10^-4). Future work will look at optimising
				 * maxpropose, maxaccept and maxfail directly. */
				break; /* end of cooling */
			}
			else{	/* system not frozen, so further decrease temp. */
				dect = LVB_TRUE;
			}
		}

		if (dect == LVB_TRUE)
		{
			t_n++;	/* originally n is 0 */

			if (cooling_schedule == 0)  /* Geometric cooling */
			{
				/* Ensure t doesn't go out of bound */
				ln_t = ((double) t_n) * log_wrapper(grad_geom) + log_wrapper(t0);
				if (ln_t < log_wrapper(LVB_EPS)) t = LVB_EPS;
				else t = pow_wrapper(grad_geom, (double) t_n) * t0; /* decrease the temperature */
			}
			else /* Linear cooling */
			{
				t = t0 - grad_linear * t_n;
				/* Make sure t doesn't go out of bounce */
				if (t < DBL_EPSILON || t <= LVB_EPS) t = LVB_EPS;
			}
			proposed = 0;
			accepted = 0;
			dect = LVB_FALSE;
		}
		iter++;
    }

    /* free "local" dynamic heap memory */
    free(p_current_tree);
    free(p_proposed_tree);
    return lenbest;

} /* end anneal() */
