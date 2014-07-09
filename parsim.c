/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** parsim.c - tree length calculation ********** */

#include "lvb.h"

long getplen(Branch *barray, const long root, const long m, const long n, const long *weights)
{
    long branch;			/* current branch number */
    const long branch_cnt = brcnt(n);	/* branch count */
    long changes = 0;			/* tree length (number of changes) */
    unsigned current_ss;		/* current state set */
    long done = 0;			/* count of branches "done" */
    long i;				/* loop counter */
    long k;				/* current character number */
    long left;				/* current left child number */
    unsigned left_ss;			/* left state set */
    long right;				/* current right child number */
    unsigned right_ss;			/* right state set */
    long todo_cnt = 0;			/* count of branches "to do" */
    static long todo_arr[MAX_BRANCHES + 1];	/* list of "dirty" branch nos */

    lvb_assert((n >= MIN_N) && (n <= MAX_N));
    lvb_assert((m >= MIN_M) && (m <= MAX_M));
    lvb_assert((root >= 0) && (root < branch_cnt));

    for (i = n; i < branch_cnt; i++) {
    	if (barray[i].sset[0] == 0U) todo_arr[todo_cnt++] = i;
    }

    /* calculate state sets and changes where not already known */
    while (done < todo_cnt) {
		for (i = 0; i < todo_cnt; i++) {
			branch = todo_arr[i];
			if (barray[branch].sset[0] == 0U)	/* "dirty" */
			{
				left = barray[branch].left;
				right = barray[branch].right;
				if ((barray[left].sset[0] != 0U) && (barray[right].sset[0] != 0U))
				{
					barray[branch].changes = 0;
					for (k = 0; k < m; k++){
						left_ss = barray[left].sset[k];
						right_ss = barray[right].sset[k];
						current_ss = left_ss & right_ss;
						if (current_ss == 0U){
							current_ss = left_ss | right_ss;
							barray[branch].changes += weights[k];
						}
						barray[branch].sset[k] = current_ss;
					}
					done++;
				}
			}
		}
    }

    /* count changes across tree */
    for (i = n; i < branch_cnt; i++) changes += barray[i].changes;

    /* root: add length for root branch structure, and also for true root which
     * lies outside the LVB tree data structure; all without altering the
     * "root" struct statesets (since these represent actual data for the
     * leaf) */
    left = barray[root].left;
    right = barray[root].right;
    for (k = 0; k < m; k++) {
		left_ss = barray[left].sset[k];
		right_ss = barray[right].sset[k];
		current_ss = left_ss & right_ss;
		if (current_ss == 0U){
			current_ss = left_ss | right_ss;
			changes += weights[k];
		}
		if ((current_ss & barray[root].sset[k]) == 0U) changes += weights[k];
    }

    lvb_assert(changes > 0);

    return changes;

} /* end getplen() */
