/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

parsim.c - tree evaluation

=cut

**********/

#include "lvb.h"

/**********

=head1 getplen - CALCULATE TREE LENGTH

=head2 SYNOPSIS

    long getplen(const Branch *const barray, const long root,
    unsigned char **bmat, unsigned char **scratch,
    const long m, const long n, const long *weights);

=head2 DESCRIPTION

Calculates the length of a tree using Fitch parsimony. In LVB, a large
amount of time is spent in this function.

=head2 PARAMETERS

=head3 INPUT

=over4

=item barray

Pointer to the first element of the array containing the tree to be
evaluated.

=item root

Index of the root of the tree to be evaluated. C<barray>[C<root>] is
the root.

=item bmat

Pointer to the first element in an array of pointers, each of which
points to a row in the binary-encoded data matrix that will be used to
calculate the length of the tree. C<bmat>[0..C<n>-1] must point to the
first elements of encoded rows for objects 0..C<n>-1 respectively.
C<bmat> is not const-qualified, since it is necessarily accessed
through non-const-qualified pointers.

=item m

The number of columns (bases) in the encoded data matrix.

=item n

The number of objects in the tree, which is also the number of rows
(objects) in the encoded data matrix.

=item weights

Pointer to first element of an m-element array giving the weights of
the sites. For each site, each change is counted the number of times
given in the corresponding element of this array. For no weighting,
each element of the array should be 1.

=back

=head3 OUTPUT

=over4

=item scratch

C<scratch> must point to an array of C<n>-2 pointers, each of which
must point to an allocated C<m>-element array of unsigned char. Arrays
accessible through C<scratch> are used for temporary storage and have
meaningless contents on return. 

=back

=head2 RETURN

Returns total changes in all characters for the tree accessed through
C<barray>, according to the data matrix accessed through C<bmat>.

=cut

**********/

long getplen(Branch *barray, const long root, const long m,
 const long n, const long *weights)
{
    const long branch_cnt = brcnt(n);	/* branch count */
    long changes = 0;			/* tree length (number of changes) */
    unsigned current_ss;		/* current state set */
    long i;				/* loop counter */
    long k;				/* current character number */
    long left;				/* current left child number */
    unsigned left_ss;			/* left state set */
    long right;				/* current right child number */
    unsigned right_ss;			/* right state set */
    long todo_cnt = 0;			/* count of branches "to do" */
    long internal_cnt = 0;		/* non-leaf non-root branches */
    long done = 0;			/* count of branches "done" */
    long branch;			/* current branch number */
    static long todo_arr[MAX_BRANCHES + 1];	/* list of "dirty" branch nos */
    static long internal_arr[MAX_BRANCHES + 1];	/* list of internal br. nos */

    lvb_assert((n >= MIN_N) && (n <= MAX_N));
    lvb_assert((m >= MIN_M) && (m <= MAX_M));
    lvb_assert((root >= 0) && (root < branch_cnt));

    for (i = 0; i < branch_cnt; i++)
    {
	if (barray[i].sset[0] == 0U)
	{
	    todo_arr[todo_cnt] = i;
	    todo_cnt++;
	}
	if (barray[i].object == UNSET)
	{
	    internal_arr[internal_cnt] = i;
	    internal_cnt++;
	}
    }

    lvb_assert(internal_cnt == (branch_cnt - n));
    lvb_assert(todo_cnt <= (n - 3));	/* max: internal branches */

    /* calculate state sets and changes where not already known */
    while (done < todo_cnt)
    {
	for (i = 0; i < todo_cnt; i++)
	{
	    branch = todo_arr[i];
	    if (barray[branch].sset[0] == 0U)	/* "dirty" */
	    {
		left = barray[branch].left;
		right = barray[branch].right;
		lvb_assert((left >= 0) && (left < branch_cnt));
		lvb_assert((right >= 0) && (right < branch_cnt));
		if ((barray[left].sset[0] != 0U)
		 && (barray[right].sset[0] != 0U))
		{
		    barray[branch].changes = 0;
		    for (k = 0; k < m; k++)
		    {
 			left_ss = barray[left].sset[k];
			right_ss = barray[right].sset[k];
			current_ss = left_ss & right_ss;
			if (current_ss == 0U)
			{
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
    for (i = 0; i < internal_cnt; i++)
    {
	branch = internal_arr[i];
	changes += barray[branch].changes;
    }
    
    /* root: add length for root branch structure, and also for true root which
     * lies outside the LVB tree data structure; all without altering the
     * "root" struct statesets (since these represent actual data for the
     * leaf) */
    left = barray[root].left;
    right = barray[root].right;
    for (k = 0; k < m; k++)
    {
	left_ss = barray[left].sset[k];
	right_ss = barray[right].sset[k];
	current_ss = left_ss & right_ss;
	if (current_ss == 0U)
	{
	    current_ss = left_ss | right_ss;
	    changes += weights[k];
	}
	if ((current_ss & barray[root].sset[k]) == 0U)
	    changes += weights[k];
    }

    return changes;

} /* end getplen() */
