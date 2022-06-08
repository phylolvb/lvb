/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and
Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2020 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
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

/* ========== Treestack.c - treestack Operations ========== */

#include "Treestack.h"

static void TreestackAllocationIncrease(Dataptr restrict MSA, TREESTACK *sp)
/* increase allocation for tree stack *sp */
{
    int i;	/* loop counter */
//    long to_copy = MSA->mssz * sizeof(int);

    sp->size++;

    /* allocate for stack itself */
    if (sp->stack == NULL)	/* 1st call, stack does not exist */
    {
        sp->stack = (TREESTACK_TREES *) alloc(sp->size * sizeof(TREESTACK_TREES), "initial best tree stack");
        sp->next = 0;
        lvb_assert(sp->size == 1);	/* was incremented above */
    }
    else {
        sp->stack = (TREESTACK_TREES *) realloc(sp->stack, sp->size * sizeof(TREESTACK_TREES));
        if (sp->stack == NULL){
            crash("out of memory: cannot increase allocation for best tree stack to %ld elements", sp->size);
        }
    }

    /* allocate space within stack */
    for (i = sp->next; i < sp->size; i++){
    	/* sp->stack[i].tree = treealloc(MSA->n); */
    	sp->stack[i].tree = treealloc(MSA, LVB_FALSE);
    	sp->stack[i].root = -1;

    	int j;
        /* set memory for sitestate */
    	sp->stack[sp->next].p_sitestate = (Objset *) alloc(MSA->nsets * sizeof(Objset), "object set object arrays");
    	for (j = 0; j < MSA->nsets; j++){
    		sp->stack[sp->next].p_sitestate[j].set = NULL; // alloc(to_copy, "object set object arrays");
    		sp->stack[sp->next].p_sitestate[j].cnt = UNSET;
    	}

    }

} /* end TreestackAllocationIncrease() */

long PushCurrentTreeToStack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray, const long root, Lvb_bool b_with_sitestate)
/* push tree in BranchArray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size) TreestackAllocationIncrease(MSA, sp);
    treecopy(MSA, sp->stack[sp->next].tree, BranchArray, b_with_sitestate);
    sp->stack[sp->next].root = root;

    /* need to copy the sitestate_2 to the sp->stack[sp->next].sitestate  */
    copy_sitestate(MSA, sp->stack[sp->next].p_sitestate);
    sp->next++;

    return 1;

} /* end PushCurrentTreeToStack() */

/**********

=head1 CountTreestack - RETURN COUNT OF TREES ON STACK

=head2 SYNOPSIS

long CountTreestack(TREESTACK s);

=head2 DESCRIPTION

Return the number of trees currently stored on the stack C<s>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item s

The tree stack whose size we wish to know.

=back

=head2 RETURN

Returns number of trees on stack C<s>.

=cut

**********/

long CountTreestack(TREESTACK s)
{
    return s.next;

} /* end CountTreestack() */

/**********

=head1 CreateNewTreestack - RETURN A NEW TREE STACK

=head2 SYNOPSIS

    TREESTACK CreateNewTreestack(void);

=head2 DESCRIPTION

Returns a new tree stack.

=head2 PARAMETERS

None.

=head2 RETURN

Returns a new, empty tree stack.

=cut

**********/

TREESTACK CreateNewTreestack(void)
{
    TREESTACK s;	/* return value */

    s.size = 0;
    s.next = 0;
    s.stack = NULL;
    return s;

} /* end CreateNewTreestack() */

/**********

=head1 CompareTreeToTreestack - PUSH TREE ONTO TREE STACK

=head2 SYNOPSIS

    long CompareTreeToTreestack(TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray,
    const long root);

=head2 DESCRIPTION

Push copy of a tree onto an existing tree stack. Will not push if its
topology is already present on the stack. The stack will increase its
own memory allocation if necessary.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item BranchArray

Pointer to first element of array containing tree to be pushed.

=item root

Root branch number of tree to be pushed.

=back

=head3 INOUT

=over 4

=item sp

Pointer to tree stack that will receive a new copy of the tree.

=back

=head2 RETURN

Returns 1 if the tree was pushed, or 0 if not.

=cut

**********/

long CompareTreeToTreestack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray, const long root, Lvb_bool b_with_sitestate)
{
#define MIN_THREAD_SEARCH_SSET 5

    long i, slice, slice_tail, new_root = 0;
    static TREESTACK_TREE_NODES *copy_2 = NULL;			/* possibly re-rooted tree 2 */
    Lvb_bool b_First = LVB_TRUE;
    Lvb_bool b_find_sset = LVB_FALSE;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
    treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0){
    	lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

    /* return before push if not a new topology */
    /* check backwards as similar trees may be discovered together */

    if (sp->next == 0){
     	makesets(MSA, copy_2, new_root /* always root zero */);
    } else{
    	if (sp->next > MIN_THREAD_SEARCH_SSET) slice = sp->next / MSA->n_threads_getplen;
    	if (sp->next > MIN_THREAD_SEARCH_SSET && slice > 0){
    		makesets(MSA, copy_2, 0 /* always root zero */);
    		slice_tail = (sp->next - (slice * MSA->n_threads_getplen));
    		omp_set_dynamic(0);	  /* disable dinamic threathing */
    		#pragma omp parallel num_threads(MSA->n_threads_getplen) private(i) shared(slice, slice_tail, b_find_sset)
    		{
    			int n_count = 0;
    			int n_end = slice * (omp_get_thread_num() + 1);
    			int n_begin = slice * omp_get_thread_num();
    			if (MSA->n_threads_getplen == (omp_get_thread_num() + 1)) n_end += slice_tail;
    			for (i = n_begin; i < n_end; i++) {
    				if (setstcmp_with_sitestate2(MSA, sp->stack[i].p_sitestate) == 0){
    					b_find_sset = LVB_TRUE;
    					break;
    				}
    				if (n_count > 0 && (n_count % 20) == 0 && b_find_sset == LVB_TRUE){
    					break;
    				}
    				n_count += 1;
    			}
    		}
    		if (b_find_sset == LVB_TRUE) return 0;
    	} else{
            for (i = sp->next - 1; i >= 0; i--) {
    		    if (TopologyComparison(MSA, sp->stack[i].p_sitestate, copy_2, b_First) == 0) return 0; // If trees are the same, return 0
                b_First = LVB_FALSE;
              }
          }
    }

    /* topology is new so must be pushed */
    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);

    return 1;

} /* end CompareTreeToTreestack() */

/**********

=head1 PullTreefromTreestack - POP TREE OFF TREE STACK

=head2 SYNOPSIS

    long PullTreefromTreestack(TREESTACK_TREE_NODES *BranchArray, long *root, TREESTACK *sp);

=head2 DESCRIPTION

Pop a tree off a tree stack.

=head2 PARAMETERS

=head3 OUTPUT

=over 4

=item BranchArray

Pointer to first element of array to contain the popped tree. There
must be sufficient space in this array prior to the call.

=item root

Pointer to scalar that will receive the index of the root in C<BranchArray>.

=back

=head3 INOUT

=over 4

=item sp

Pointer to tree stack from which we will pop the tree.

=back

=head2 RETURN

Returns 1 if a tree was popped, or 0 if the stack was empty.

=cut

**********/

long PullTreefromTreestack(Dataptr MSA, TREESTACK_TREE_NODES *BranchArray, long *root, TREESTACK *sp, Lvb_bool b_with_sitestate)
{
    long val;	/* return value */

    if (sp->next >= 1){
        sp->next--;
        treecopy(MSA, BranchArray, sp->stack[sp->next].tree, b_with_sitestate);
        *root = sp->stack[sp->next].root;

        val = 1;
    }
    else{
        val = 0;
    }

    return val;

} /* end PullTreefromTreestack() */

int PrintTreestack(Dataptr MSA, TREESTACK *sp, FILE *const outfp, Lvb_bool onerandom)
{
    const int d_obj1 = 0L;	/* 1st obj. for output trees */
    long root;			/* root of current tree */
    int i;			/* loop counter */
    int lower;			/* lowest index of trees to print */
    int upper;			/* 1 + upper index of trees to print */
    TREESTACK_TREE_NODES *BranchArray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    BranchArray = treealloc(MSA, LVB_FALSE);

    if (onerandom == LVB_TRUE)	/* choose one random tree to print */
    {
		lower = randpint(sp->next - 1);
		upper = lower + 1;
    } else {
		lower = 0;
		upper = sp->next;
    }

    for (i = lower; i < upper; i++) {
        treecopy(MSA, BranchArray, sp->stack[i].tree, LVB_FALSE);
        if (sp->stack[i].root != d_obj1) lvb_reroot(MSA, BranchArray, sp->stack[i].root, d_obj1, LVB_FALSE);
        root = d_obj1;
        lvb_treeprint(MSA, outfp, BranchArray, root);
    }
    if (fflush(outfp) != 0)
    	crash("file write error when writing best trees");

    /* deallocate "local" dynamic heap memory */
    free(BranchArray);
    return upper - lower;	/* number of trees printed */

} /* end PrintTreestack() */

/**********

=head1 FreeTreestackMemory - DEALLOCATE TREE STACK

=head2 SYNOPSIS

    void FreeTreestackMemory(TREESTACK *sp);

=head2 DESCRIPTION

Clear a tree stack and deallocate dynamically allocated heap
memory associated with it.

=head2 PARAMETERS

=head3 INOUT

=over 4

=item sp

The stack to be emptied and deallocated.

=back

=head2 RETURN

None.

=cut

**********/

void FreeTreestackMemory(Dataptr restrict MSA, TREESTACK *sp)
/* free all memory in tree stack *sp */
{
	int i;	/* loop counter */

    for (i = 0; i < sp->size; i++){
    	if (sp->stack[i].tree != NULL) free(sp->stack[i].tree);
        sp->stack[i].tree = NULL;
        sp->stack[i].root = -1;

        int j;
        if (sp->stack[i].p_sitestate != NULL){
        	for (j = 0; j < (int)MSA->nsets; j++){
				if (sp->stack[i].p_sitestate[j].set != NULL) free(sp->stack[i].p_sitestate[j].set);
				sp->stack[i].p_sitestate[j].cnt = UNSET;
			}
        	free(sp->stack[i].p_sitestate);
        }
    }
    free(sp->stack);
    sp->next = 0;
    sp->size = 0;
    sp->stack = NULL;

} /* end bstfree() */

/**********

=head1 ClearTreestack - EMPTY TREE STACK

=head2 SYNOPSIS

    void ClearTreestack(TREESTACK *sp);

=head2 DESCRIPTION

Empty a tree stack but do not deallocate memory associated with it.
This memory will be available for re-use when trees are pushed onto the
stack again.

=head2 PARAMETERS

=head3 INOUT

=over 4

=item sp

The stack to be emptied.

=back

=head2 RETURN

None.

=cut

**********/

void ClearTreestack(TREESTACK *sp)
/* clear stack *sp; note its allocation is not changed */
{
    sp->next = 0;	/* clear stack */

} /* end ClearTreestack() */
