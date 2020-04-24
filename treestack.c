/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== sops.c - string manipulation functions ========== */

#include "lvb.h"

static void upsize(Dataptr restrict matrix, Treestack *sp)
/* increase allocation for tree stack *sp */
{
    int i;	/* loop counter */
//    long to_copy = matrix->mssz * sizeof(int);

    sp->size++;
 
    /* allocate for stack itself */
    if (sp->stack == NULL)	/* 1st call, stack does not exist */
    {
        sp->stack = (Treestack_element *) Alloc(sp->size * sizeof(Treestack_element), "initial best tree stack");
        sp->next = 0;
        lvb_assert(sp->size == 1);	/* was incremented above */
    }
    else {
        sp->stack = (Treestack_element *) realloc(sp->stack, sp->size * sizeof(Treestack_element));
        if (sp->stack == NULL){
            CrashVerbosely("out of memory: cannot increase allocation for best tree stack to %ld elements", sp->size);
        }
    }

    /* allocate space within stack */
    for (i = sp->next; i < sp->size; i++){
    	/* sp->stack[i].tree = AllocBlankTreeArray(matrix->n); */
    	sp->stack[i].tree = AllocBlankTreeArray(matrix, LVB_FALSE);
    	sp->stack[i].root = -1;

    	int j;
        /* set memory for sset */
    	sp->stack[sp->next].p_sset = (Objset *) Alloc(matrix->nsets * sizeof(Objset), "object set object arrays");
    	for (j = 0; j < matrix->nsets; j++){
    		sp->stack[sp->next].p_sset[j].set = NULL; // Alloc(to_copy, "object set object arrays");
    		sp->stack[sp->next].p_sset[j].cnt = UNSET;
    	}

    }
 
} /* end upsize() */

static void dopush(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset)
/* push tree in barray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size) upsize(matrix, sp);
    CopyTreeArray(matrix, sp->stack[sp->next].tree, barray, b_with_sset);
    sp->stack[sp->next].root = root;

    /* need to copy the sset_2 to the sp->stack[sp->next].sset  */
    CopySingleStates(matrix, sp->stack[sp->next].p_sset);
    sp->next++;
 
} /* end dopush() */

/**********

=head1 NewTreestack - RETURN A NEW TREE STACK

=head2 SYNOPSIS

    Treestack NewTreestack(void);

=head2 DESCRIPTION

Returns a new tree stack.

=head2 PARAMETERS

None.

=head2 RETURN

Returns a new, empty tree stack.

=cut

**********/ 

Treestack NewTreestack(void)
{
    Treestack s;	/* return value */

    s.size = 0;
    s.next = 0;
    s.stack = NULL;
    return s;

} /* end NewTreestack() */

/**********

=head1 PushTreeOntoTreestack - PUSH TREE ONTO TREE STACK

=head2 SYNOPSIS

    long PushTreeOntoTreestack(Treestack *sp, const Branch *const barray,
    const long root);

=head2 DESCRIPTION

Push copy of a tree onto an existing tree stack. Will not push if its
topology is already present on the stack. The stack will increase its
own memory allocation if necessary.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item barray

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



long PushTreeOntoTreestack(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset)
{
#ifdef LVB_MAPREDUCE 
    long i, new_root = root;			/* loop counter */
	long stackroot;		/* root of current tree */
	static Branch *copy_2 = NULL;	/* possibly re-rooted tree 2 */
	Lvb_bool b_First = LVB_TRUE;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) {
		copy_2 = AllocBlankTreeArray(matrix, b_with_sset);
	}
	CopyTreeArray(matrix, copy_2, barray, b_with_sset);

	/* return before push if not a new topology */
	/* check backwards as similar trees may be discovered together */
	for (i = sp->next - 1; i >= 0; i--) {
		stackroot = sp->stack[i].root;
		if(stackroot != new_root){
			RerootTree(matrix, copy_2, new_root, stackroot, b_with_sset);
			new_root = stackroot;
			b_First = LVB_TRUE;
		}
		if (CompareTrees(matrix, sp->stack[i].p_sset, copy_2, b_First) == 0) return 0;
		b_First = LVB_FALSE;
	}

    #else
#define MIN_THREAD_SEARCH_SSET		5

    long i, slice = 0, slice_tail, new_root = 0;
    static Branch *copy_2 = NULL;			/* possibly re-rooted tree 2 */
    Lvb_bool b_First = LVB_TRUE;
    Lvb_bool b_find_sset = LVB_FALSE;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) copy_2 = AllocBlankTreeArray(matrix, b_with_sset);
    CopyTreeArray(matrix, copy_2, barray, b_with_sset);
    if (root != 0){
    	RerootTree(matrix, copy_2, root, new_root, b_with_sset);
    }

    /* return before push if not a new topology */
    /* check backwards as similar trees may be discovered together */
    if (sp->next == 0){
    	MakeArrayOfStateSets(matrix, copy_2, new_root /* always root zero */);
    }
    else{
    	if (sp->next > MIN_THREAD_SEARCH_SSET) slice = sp->next / matrix->n_threads_getplen;
    	if (sp->next > MIN_THREAD_SEARCH_SSET && slice > 0){
    		MakeArrayOfStateSets(matrix, copy_2, 0 /* always root zero */);
    		slice_tail = (sp->next - (slice * matrix->n_threads_getplen));
    		omp_set_dynamic(0);	  /* disable dinamic threathing */
    		#pragma omp parallel num_threads(matrix->n_threads_getplen) private(i) shared(slice, slice_tail, b_find_sset)
    		{
    			int n_count = 0;
    			int n_end = slice * (omp_get_thread_num() + 1);
    			int n_begin = slice * omp_get_thread_num();
    			if (matrix->n_threads_getplen == (omp_get_thread_num() + 1)) n_end += slice_tail;
    			for (i = n_begin; i < n_end; i++) {
    				if (CompareObjsectStatesWithSecondStates(matrix, sp->stack[i].p_sset) == 0){
						#pragma omp atomic
    					b_find_sset |= LVB_TRUE;
    					break;
    				}
    				if (n_count > 0 && (n_count % 20) == 0 && b_find_sset == LVB_TRUE){
    					break;
    				}
    				n_count += 1;
    			}
    		}
    		if (b_find_sset == LVB_TRUE) return 0;
    	}
    	else{
    		for (i = sp->next - 1; i >= 0; i--) {
    			if (CompareTrees(matrix, sp->stack[i].p_sset, copy_2, b_First) == 0) return 0;
    			b_First = LVB_FALSE;
    		}
    	}
    }

    #endif

    /* topology is new so must be pushed */
    lvb_assert(root < matrix->n);
    dopush(matrix, sp, barray, root, b_with_sset);
    return 1;

} /* end PushTreeOntoTreestack() */

long PushTreeOntoTreestackOnly(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset)
{
    dopush(matrix, sp, barray, root, b_with_sset);
    return 1;

} /* end PushTreeOntoTreestack() */

/**********

=head1 PullTreeOffTreestack - POP TREE OFF TREE STACK

=head2 SYNOPSIS

    long PullTreeOffTreestack(Branch *barray, long *root, Treestack *sp);
    
=head2 DESCRIPTION

Pop a tree off a tree stack.

=head2 PARAMETERS

=head3 OUTPUT

=over 4

=item barray

Pointer to first element of array to contain the popped tree. There
must be sufficient space in this array prior to the call.

=item root

Pointer to scalar that will receive the index of the root in C<barray>.

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

long PullTreeOffTreestack(Dataptr matrix, Branch *barray, long *root, Treestack *sp, Lvb_bool b_with_sset)
{
    long val;	/* return value */

    if (sp->next >= 1){
        sp->next--;
        CopyTreeArray(matrix, barray, sp->stack[sp->next].tree, b_with_sset);
        *root = sp->stack[sp->next].root;

        val = 1;
    }
    else{
        val = 0;
    }

    return val;

} /* end PullTreeOffTreestack() */

int PrintTreestack(Dataptr matrix, Treestack *sp, FILE *const outfp, Lvb_bool onerandom)
{
    const int d_obj1 = 0L;	/* 1st obj. for output trees */
    long root;			/* root of current tree */
    int i;			/* loop counter */
    int lower;			/* lowest index of trees to print */
    int upper;			/* 1 + upper index of trees to print */
    Branch *barray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    barray = AllocBlankTreeArray(matrix, LVB_FALSE);

    if (onerandom == LVB_TRUE)	/* choose one random tree to print */
    {
		lower = RandomNumberGenerator(sp->next - 1);
		upper = lower + 1;
    } else {
		lower = 0;
		upper = sp->next;
    }

    for (i = lower; i < upper; i++) {
        CopyTreeArray(matrix, barray, sp->stack[i].tree, LVB_FALSE);
        if (sp->stack[i].root != d_obj1) RerootTree(matrix, barray, sp->stack[i].root, d_obj1, LVB_FALSE);
        root = d_obj1;
        PrintCurrentTree(matrix, outfp, barray, root);
    }
    if (fflush(outfp) != 0)
    	CrashVerbosely("file write error when writing best trees");

    /* deallocate "local" dynamic heap memory */
    free(barray);
    return upper - lower;	/* number of trees printed */

} /* end PrintTreestack() */

/**********

=head1 FreeTreestack - DEALLOCATE TREE STACK

=head2 SYNOPSIS

    void FreeTreestack(Treestack *sp);

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

void FreeTreestack(Dataptr restrict matrix, Treestack *sp)
/* free all memory in tree stack *sp */
{
	int i;	/* loop counter */

    for (i = 0; i < sp->size; i++){
    	if (sp->stack[i].tree != NULL) free(sp->stack[i].tree);
        sp->stack[i].tree = NULL;
        sp->stack[i].root = -1;

        int j;
        if (sp->stack[i].p_sset != NULL){
        	for (j = 0; j < (int)matrix->nsets; j++){
				if (sp->stack[i].p_sset[j].set != NULL) free(sp->stack[i].p_sset[j].set);
				sp->stack[i].p_sset[j].cnt = UNSET;
			}
        	free(sp->stack[i].p_sset);
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

    void ClearTreestack(Treestack *sp);

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

void ClearTreestack(Treestack *sp)
/* clear stack *sp; note its allocation is not changed */
{
    sp->next = 0;	/* clear stack */

} /* end ClearTreestack() */
