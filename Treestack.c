#ifdef LVB_NP

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

/* ========== Treestack.c - TREESTACK Operations ========== */

#include "LVB.h"

static void TreestackAllocationIncrease(Dataptr restrict MSA, TREESTACK *sp)
/* increase allocation for tree stack *sp */
{
    int i;	/* loop counter */
//    long to_copy = MSA->mssz * sizeof(int);

    sp->size++;
 
    /* allocate for stack itself */
    if (sp->stack == NULL)	/* 1st call, stack does not exist */
    {
        sp->stack = (TREESTACK_TREE *) alloc(sp->size * sizeof(TREESTACK_TREE), "initial best tree stack");
        sp->next = 0;
        lvb_assert(sp->size == 1);	/* was incremented above */
    }
    else {
        sp->stack = (TREESTACK_TREE *) realloc(sp->stack, sp->size * sizeof(TREESTACK_TREE));
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

long PushCurrentTreeToStack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
/* push tree in BranchArray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size) TreestackAllocationIncrease(MSA, sp);
    CopyCurrentTree(MSA, sp->stack[sp->next].tree, BranchArray, b_with_sitestate);
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

    long CompareTreeToTreestack(TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray,
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

long CompareTreeToTreestack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
{
    long i = 0, new_root = 0;
    static TREESTACK_TREE_BRANCH *current_tree_copy_ptr = NULL;			/* possibly re-rooted tree 2 */
    Lvb_bool b_First = LVB_TRUE;

    #ifndef LVB_MAPREDUCE
    
    #define MIN_THREAD_SEARCH_SSET		5
    long slice = 0, slice_tail = 0;
    Lvb_bool b_find_sitestate = LVB_FALSE;

    #endif

	/* allocate "local" static heap memory - static - do not free! */
	if (current_tree_copy_ptr == NULL) current_tree_copy_ptr = treealloc(MSA, b_with_sitestate);
    CopyCurrentTree(MSA, current_tree_copy_ptr, BranchArray, b_with_sitestate);
    /* CopyCurrentTree(dataset, from, to, sitestates)*/
    if (root != 0){
    	RerootCurrentTree(MSA, current_tree_copy_ptr, root, new_root, b_with_sitestate);
        /* RerootCurrentTree(data, current tree, oldroot, newroot, sitestates) */
    }

    /* return before push if not a new topology */
    /* check backwards as similar trees may be discovered together */
    if (sp->next == 0){
    	makesets(MSA, current_tree_copy_ptr, new_root /* always root zero */);
    }
    #ifndef LVB_MAPREDUCE
    // Multithreading
    else{ 
    	if (sp->next > MIN_THREAD_SEARCH_SSET) slice = sp->next / MSA->n_threads_getplen;
    	if (sp->next > MIN_THREAD_SEARCH_SSET && slice > 0){
    		makesets(MSA, current_tree_copy_ptr, 0 /* always root zero */);
    		slice_tail = (sp->next - (slice * MSA->n_threads_getplen));
    		omp_set_dynamic(0);	  /* disable dinamic threathing */
    		#pragma omp parallel num_threads(MSA->n_threads_getplen) private(i) shared(slice, slice_tail, b_find_sitestate)
    		{
    			int n_count = 0;
    			int n_end = slice * (omp_get_thread_num() + 1);
    			int n_begin = slice * omp_get_thread_num();
    			if (MSA->n_threads_getplen == (omp_get_thread_num() + 1)) n_end += slice_tail;
    			for (i = n_begin; i < n_end; i++) {
    				if (setstcmp_with_sitestate2(MSA, sp->stack[i].p_sitestate) == 0){
						// #pragma omp atomic
    					b_find_sitestate = LVB_TRUE;
    					break;
    				}
    				if (n_count > 0 && (n_count % 20) == 0 && b_find_sitestate == LVB_TRUE){
    					break;
    				}
    				n_count += 1;
    			}
    		}
    		if (b_find_sitestate == LVB_TRUE) return 0;
    	}
    #endif
    	else{
    		for (i = sp->next - 1; i >= 0; i--) {
    			if (TopologyComparison(MSA, sp->stack[i].p_sitestate, current_tree_copy_ptr, b_First) == 0) return 0;
    			b_First = LVB_FALSE;
    		}
    	}
    #ifndef LVB_MAPREDUCE
    }
    #endif

    /* topology is new so must be pushed */
    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);

    FILE *printcurrenttree;
    FILE *printcurrenttreehashcomparison;
    printcurrenttree = fopen("PrintCurrentTree", "w");
    printcurrenttreehashcomparison = fopen("PrintCurrentTreeHashComparison", "a+");
       
    TopologyHashing(MSA, printcurrenttree, BranchArray, root);
    TopologyHashing(MSA, printcurrenttreehashcomparison, BranchArray, root);

    fclose(printcurrenttree);
    fclose(printcurrenttreehashcomparison);

    HashCurrentTree();

    
    return 1;

} /* end CompareTreeToTreestack() */

/**********

=head1 PullTreefromTreestack - POP TREE OFF TREE STACK

=head2 SYNOPSIS

    long PullTreefromTreestack(TREESTACK_TREE_BRANCH *BranchArray, long *root, TREESTACK *sp);
    
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

long PullTreefromTreestack(Dataptr MSA, TREESTACK_TREE_BRANCH *BranchArray, long *root, TREESTACK *sp, Lvb_bool b_with_sitestate)
{
    long val;	/* return value */

    if (sp->next >= 1){
        sp->next--;
        CopyCurrentTree(MSA, BranchArray, sp->stack[sp->next].tree, b_with_sitestate);
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
    TREESTACK_TREE_BRANCH *BranchArray;		/* current unpacked tree */

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
        CopyCurrentTree(MSA, BranchArray, sp->stack[i].tree, LVB_FALSE);
        if (sp->stack[i].root != d_obj1) RerootCurrentTree(MSA, BranchArray, sp->stack[i].root, d_obj1, LVB_FALSE);
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

#elif LVB_PARALLEL_SEARCH

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

/* ========== Treestack.c - TREESTACK Operations ========== */

#include "LVB.h"

#ifdef MPI_Implementation

#include "StoreStates.h"

static void TreestackAllocationIncrease(Dataptr MSA, TREESTACK *sp)
/* increase allocation for tree stack *sp */
{
    long i;	/* loop counter */

    sp->size++;
 
    /* allocate for stack itself */
    if (sp->stack == NULL)	/* 1st call, stack does not exist */
    {
        sp->stack = (TREESTACK_TREE *) alloc(sp->size * sizeof(TREESTACK_TREE), "initial best tree stack");
        sp->next = 0;
        lvb_assert(sp->size == 1);	/* was incremented above */
    }
    else {
        sp->stack = (TREESTACK_TREE *) realloc(sp->stack, sp->size * sizeof(TREESTACK_TREE));
        if (sp->stack == NULL)
            crash("out of memory: cannot increase allocation for\n"
            		"best tree stack to %ld elements", sp->size);
    }

    /* allocate space within stack */
    for (i = sp->next; i < sp->size; i++){
    	sp->stack[i].tree = treealloc(MSA, LVB_FALSE);
    	sp->stack[i].root = -1;
    }
 
} /* end TreestackAllocationIncrease() */

static void PushCurrentTreeToStack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
/* push tree in BranchArray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size) TreestackAllocationIncrease(MSA, sp);
    CopyCurrentTree(MSA, sp->stack[sp->next].tree, BranchArray, b_with_sitestate);
    sp->stack[sp->next].root = root;
    sp->next++;
 
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

/* ********** checkpoint_treestack() - checkpoint the tree stack ********** */

/* The layout of the treestack, when dumped to file, is as follows.
 * o long int: the number of trees on the stack; and
 * o trees, each stored as:
 *     - long int, the root of the tree; followed by
 *     - array of TREESTACK_TREE_BRANCH, the contents of the tree. */

void checkpoint_treestack(FILE *fp, TREESTACK *s, Dataptr MSA, Lvb_bool b_with_sitestate)
{
    long i;					/* loop counter */
    TREESTACK_TREE *stack = s->stack;	/* actual stack */
    TREESTACK_TREE current_element;		/* current element of stack */
    long trees = s->next;			/* number of trees stored */
    long size = s->size;			/* size */

    unsigned long n_bytes_to_write = 2 * sizeof(long) + sizeof(unsigned short);
    /* size of tree stack element */
    if (b_with_sitestate == LVB_TRUE) n_bytes_to_write += trees * (MSA->tree_bytes + sizeof(long));
    else n_bytes_to_write += trees * (MSA->tree_bytes_without_sitestate + sizeof(long));
    unsigned long checksum = 0;
    unsigned short type_block = STATE_BLOCK_TREESTACK;
    fwrite(&n_bytes_to_write, sizeof(n_bytes_to_write), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_write), (unsigned char *) &n_bytes_to_write, checksum);
    fwrite(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);

    /* the number of trees on the stack */
    fwrite(&trees, sizeof(trees), 1, fp); checksum = CalculateBlockCRC32(sizeof(trees), (unsigned char *) &trees, checksum);
    fwrite(&size, sizeof(size), 1, fp); checksum = CalculateBlockCRC32(sizeof(size), (unsigned char *) &size, checksum);

    /* trees */
    for (i = 0; i < trees; i++) {
		current_element = stack[i];
		fwrite(&current_element.root, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(trees), (unsigned char *) &current_element.root, checksum);

		if (b_with_sitestate == LVB_TRUE){
			fwrite(current_element.tree, MSA->tree_bytes, 1, fp);
			checksum = CalculateBlockCRC32(MSA->tree_bytes, (unsigned char *) current_element.tree, checksum);
		}
		else{
			fwrite(current_element.tree, MSA->tree_bytes_without_sitestate, 1, fp);
			checksum = CalculateBlockCRC32(MSA->tree_bytes_without_sitestate, (unsigned char *) current_element.tree, checksum);
		}
    }
    fwrite(&checksum, sizeof(unsigned long), 1, fp);
//   print_information_checkpoint("Save data treestack", n_bytes_to_write, checksum);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(fflush(fp) == 0);
}

/* ********** restore_treestack() - restore the treestack ********** */

void restore_treestack(FILE *fp, TREESTACK *sp, Dataptr MSA, Lvb_bool b_with_sitestate)
{
    long i;			/* loop counter */
    long current_root;		/* current root read in from file */
    TREESTACK_TREE_BRANCH *p_current_tree;	/* current tree read in from file */
    long trees;			/* number of trees */
    long size;			/* number of trees */
	unsigned long n_bytes_to_write = 2 * sizeof(long) + sizeof(unsigned short), n_bytes_to_read = 0;
    unsigned long checksum = 0, checksum_read, n_read_values;
    unsigned short type_block;

    /* 'local' heap memory */
    p_current_tree = treealloc(MSA, b_with_sitestate);

    n_read_values = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_read), (unsigned char *) &n_bytes_to_read, checksum);
    n_read_values = fread(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
    n_read_values = fread(&trees, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(trees), (unsigned char *) &trees, checksum);
    n_read_values = fread(&size, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(size), (unsigned char *) &size, checksum);
    if (b_with_sitestate == LVB_TRUE) n_bytes_to_write += trees * (MSA->tree_bytes + sizeof(long));
    else n_bytes_to_write += trees * (MSA->tree_bytes_without_sitestate + sizeof(long));
    for (i = 0; i < trees; i++) {
    	n_read_values = fread(&current_root, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) &current_root, checksum);
    	if (b_with_sitestate == LVB_TRUE){
    		n_read_values = fread(p_current_tree, MSA->tree_bytes, 1, fp);
    		checksum = CalculateBlockCRC32(MSA->tree_bytes, (unsigned char *) p_current_tree, checksum);
    	}
    	else{
    		n_read_values = fread(p_current_tree, MSA->tree_bytes_without_sitestate, 1, fp);
    		checksum = CalculateBlockCRC32(MSA->tree_bytes_without_sitestate, (unsigned char *) p_current_tree, checksum);
    	}
    	lvb_assert(n_read_values == 1);
		PushCurrentTreeToStack(MSA, sp, p_current_tree, current_root, b_with_sitestate);
    }
    n_read_values = fread(&checksum_read, sizeof(unsigned long), 1, fp);
	lvb_assert(n_bytes_to_read == n_bytes_to_write);
	lvb_assert(checksum_read == checksum);
	lvb_assert(ferror(fp) == 0);
    /* free 'local' heap memory */
    free(p_current_tree);
}


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

TREESTACK * CreateNewTreestack(void)
{
    TREESTACK * p_stack;
    p_stack = (TREESTACK *) alloc(sizeof(TREESTACK), "alloc TREESTACK structure");

    p_stack->size = 0;
    p_stack->next = 0;
    p_stack->stack = NULL;

    return p_stack;
} /* end CreateNewTreestack() */

/**********

=head1 CompareTreeToTreestack - PUSH TREE ONTO TREE STACK

=head2 SYNOPSIS

    long CompareTreeToTreestack(TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray,
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

long CompareTreeToTreestack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
{
	long i, new_root = root;			/* loop counter */
	long stackroot;		/* root of current tree */
	static TREESTACK_TREE_BRANCH *current_tree_copy_ptr = NULL;	/* possibly re-rooted tree 2 */
	Lvb_bool b_First = LVB_TRUE;

	/* allocate "local" static heap memory - static - do not free! */
	if (current_tree_copy_ptr == NULL) {
		current_tree_copy_ptr = treealloc(MSA, b_with_sitestate);
	}
	CopyCurrentTree(MSA, current_tree_copy_ptr, BranchArray, b_with_sitestate);

	/* return before push if not a new topology */
	/* check backwards as similar trees may be discovered together */
	for (i = sp->next - 1; i >= 0; i--) {
		stackroot = sp->stack[i].root;
		if(stackroot != new_root){
			RerootCurrentTree(MSA, current_tree_copy_ptr, new_root, stackroot, b_with_sitestate);
			new_root = stackroot;
			b_First = LVB_TRUE;
		}
		if (TopologyComparison(MSA, sp->stack[i].tree, current_tree_copy_ptr, stackroot, b_First) == 0) return 0;
		b_First = LVB_FALSE;
	}

	/* topology is new so must be pushed */
	lvb_assert(root < MSA->n);

    /* topology is new so must be pushed */
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);
    return 1;

} /* end CompareTreeToTreestack() */

//uint64_t mrStack_push(Dataptr MSA, TREESTACK *sp, TREESTACK_TREE_BRANCH *BranchArray, const long root, MapReduce *mrObj, MISC *misc)
//{
//    	uint64_t nKV = tree_setpush(MSA, BranchArray, root, mrObj, misc);
//    	return nKV;
//}




/**********

=head1 PullTreefromTreestack - POP TREE OFF TREE STACK

=head2 SYNOPSIS

    long PullTreefromTreestack(TREESTACK_TREE_BRANCH *BranchArray, long *root, TREESTACK *sp);
    
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

long PullTreefromTreestack(Dataptr MSA, TREESTACK_TREE_BRANCH *BranchArray, long *root, TREESTACK *sp, Lvb_bool b_with_sitestate)
{
    long val;	/* return value */

    if (sp->next >= 1){
        sp->next--;
        CopyCurrentTree(MSA, BranchArray, sp->stack[sp->next].tree, b_with_sitestate);
        *root = sp->stack[sp->next].root;
        val = 1;
    }
    else{
    	val = 0;
    }
    return val;

} /* end PullTreefromTreestack() */


long PrintTreestack(Dataptr MSA, TREESTACK *sp, FILE *const outfp, Lvb_bool onerandom)
{
    const long d_obj1 = 0L;	/* 1st obj. for output trees */
    long root;			/* root of current tree */
    long i;			/* loop counter */
    long lower;			/* lowest index of trees to print */
    long upper;			/* 1 + upper index of trees to print */
    TREESTACK_TREE_BRANCH *BranchArray;		/* current unpacked tree */

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
        CopyCurrentTree(MSA, BranchArray, sp->stack[i].tree, LVB_FALSE);
        if (sp->stack[i].root != d_obj1) RerootCurrentTree(MSA, BranchArray, sp->stack[i].root, d_obj1, LVB_FALSE);
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

void FreeTreestackMemory(TREESTACK *sp)
/* free all memory in tree stack *sp */
{
    long i;	/* loop counter */

    for (i = 0; i < sp->size; i++){
    	if (sp->stack[i].tree != NULL) free(sp->stack[i].tree);
        sp->stack[i].tree = NULL;
        sp->stack[i].root = -1;
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

#endif


#endif