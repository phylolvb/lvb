/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
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

/**********

=head1 NAME

treestack.c - tree stack functions

=head1 DESCRIPTION

Provides operations for accessing and maintaining a stack of tree
topologies.

=cut

**********/

#include "lvb.h"

static void upsize(Dataptr matrix, Treestack *sp)
/* increase allocation for tree stack *sp */
{
    long i;	/* loop counter */

    sp->size++;
 
    /* allocate for stack itself */
    if (sp->stack == NULL)	/* 1st call, stack does not exist */
    {
        sp->stack = (Treestack_element *) alloc(sp->size * sizeof(Treestack_element), "initial best tree stack");
        sp->next = 0;
        lvb_assert(sp->size == 1);	/* was incremented above */
    }
    else {
        sp->stack = (Treestack_element *) realloc(sp->stack, sp->size * sizeof(Treestack_element));
        if (sp->stack == NULL)
            crash("out of memory: cannot increase allocation for\n"
            		"best tree stack to %ld elements", sp->size);
    }

    /* MIGUEL */
    /* allocate space within stack */
    for (i = sp->next; i < sp->size; i++){
    	/* sp->stack[i].tree = treealloc(matrix->n); */
    	sp->stack[i].tree = treealloc(matrix, LVB_FALSE); /* MIGUEL */
    	sp->stack[i].root = -1;
    }
 
} /* end upsize() */

static void dopush(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root)
/* push tree in barray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size) upsize(matrix, sp);
    treecopy(matrix, sp->stack[sp->next].tree, barray, LVB_FALSE); /* MIGUEL */
    sp->stack[sp->next].root = root;
    sp->next++;
 
} /* end dopush() */

/**********

=head1 treestack_cnt - RETURN COUNT OF TREES ON STACK

=head2 SYNOPSIS

long treestack_cnt(Treestack s);

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

long treestack_cnt(Treestack s)
{
    return s.next;

} /* end treestack_cnt() */

/**********

=head1 treestack_new - RETURN A NEW TREE STACK

=head2 SYNOPSIS

    Treestack treestack_new(void);

=head2 DESCRIPTION

Returns a new tree stack.

=head2 PARAMETERS

None.

=head2 RETURN

Returns a new, empty tree stack.

=cut

**********/ 

Treestack treestack_new(void)
{
    Treestack s;	/* return value */

    s.size = 0;
    s.next = 0;
    s.stack = NULL;

    return s;

} /* end treestack_new() */

/**********

=head1 treestack_push - PUSH TREE ONTO TREE STACK

=head2 SYNOPSIS

    long treestack_push(Treestack *sp, const Branch *const barray,
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

long treestack_push(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root)
{
    long i;			/* loop counter */
    Branch *stacktree = NULL;	/* current tree on stack */
    long stackroot;		/* root of current tree */

    /* return before push if not a new topology */
    /* check backwards as similar trees may be discovered together */
    for (i = sp->next - 1; i >= 0; i--) {
        stacktree = sp->stack[i].tree;
        stackroot = sp->stack[i].root;
        if (treecmp(matrix, stacktree, stackroot, barray, root) == 0) return 0;
    }

    /* topology is new so must be pushed */
    dopush(matrix, sp, barray, root);
    return 1;

} /* end treestack_push() */


long treestack_push_only(Dataptr matrix, Treestack *sp, const Branch *const barray, const long root)
{

    dopush(matrix, sp, barray, root);
    return 1;
}

//uint64_t mrStack_push(Dataptr matrix, Treestack *sp, Branch *barray, const long root, MapReduce *mrObj, MISC *misc)
//{
//    	uint64_t nKV = tree_setpush(matrix, barray, root, mrObj, misc);
//    	return nKV;
//}




/**********

=head1 treestack_pop - POP TREE OFF TREE STACK

=head2 SYNOPSIS

    long treestack_pop(Branch *barray, long *root, Treestack *sp);
    
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

long treestack_pop(Dataptr matrix, Branch *barray, long *root, Treestack *sp)
{
    long val;	/* return value */

    if (sp->next >= 1){
        sp->next--;
        treecopy(matrix, barray, sp->stack[sp->next].tree, LVB_FALSE);
        *root = sp->stack[sp->next].root;

        val = 1;
    }
    else{
        val = 0;
    }

    return val;

} /* end treestack_pop() */

long treestack_print(Dataptr matrix, Treestack *sp, FILE *const outfp, Lvb_bool onerandom)
{
    const long d_obj1 = 0L;	/* 1st obj. for output trees */
    long root;			/* root of current tree */
    long i;			/* loop counter */
    long lower;			/* lowest index of trees to print */
    long upper;			/* 1 + upper index of trees to print */
    Branch *barray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    barray = treealloc(matrix, LVB_FALSE);

    if (onerandom == LVB_TRUE)	/* choose one random tree to print */
    {
		lower = randpint(sp->next - 1);
		upper = lower + 1;
    } else {
		lower = 0;
		upper = sp->next;
    }

    for (i = lower; i < upper; i++) {
        treecopy(matrix, barray, sp->stack[i].tree, LVB_FALSE);
        if (sp->stack[i].root != d_obj1) lvb_reroot(matrix, barray, sp->stack[i].root, d_obj1, LVB_FALSE);
        root = d_obj1;
        lvb_treeprint(matrix, outfp, barray, root);
    }
    if (fflush(outfp) != 0)
    	crash("file write error when writing best trees");

    /* deallocate "local" dynamic heap memory */
    free(barray);

    return upper - lower;	/* number of trees printed */

} /* end treestack_print() */

/**********

=head1 treestack_dump - DUMP AND CLEAR TREE STACK

=head2 SYNOPSIS

    long treestack_dump(Treestack *sp, FILE *const outfp);

=head2 DESCRIPTION

Pop all trees off a stack and dump them.

=head2 PARAMETERS

=head3 OUTPUT

=over4

=item outfp

Pointer to the file to which we wish to output trees.

=back

=head3 INOUT

=over 4

=item sp

The stack to be emptied and dumped.

=back

=head2 RETURN

Returns the number of trees dumped.

=cut

**********/

long treestack_dump(Dataptr matrix, Treestack *sp, FILE *const outfp)
/* pop all trees on stack *sp and dump them to file outfp;
 * first branch (number 0); return number of trees dumped */
{
    long cnt = 0;		/* tree count */
    long root;			/* number of root branch */
    Branch *barray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    barray = treealloc(matrix, LVB_FALSE);

    while ((treestack_pop(matrix, barray, &root, sp)) != 0){
		treedump(matrix, outfp, barray, LVB_FALSE);
		cnt++;
    }

    if (fflush(outfp) != 0) crash("file write error when dumping best trees");

    /* free "local" dynamic heap memory */
    free(barray);

    return cnt;

} /* end bstdump() */

/**********

=head1 treestack_free - DEALLOCATE TREE STACK

=head2 SYNOPSIS

    void treestack_free(Treestack *sp);

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

void treestack_free(Treestack *sp)
/* free all memory in tree stack *sp */
{
    long i;	/* loop counter */

    for (i = 0; i < sp->size; i++){
    	free(sp->stack[i].tree);
        sp->stack[i].tree = NULL;
        sp->stack[i].root = -1;
    }
    free(sp->stack);
    sp->next = 0;
    sp->size = 0;
    sp->stack = NULL;
 
} /* end bstfree() */

/**********

=head1 treestack_clear - EMPTY TREE STACK

=head2 SYNOPSIS

    void treestack_clear(Treestack *sp);

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

void treestack_clear(Treestack *sp)
/* clear stack *sp; note its allocation is not changed */
{
    sp->next = 0;	/* clear stack */

} /* end treestack_clear() */


/**********

=head1 treestack_transfer - TRANSFER TREES BETWEEN TREE STACKS

=head2 SYNOPSIS

    void treestack_transfer(Treestack *destp, Treestack *sourcep);

=head2 DESCRIPTION

Transfer one tree stack in its entirety to another. Order of the transferred
trees is not preserved. Current contents and current order of the destination
stack are preserved. Duplicate topologies are not transferred. The source
stack is emptied but not deallocated.

=head2 PARAMETERS

=head3 INOUT

=over 4

=item destp

Pointer to the stack to be added to.

=item sourcep

Pointer to the stack to be transferred to C<destp> and cleared.

=back

=head2 RETURN

Number of trees actually transferred (excluding duplicates which
are not transferred).

=cut

**********/

long treestack_transfer(Dataptr matrix, Treestack *destp, Treestack *sourcep)
{
    Branch *barray;		/* current tree, in transit */
    long root;			/* number of root branch */
    long pushed = 0;		/* number of trees transferred */

    /* "local" dynamic heap memory */
    barray = treealloc(matrix, LVB_FALSE);
    while (treestack_pop(matrix, barray, &root, sourcep) == 1) {
        pushed += treestack_push(matrix, destp, barray, root);
    }

    /* free "local" dynamic heap memory */
    free(barray);
    return pushed;
}
