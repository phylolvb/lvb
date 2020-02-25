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

/* ========== treestack_handling.c - Treestack Manipulation Functions ========== */

#include "lvb.h"

#ifndef NP_Implementation
#include "parallel_checkpointing.h"
#else
#include "optimise_tree.h"
#endif

static void upsize(Dataptr restrict matrix, Treestack *sp)
/* increase allocation for tree stack *sp */
{
    long i, j;	/* loop counter */
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

    /* allocate space within stack */
    for (i = sp->next; i < sp->size; i++){
    	sp->stack[i].tree = treealloc(matrix, LVB_FALSE);
    	sp->stack[i].root = -1;

        /* set memory for sset */
    	sp->stack[sp->next].p_sset = (Objset *) alloc(matrix->nsets * sizeof(Objset), "object set object arrays");
    	for (j = 0; j < matrix->nsets; j++){
    		sp->stack[sp->next].p_sset[j].set = NULL; // alloc(to_copy, "object set object arrays");
    		sp->stack[sp->next].p_sset[j].cnt = UNSET;
        }
    }
 
} /* end upsize() */

static void dopush(Dataptr restrict matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset)
/* push tree in barray (of root root) on to stack *sp */
{
    lvb_assert(sp->next <= sp->size);
    if (sp->next == sp->size) upsize(matrix, sp);
    treecopy(matrix, sp->stack[sp->next].tree, barray, b_with_sset);
    sp->stack[sp->next].root = root;
    /* need to copy the sset_2 to the sp->stack[sp->next].sset  */
    copy_sset(matrix, sp->stack[sp->next].p_sset);
    sp->next++;
 
} /* end dopush() */

#ifndef NP_Implementation
/* ********** checkpoint_treestack() - checkpoint the tree stack ********** */

/* The layout of the treestack, when dumped to file, is as follows.
  o long int: the number of trees on the stack; and
  o trees, each stored as:
    - long int, the root of the tree; followed by
    - array of Branch, the contents of the tree. */

void checkpoint_treestack(FILE *fp, Treestack *s, Dataptr restrict matrix, Lvb_bool b_with_sset)
{
    long i;					/* loop counter */
    Treestack_element *stack = s->stack;	/* actual stack */
    Treestack_element current_element;		/* current element of stack */
    long trees = s->next;			/* number of trees stored */
    long size = s->size;			/* size */

    unsigned long n_bytes_to_write = 2 * sizeof(long) + sizeof(unsigned short);
    /* size of tree stack element */
    if (b_with_sset == LVB_TRUE) n_bytes_to_write += trees * (matrix->tree_bytes + sizeof(long));
    else n_bytes_to_write += trees * (matrix->tree_bytes_without_sset + sizeof(long));
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

		if (b_with_sset == LVB_TRUE){
			fwrite(current_element.tree, matrix->tree_bytes, 1, fp);
			checksum = CalculateBlockCRC32(matrix->tree_bytes, (unsigned char *) current_element.tree, checksum);
		}
		else{
			fwrite(current_element.tree, matrix->tree_bytes_without_sset, 1, fp);
			checksum = CalculateBlockCRC32(matrix->tree_bytes_without_sset, (unsigned char *) current_element.tree, checksum);
		}
    }
    fwrite(&checksum, sizeof(unsigned long), 1, fp);
//   print_information_checkpoint("Save data treestack", n_bytes_to_write, checksum);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(fflush(fp) == 0);
}

/* ********** restore_treestack() - restore the treestack ********** */

void restore_treestack(FILE *fp, Treestack *sp, Dataptr restrict matrix, Lvb_bool b_with_sset)
{
    long i;			/* loop counter */
    long current_root;		/* current root read in from file */
    Branch *p_current_tree;	/* current tree read in from file */
    long trees;			/* number of trees */
    long size;			/* number of trees */
	unsigned long n_bytes_to_write = 2 * sizeof(long) + sizeof(unsigned short), n_bytes_to_read = 0;
    unsigned long checksum = 0, checksum_read, n_read_values;
    unsigned short type_block;

    /* 'local' heap memory */
    p_current_tree = treealloc(matrix, b_with_sset);

    n_read_values = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_read), (unsigned char *) &n_bytes_to_read, checksum);
    n_read_values = fread(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
    n_read_values = fread(&trees, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(trees), (unsigned char *) &trees, checksum);
    n_read_values = fread(&size, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(size), (unsigned char *) &size, checksum);
    if (b_with_sset == LVB_TRUE) n_bytes_to_write += trees * (matrix->tree_bytes + sizeof(long));
    else n_bytes_to_write += trees * (matrix->tree_bytes_without_sset + sizeof(long));
    for (i = 0; i < trees; i++) {
    	n_read_values = fread(&current_root, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) &current_root, checksum);
    	if (b_with_sset == LVB_TRUE){
    		n_read_values = fread(p_current_tree, matrix->tree_bytes, 1, fp);
    		checksum = CalculateBlockCRC32(matrix->tree_bytes, (unsigned char *) p_current_tree, checksum);
    	}
    	else{
    		n_read_values = fread(p_current_tree, matrix->tree_bytes_without_sset, 1, fp);
    		checksum = CalculateBlockCRC32(matrix->tree_bytes_without_sset, (unsigned char *) p_current_tree, checksum);
    	}
    	lvb_assert(n_read_values == 1);
		treestack_push_only(matrix, sp, p_current_tree, current_root, b_with_sset);
    }
    n_read_values = fread(&checksum_read, sizeof(unsigned long), 1, fp);
	lvb_assert(n_bytes_to_read == n_bytes_to_write);
	lvb_assert(checksum_read == checksum);
	lvb_assert(ferror(fp) == 0);
    /* free 'local' heap memory */
    free(p_current_tree);
}
#endif

Treestack * treestack_new(void)
{
    Treestack * p_stack;
    p_stack = (Treestack *) alloc(sizeof(Treestack), "alloc Treestack structure");

    p_stack->size = 0; // s.size = 0;
    p_stack->next = 0; // s.next = 0;
    p_stack->stack = NULL; // s.stack = NULL;

    return p_stack; //return s;
} /* end treestack_new() */

//needs reducing
long treestack_push(Dataptr restrict matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset)
{
    #define MIN_THREAD_SEARCH_SSET		5

#ifndef NP_Implementation
	long i, new_root = root;			/* loop counter */
	long stackroot;		/* root of current tree */
	static Branch *copy_2 = NULL;	/* possibly re-rooted tree 2 */
	Lvb_bool b_First = LVB_TRUE;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) {
		copy_2 = treealloc(matrix, b_with_sset);
	}
	treecopy(matrix, copy_2, barray, b_with_sset);

	/* return before push if not a new topology */
	/* check backwards as similar trees may be discovered together */
	for (i = sp->next - 1; i >= 0; i--) {
		stackroot = sp->stack[i].root;
		if(stackroot != new_root){
			lvb_reroot(matrix, copy_2, new_root, stackroot, b_with_sset);
			new_root = stackroot;
			b_First = LVB_TRUE;
		}
		if (treecmp(matrix, sp->stack[i].tree, copy_2, stackroot, b_First) == 0) return 0;
		b_First = LVB_FALSE;
	}
#else
    int i, slice = 0, slice_tail, new_root = 0;
    static Branch *copy_2 = NULL;			/* possibly re-rooted tree 2 */
    Lvb_bool b_First = LVB_TRUE;
    Lvb_bool b_find_sset = LVB_FALSE;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) copy_2 = treealloc(matrix, b_with_sset);
    treecopy(matrix, copy_2, barray, b_with_sset);
    if (root != 0){
    	lvb_reroot(matrix, copy_2, root, new_root, b_with_sset);
    }

    /* return before push if not a new topology */
    /* check backwards as similar trees may be discovered together */
    if (sp->next == 0){
    	makesets(matrix, copy_2, new_root /* always root zero */);
    }
    else{
    	if (sp->next > MIN_THREAD_SEARCH_SSET) slice = sp->next / matrix->n_threads_getplen;
    	if (sp->next > MIN_THREAD_SEARCH_SSET && slice > 0){
    		makesets(matrix, copy_2, 0 /* always root zero */);
    		slice_tail = (sp->next - (slice * matrix->n_threads_getplen));
    		omp_set_dynamic(0);	  /* disable dinamic threathing */
    		#pragma omp parallel num_threads(matrix->n_threads_getplen) private(i) shared(slice, slice_tail, b_find_sset)
    		{
    			int n_count = 0;
    			int n_end = slice * (omp_get_thread_num() + 1);
    			int n_begin = slice * omp_get_thread_num();
    			if (matrix->n_threads_getplen == (omp_get_thread_num() + 1)) n_end += slice_tail;
    			for (i = n_begin; i < n_end; i++) {
    				if (setstcmp_with_sset2(matrix, sp->stack[i].p_sset) == 0){
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
    			if (treecmp(matrix, sp->stack[i].p_sset, copy_2, b_First) == 0) return 0;
    			b_First = LVB_FALSE;
    		}
    	}
    }
#endif
    /* topology is new so must be pushed */
    lvb_assert(root < matrix->n);
    dopush(matrix, sp, barray, root, b_with_sset);
    return 1;

} /* end treestack_push() */

#ifndef NP_Implementation
long treestack_push_only(Dataptr restrict matrix, Treestack *sp, const Branch *const barray, const long root, Lvb_bool b_with_sset)
{
    dopush(matrix, sp, barray, root, b_with_sset);
    return 1;
}
#endif

long treestack_pop(Dataptr restrict matrix, Branch *barray, long *root, Treestack *sp, Lvb_bool b_with_sset)
{
    long val;	/* return value */

    if (sp->next >= 1){
        sp->next--;
        treecopy(matrix, barray, sp->stack[sp->next].tree, b_with_sset);
        *root = sp->stack[sp->next].root;
        val = 1;
    }
    else{
    	val = 0;
    }
    return val;

} /* end treestack_pop() */

long treestack_print(Dataptr restrict matrix, Treestack *sp, FILE *const outfp, Lvb_bool onerandom)
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

long treestack_dump(Dataptr restrict matrix, Treestack *sp, FILE *const outfp)
/* pop all trees on stack *sp and dump them to file outfp;
 * first branch (number 0); return number of trees dumped */
{
    long cnt = 0;		/* tree count */
    long root;			/* number of root branch */
    Branch *barray;		/* current unpacked tree */

    /* "local" dynamic heap memory */
    barray = treealloc(matrix, LVB_FALSE);

    while ((treestack_pop(matrix, barray, &root, sp, LVB_FALSE)) != 0){
		treedump(matrix, outfp, barray, LVB_FALSE);
		cnt++;
    }

    if (fflush(outfp) != 0) crash("file write error when dumping best trees");

    /* free "local" dynamic heap memory */
    free(barray);

    return cnt;

} /* end treestack_dump() */

//needs reducing
void treestack_free(Dataptr restrict matrix, Treestack *sp)
/* free all memory in tree stack *sp */
{
    long i, j; // loop counter

    for (i = 0; i < sp->size; i++){
    	if (sp->stack[i].tree != NULL) free(sp->stack[i].tree);
        sp->stack[i].tree = NULL;
        sp->stack[i].root = -1;

        if (sp->stack[i].p_sset != NULL){
        	for (j = 0; j < (long)matrix->nsets; j++){
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

void treestack_clear(Treestack *sp)
/* clear stack *sp; note its allocation is not changed */
{
    sp->next = 0;	/* clear stack */

} /* end treestack_clear() */

long treestack_transfer(Dataptr restrict matrix, Treestack *destp, Treestack *sourcep, Lvb_bool b_with_sset)
{
    Branch *barray;		/* current tree, in transit */
    long root;			/* number of root branch */
    long pushed = 0;		/* number of trees transferred */

    /* "local" dynamic heap memory */
    barray = treealloc(matrix, b_with_sset);
    while (treestack_pop(matrix, barray, &root, sourcep, b_with_sset) == 1) {
        pushed += treestack_push(matrix, destp, barray, root, b_with_sset);
    }

    /* free "local" dynamic heap memory */
    free(barray);
    return pushed;
}
