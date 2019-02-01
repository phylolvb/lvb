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

/* ********** trops.c - tree operations ********** */

#include "lvb.h"

#define CLADESEP ","	/* clade separator for trees */


/* static void cr_nbo(const Branch *const barray, const long obj); */
/* static void cr_tbo(const Branch *const barray, const long obj); */

static void cr_bpnc(const Branch *const barray, const long branch);
static void cr_chaf(const Branch *const barray, const long destination, const long newchild);
static void cr_uxe(FILE *const stream, const char *const msg);
static void fillsets(Dataptr, Objset *const sstruct, const Branch *const tree, const int root);
static void getobjs(Dataptr, const Branch *const barray, const int root, int *const objarr, int *const cnt);
static int getsister(const Branch *const barray, const long branch);
static int osetcmp(const void *oset1, const void *oset2);
static int *randleaf(Dataptr restrict, Branch *const barray,const Lvb_bool *const leafmask);
static void realgetobjs(Dataptr restrict, const Branch *const barray, const int root, int *const objarr, int *const cnt);
static Lvb_bool *randtopology(Dataptr restrict, Branch *const barray, const long nobjs);
static long setstcmp(Dataptr restrict, Objset *const oset_1, Objset *const oset_2);
static void sort(Dataptr matrix, Objset *const oset_2, const long nels);
static void ssarralloc(Dataptr restrict matrix, Objset *nobjset_2);
static void tree_make_canonical(Dataptr restrict, Branch *const barray, int *objnos);
static void ur_print(Dataptr restrict, FILE *const stream, const Branch *const barray, const long root);

/* object sets for tree 2 in comparison */
static Objset sset_2[MAX_N - 3] = { { NULL, 0 } };

void nodeclear(Branch *const barray, const long brnch)
/* Initialize all scalars in branch brnch to UNSET or zero as appropriate,
 * and mark it "dirty" */
{
    barray[brnch].left = UNSET;
    barray[brnch].right = UNSET;
    barray[brnch].parent = UNSET;
    barray[brnch].changes = UNSET;
    barray[brnch].sset[0] = 0U;		/* "dirty" */

} /* end nodeclear() */

long tree_bytes(Dataptr restrict matrix)
/* return bytes required for contiguous allocation of a tree for the data
 * accessible by matrix, if branches and their statesets are allocated
 * as one contiguous array */
{
	/* return in bytes */
    return (matrix->nbranches * sizeof(Branch)) + (matrix->nbranches * matrix->bytes);
} /* end tree_bytes() */

long tree_bytes_without_sset(Dataptr restrict matrix)
/* return bytes required for contiguous allocation of a tree for the data
 * accessible by matrix, */
{
	/* return in bytes */
    return (matrix->nbranches * sizeof(Branch));
} /* end tree_bytes() */



void treeclear(Dataptr restrict matrix, Branch *const barray)
/* clear all branches in array barray, on the assumption that its size fits
 * the data matrix; mark all branches dirty */
{
    long i;					/* loop counter */
    for (i = 0; i < matrix->nbranches; i++) nodeclear(barray, i);

} /* end treeclear() */


static void make_dirty_below(Dataptr restrict matrix, Branch *tree, long dirty_node)
/* mark nodes "dirty" from branch dirty_node, which must not be the root,
 * down to (but not including) the root branch of the tree tree; the true
 * root lies outside the LVB tree data structure so cannot be marked
 * dirty, but will always be dirty after any rearrangement */
{

    lvb_assert(dirty_node >= matrix->n);	/* not leaf/root */
    lvb_assert(tree[dirty_node].parent != UNSET);
    do {
		tree[dirty_node].sset[0] = 0U;	/* " make dirty" */
		dirty_node = tree[dirty_node].parent;
    } while (tree[dirty_node].parent != UNSET);

} /* end make_dirty_below() */

/* static void make_dirty_tree(Dataptr restrict matrix, Branch *tree) */
/* mark all branches in tree tree as dirty: internal, external and root */
// {
//    long i, j;					/* loop counter */
//    for (i = 0; i < matrix->nbranches; i++){
//    	for (j = 0; j < matrix->nwords; j++){ /* overkill beyond j=0, but harmless */
//    		tree[i].sset[j] = 0U;
//    	}
//    }
// } /* end make_dirty_tree() */

void mutate_deterministic(Dataptr restrict matrix, Branch *const desttree,
    const Branch *const sourcetree, long root, long p, Lvb_bool left)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a change in topology, the change being caused by nearest
 * neighbour interchange (NNI) rearrangement at branch p, involving the
 * right node if right is LVB_TRUE, otherwise the left node; N.B. code
 * is largely copied from mutate_nni() */
{
    Branch *tree;
    long u, v, a, b, c;

    lvb_assert(p != root);
    lvb_assert(p >= matrix->n);

    /* for ease of reading, make alias of desttree, tree */
    tree = desttree;
    treecopy(matrix, tree, sourcetree, LVB_TRUE);

    u = p;
    v = tree[u].parent;
    a = tree[u].left;
    b = tree[u].right;
 
    if (tree[v].left == u) c = tree[v].right;
    else c = tree[v].left;

    if (left != LVB_TRUE)
    {
		if (tree[v].left == u) tree[v].right = b;
		else tree[v].left = b;
		tree[u].left = a;
		tree[u].right = c;
		tree[a].parent = tree[c].parent = u;
		tree[b].parent = v;
    }
    else
    {
		if (tree[v].left == u) tree[v].right = a;
		else tree[v].left = a;
		tree[u].left = b;
		tree[u].right = c;
		tree[b].parent = tree[c].parent = u;
		tree[a].parent = v;
    }

    make_dirty_below(matrix, tree, u);

} /* end mutate_nni() */

void mutate_nni(Dataptr restrict matrix, Branch *const desttree, const Branch *const sourcetree, int root)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a random change in topology, the change being caused by nearest
 * neighbour interchange (NNI) rearrangement; N.B. the code is mostly
 * the same in mutate_deterministic() */
{
    Branch *tree;
    int u, v, a, b, c;

    /* for ease of reading, make alias of desttree, tree */
    tree = desttree;
    treecopy(matrix, tree, sourcetree, LVB_TRUE);

    /* get a random internal branch */
    u = randpint(matrix->nbranches - matrix->n - 1) + matrix->n;
    v = tree[u].parent;
    a = tree[u].left;
    b = tree[u].right;
 
    if (tree[v].left == u) c = tree[v].right;
    else c = tree[v].left;

    if (uni() < 0.5)  {
		if (tree[v].left == u) tree[v].right = b;
		else tree[v].left = b;
		tree[u].left = a;
		tree[u].right = c;
		tree[a].parent = tree[c].parent = u;
		tree[b].parent = v;
    }
    else {
		if (tree[v].left == u) tree[v].right = a;
		else tree[v].left = a;
		tree[u].left = b;
		tree[u].right = c;
		tree[b].parent = tree[c].parent = u;
		tree[a].parent = v;
    }

    make_dirty_below(matrix, tree, u);

} /* end mutate_nni() */

static Lvb_bool is_descendant(Branch *tree, int root, int ancestor, int candidate)
/* return LVB_TRUE if candidate is among the branches in the clade descending
from ancestor, LVB_FALSE otherwise */
{
    Lvb_bool val = LVB_FALSE;		/* return value */
    int par = tree[candidate].parent;	/* current parent */
    int newpar;			/* next parent */

    while (par != UNSET){
		if (par == ancestor){
			val = LVB_TRUE;
			par = UNSET;	/* exit the loop */
		}
		else{
			newpar = tree[par].parent;
			par = newpar;
		}
    }
    return val;

} /* end is_descendant() */

void mutate_spr(Dataptr restrict matrix, Branch *const desttree, const Branch *const sourcetree, int root)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a random change in topology, the change being caused by subtree
 * pruning and regrafting (SPR) rearrangement */
{
    int src;				/* branch to move */
    int dest;				/* destination of branch to move */
    int dest_parent;			/* parent of destination branch */
    int src_parent;			/* parent of branch to move */
    int excess_br;			/* branch temporarily excised */
    int orig_child = UNSET;		/* original child of destination */
    int parents_par;			/* parent of parent of br. to move */
    int src_sister;			/* sister of branch to move */
    Branch *tree;			/* destination tree */

    /* for ease of reading, make alias of desttree, tree */
    tree = desttree;
    treecopy(matrix, tree, sourcetree, LVB_TRUE);

    /* get random branch but not root and not root's immediate descendant */
    do {
    	src = randpint(matrix->nbranches - 1);
    } while ((src == root) || (src == tree[root].left) || (src == tree[root].right));

    src_parent = tree[src].parent;
    lvb_assert(src_parent != UNSET);
    src_sister = getsister(tree, src);
    lvb_assert(src_sister != UNSET);

    /* get destination that is not source or its parent, sister or descendant
     * or the root */
    do {
    	dest = randpint(matrix->nbranches - 1);
    } while ((dest == src) || (dest == src_parent) || (dest == src_sister)
       || (dest == root) || is_descendant(tree, root, src, dest));

    /* excise source branch, leaving a damaged data structure */
    if (tree[src_parent].left == src) {
    	tree[src_parent].left = UNSET;
    }
    else if (tree[src_parent].right == src) {
    	tree[src_parent].right = UNSET;
    }
    else {
    	cr_bpnc(tree, src);
    }
    tree[src].parent = UNSET;

    /* fix data structure by "freeing" the excess branch */
    parents_par = tree[src_parent].parent;
    lvb_assert(parents_par != UNSET);
    if (tree[parents_par].left == src_parent) {
    	tree[parents_par].left = src_sister;
    }
    else {
    	tree[parents_par].right = src_sister;
    }
    tree[src_sister].parent = parents_par;

    excess_br = src_parent;	/* for ease of human understanding */
    nodeclear(tree, excess_br);

    /* make space at destination, re-using the excess branch */
    dest_parent = tree[dest].parent;
    if (tree[dest_parent].left == dest) {
    	orig_child = tree[dest_parent].left;
    	tree[dest_parent].left = excess_br;
    }
    else if (tree[dest_parent].right == dest) {
    	orig_child = tree[dest_parent].right;
    	tree[dest_parent].right = excess_br;
    }
    else {
    	crash("destination %ld is not a child of it's parent %ld\n", dest, dest_parent);
    }
    tree[excess_br].parent = dest_parent;
    tree[excess_br].left = dest;
    lvb_assert(orig_child != UNSET);
    tree[orig_child].parent = excess_br;

    /* add source branch to this new location */
    tree[excess_br].right = src;
    tree[src].parent = excess_br;

    /* ensure recalculation of lengths where necessary */
    make_dirty_below(matrix, tree, excess_br);
    if (parents_par != root){
    	make_dirty_below(matrix, tree, parents_par);
    }
} /* end mutate_spr() */

void mutate_tbr(Dataptr restrict matrix, Branch *const desttree, const Branch *const sourcetree, int root)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a random change in topology, the change being caused by subtree
 * pruning and regrafting (SPR) rearrangement */
{
    int src;				/* branch to move */
    int dest;				/* destination of branch to move */
    int dest_parent;			/* parent of destination branch */
    int src_parent;			/* parent of branch to move */
    int excess_br;			/* branch temporarily excised */
    int orig_child = UNSET;		/* original child of destination */
    int parents_par;			/* parent of parent of br. to move */
    int src_sister;			/* sister of branch to move */
    Branch *tree;			/* destination tree */

		int oldroot;	
		int current;							/* current branch */
		int parnt;								/* parent of current branch */
		int sister = UNSET;					/* sister of current branch */
		int oldsister = UNSET;
		int tempsister = UNSET;
		int previous = UNSET;							/* previous branch */
		int newroot;	
		static int *oldparent = NULL;			/* element i was old static parent of i */

    /* for ease of reading, make alias of desttree, tree */
    tree = desttree;
    treecopy(matrix, tree, sourcetree, LVB_TRUE);

    /* get random branch but not root and not root's immediate descendant */
    do {
    	src = randpint(matrix->nbranches - 1);
    } while ((src == root) || (src == tree[root].left) || (src == tree[root].right));

    src_parent = tree[src].parent;
    lvb_assert(src_parent != UNSET);
    src_sister = getsister(tree, src);
    lvb_assert(src_sister != UNSET);

    /* get destination that is not source or its parent, sister or descendant
     * or the root */
    do {
    	dest = randpint(matrix->nbranches - 1);
    } while ((dest == src) || (dest == src_parent) || (dest == src_sister)
       || (dest == root) || is_descendant(tree, root, src, dest));

    /* excise source branch, leaving a damaged data structure */
    if (tree[src_parent].left == src) {
    	tree[src_parent].left = UNSET;
    }
    else if (tree[src_parent].right == src) {
    	tree[src_parent].right = UNSET;
    }
    else {
    	cr_bpnc(tree, src);
    }
    tree[src].parent = UNSET;

    /* fix data structure by "freeing" the excess branch */
    parents_par = tree[src_parent].parent;
    lvb_assert(parents_par != UNSET);
    if (tree[parents_par].left == src_parent) {
    	tree[parents_par].left = src_sister;
    }
    else {
    	tree[parents_par].right = src_sister;
    }
    tree[src_sister].parent = parents_par;

    excess_br = src_parent;	/* for ease of human understanding */
    nodeclear(tree, excess_br);

    /* make space at destination, re-using the excess branch */
    dest_parent = tree[dest].parent;
    if (tree[dest_parent].left == dest) {
    	orig_child = tree[dest_parent].left;
    	tree[dest_parent].left = excess_br;
    }
    else if (tree[dest_parent].right == dest) {
    	orig_child = tree[dest_parent].right;
    	tree[dest_parent].right = excess_br;
    }
    else {
    	crash("destination %ld is not a child of it's parent %ld\n", dest, dest_parent);
    }
    tree[excess_br].parent = dest_parent;
    tree[excess_br].left = dest;
    lvb_assert(orig_child != UNSET);
    tree[orig_child].parent = excess_br;

		if (oldparent == NULL) oldparent = alloc(matrix->nbranches * sizeof(int), "old parent alloc");

		int size = count(tree, src);
		int *arr=NULL;
		if (arr == NULL) arr = malloc(size * sizeof(*arr));

		int *mid_nodes=NULL;
		if (mid_nodes == NULL) mid_nodes = malloc(size * sizeof(*mid_nodes));
		int i = 0;

	/*XXXXX reroot source branch (only if size of subtree > than 2) */
	if (size > 2) {
		oldroot = src;
		for (current = 0; current < matrix->nbranches; current++)
	    	oldparent[current] = tree[current].parent;

     		addtoarray(tree, src, arr, 0);		
		
		 do {
		newroot = arr[randpint(size - 1)];
    } while (newroot == tree[oldroot].left || newroot == tree[oldroot].right);

		/* update the newroot */
		parnt = tree[newroot].parent;
		if (tree[parnt].left == newroot) sister = tree[parnt].right;
		else if (tree[parnt].right == newroot) sister = tree[parnt].left;
		
		tree[parnt].parent = previous;
		tree[parnt].left = oldparent[parnt];
		tree[parnt].right = newroot;
		
		oldsister = sister;
		previous = parnt;
		current = tree[parnt].left;
		lvb_assert(current != UNSET);

		/* loop for changing nodes between the newroot and oldroot */
		
		while (current != oldroot) {
			mid_nodes[i] = current;
			i++;
			
			parnt = oldparent[current];

			if (tree[current].left == previous) tempsister = tree[current].right;
			else if (tree[current].right == previous) tempsister = tree[current].left;
			lvb_assert(tempsister != UNSET);

			tree[current].right = oldsister;
			tree[current].parent = previous;
			tree[current].left = parnt;
			tree[parnt].parent = current;
			tree[oldsister].parent = current;
			
			oldsister = tempsister;
			previous = current;
			current = parnt;
		}
		
		/* updating the oldroot */
		tree[oldroot].parent = previous;
		if (tree[current].left == tree[oldroot].parent) tree[oldroot].left = oldsister;
    		else if (tree[current].right == tree[oldroot].parent) tree[oldroot].right = oldsister;
		tree[oldsister].parent = current;
		src = tree[newroot].parent;
	}

    /* add source branch to this new location */
    tree[excess_br].right = src;
    tree[src].parent = excess_br;

    /* ensure recalculation of lengths where necessary */
	if (size > 2) {
		make_dirty_below(matrix, tree, oldroot);
		if (i > 0) {
			int j;
			for (j = i-1; j >= 0; j--) {	
    				make_dirty_below(matrix, tree, mid_nodes[j]);
			}
		}
    		make_dirty_below(matrix, tree, src);
	}
	make_dirty_below(matrix, tree, excess_br);
	if (parents_par != root){
	 	make_dirty_below(matrix, tree, parents_par);
	}
	
	
	free(arr);
	free(mid_nodes);
} /* end mutate_tbr() */

/* Count total number of nodes of the tbr subtree */
int count(Branch *const tree, int current)
{

    if (current == UNSET)
        return 0;
    if ((tree[current].left == UNSET) && (tree[current].right == UNSET)) {
	return 1;
    } else {
	return count(tree, tree[current].left) + count(tree, tree[current].right);
    }
}

/* Generate an array of nodes from the tbr subtree */
int addtoarray(Branch *const tree, int current, int arr[], int i)
{
	if (current == UNSET)
	return 0;
	if (tree[current].left == UNSET && tree[current].right == UNSET) {
		arr[i] = current;
		i++;
	}
	if (tree[current].left != UNSET)
		i = addtoarray(tree, tree[current].left, arr, i);
	if (tree[current].right != UNSET)
		i = addtoarray(tree, tree[current].right, arr, i);
	return i;
}


long lvb_reroot(Dataptr restrict matrix, Branch *const barray, const int oldroot, const int newroot, Lvb_bool b_with_sset)
/* Change the root of the tree in barray from oldroot to newroot, which
 * must not be the same. Mark all internal nodes (everything but the leaves
 * and root) as "dirty". Return oldroot. */
{
	int current;							/* current branch */
	int parnt;								/* parent of current branch */
	int sister = UNSET;					/* sister of current branch */
	int previous;							/* previous branch */
	static int *oldparent = NULL;			/* element i was old static parent of i */

    /* check new root is a leaf but not the current root */
    lvb_assert(newroot < matrix->n);
    lvb_assert(newroot != oldroot);
    if (oldparent == NULL) oldparent = alloc(matrix->nbranches * sizeof(int), "old parent alloc");

    /* create record of parents as they are now */
    for (current = 0; current < matrix->nbranches; current++)
    	oldparent[current] = barray[current].parent;

    current = newroot;
    previous = UNSET;
    while (current != oldroot)
    {
		lvb_assert(current != UNSET);
		parnt = oldparent[current];		/* original parent */
		if (current == barray[parnt].left) sister = barray[parnt].right;
		else if (current == barray[parnt].right) sister = barray[parnt].left;
		else	/* error in tree structure */
			crash("internal error in function lvb_reroot(): current\n"
			 "branch %ld has old parent %ld, but old parent does not\n"
			 "have it as a child", current, parnt);
		barray[current].parent = previous;	/* now chld of prev. */

		/* make former parent the new left child, and former sister the
		 * new right child of the current branch */
		barray[current].left = parnt;
		barray[current].right = sister;
		barray[parnt].parent = current;
		barray[sister].parent = current;

		/* move towards original root, i.e. to original parent of
		 * current branch */
		previous = current;
		current = parnt;
    }

    /* former root is now a normal leaf, without descendants */
    barray[oldroot].left = UNSET;
    barray[oldroot].right = UNSET;

    if (b_with_sset){
    	for (current = matrix->n; current < matrix->nbranches; current++)
    		barray[current].sset[0] = 0U;
    }
    return oldroot;
} /* end lvb_reroot() */


int arbreroot(Dataptr matrix, Branch *const tree, const int oldroot)
/* Change tree's root arbitrarily, to a leaf other than oldroot.
 * Mark all nodes other than the leaves and root "dirty".
 * Return the number of the new root. */
{
	int newroot;		/* new root */

    /* find a leaf that is not the current root */
	int ugg_minus_1 = matrix->n - 1;
    do {
    	newroot = randpint(ugg_minus_1);
    } while (newroot == oldroot);

    lvb_reroot(matrix, tree, oldroot, newroot, LVB_TRUE);
    return newroot;

} /* end arbreroot() */

static int getsister(const Branch *const barray, const long branch)
/* return number of sister of branch branch in tree in barray, or UNSET if
 * branch has none */
{
    int parnt;		/* parent of current branch */

    parnt = barray[branch].parent;
    if (parnt == UNSET) return UNSET;
    if (branch == barray[parnt].left) return barray[parnt].right;
    else if (branch == barray[parnt].right) return barray[parnt].left;
    else	/* error in tree structure */
    {
		cr_bpnc(barray, branch);
		return 0;	/* NEVER reached but it shuts up compilers */
    }

} /* end getsister() */

int childadd(Branch *const tree, const int destination, const int newchild)
/* replace unset child of destination with newchild, and return
 * destination */
{
    if (tree[destination].right == UNSET) tree[destination].right = newchild;
    else if (tree[destination].left == UNSET) tree[destination].left = newchild;
    else	/* error: destination already has 2 children */
    	cr_chaf(tree, destination, newchild);

    tree[newchild].parent = destination;
    return destination;

} /* end childadd() */

static void cr_chaf(const Branch *const barray, const long destination, const long newchild)
/* crash because we want to add branch newchild to the children of
 * branch destination in tree in barray, but it already has two so
 * there is no room */
{
    crash("internal error in tree array %p: cannot make branch %ld a\n"
     "child of branch %ld since this already has 2 children (left is\n"
     "branch %ld, right is branch %ld)", (const void *) barray,
     newchild, destination, barray[destination].left,
     barray[destination].right);

} /* end cr_chaf() */

static void cr_bpnc(const Branch *const barray, const long branch)
/* crash because branch branch in tree in barray is not connected to
 * its parent, i.e. it is not a child of the branch it claims as
 * parent, according to that 'parent's' record of its own children */
{
    const long parnt = barray[branch].parent;	/* parent of branch */

    crash("internal error in tree array %p: branch record %ld says\n"
     "it has parent %ld, but branch record %ld has left child %ld\n"
     "and right child %ld", (const void *) barray, branch, parnt,
     parnt, barray[parnt].left, barray[parnt].right);

}	/* end cr_bpnc() */

Branch *mvBranch(long nbranches, Branch *const dest, const Branch *const src)
{
	long i;
	Lvb_bit_lentgh *tmp_sset;
	for(i = 0; i < nbranches; i++){
		tmp_sset = dest[i].sset;
		dest[i] = src[i];
		dest[i].sset = tmp_sset;
	}
	return dest;
}



void treecopy(Dataptr restrict matrix, Branch *const dest, const Branch *const src, Lvb_bool b_with_sset)
/* copy tree from src to dest; dest must be totally distinct from source
 * in memory, and have enough space; the approach used below may fail if
 * treealloc() is changed */
{
    long i;				/* loop counter */
    Lvb_bit_lentgh *tmp_sset;		/* temporary variable used in copy */
    unsigned char *src_statesets_all;	/* start of source's statesets */
    unsigned char *dest_statesets_all;	/* start of dest's statesets */

    if (b_with_sset){
		/* scalars */
		for (i = 0; i < matrix->nbranches; i++){
			tmp_sset = dest[i].sset;
			dest[i] = src[i];
			dest[i].sset = tmp_sset;	/* keep dest's stateset arrs for dest */
		}

		/* stateset arrays */
		src_statesets_all = ((unsigned char *) src) + matrix->nbranches * sizeof(Branch);
		dest_statesets_all = ((unsigned char *) dest) + matrix->nbranches * sizeof(Branch);
		memcpy(dest_statesets_all, src_statesets_all, matrix->nbranches * matrix->bytes);
    }else{
    	/* only the scalars */
		for (i = 0; i < matrix->nbranches; i++){
			tmp_sset = dest[i].sset;
			dest[i] = src[i];
			dest[i].sset = tmp_sset;
		}
    }

} /* end treecopy() */

void copy_sset(Dataptr restrict matrix, Objset *p_sset_1){

	long to_copy;
	for (long i = 0; i < matrix->nsets; i++){
		to_copy = sset_2[i].cnt * sizeof(int);
		if (p_sset_1[i].set == NULL){	// need to alloc memory
			p_sset_1[i].set = alloc(to_copy, "object set object arrays");

		}else if (p_sset_1[i].cnt != sset_2[i].cnt){
			p_sset_1[i].set = realloc(p_sset_1[i].set, to_copy);
			if (p_sset_1[i].set == NULL){
			     crash("out of memory: cannot increase allocation for best sset %ld elements", to_copy);
			}
		}
		memcpy(p_sset_1[i].set, sset_2[i].set, to_copy);
		p_sset_1[i].cnt = sset_2[i].cnt;
	}
}

void randtree(Dataptr matrix, Branch *const barray)
/* fill barray with a random tree, where barray[0] is the root; all branches
 * in this random tree are marked as "dirty" */
{
    Lvb_bool *leafmask;		/* LVB_TRUE where branch in array is a leaf */
    int *objnos;		/* element i is obj associated with branch i */

    treeclear(matrix, barray);
    leafmask = randtopology(matrix, barray, matrix->n);
    objnos = randleaf(matrix, barray, leafmask);
    tree_make_canonical(matrix, barray, objnos);

} /* end randtree() */

static void wherever_was_now_say(Dataptr restrict matrix, Branch *const barray, int was, int now)
{
    long branchno;			/* loop counter */

    for (branchno = 0; branchno < matrix->nbranches; branchno++){
		if(barray[branchno].parent == was) barray[branchno].parent = now;
		if(barray[branchno].left == was) barray[branchno].left = now;
		if(barray[branchno].right == was) barray[branchno].right = now;
    }

} /* end wherever_was_now_say() */

static void tree_make_canonical(Dataptr restrict matrix, Branch *const barray, int *objnos)
/* ensure that objects 0, 1, 2, ... n-1 are associated with branches 0, 1, 2,
 * ... n-1, respectively; objnos indicates for each branch the currently
 * assigned object or UNSET for internal branches */
{
	int i;				/* loop counter */
    int obj_no;			/* current object number */
    long n_lines = matrix->n;
    int impossible_1 = (int) matrix->nbranches;	/* an out-of-range branch index */
    int impossible_2 = (int) matrix->nbranches + 1;	/* an out-of-range branch index */
    int root = UNSET;			/* root branch index */
    Branch tmp_1, tmp_2;		/* temporary branches for swapping */
    unsigned char *ss0_start = (unsigned char *) barray[0].sset;	/* start of state set memory */
    Lvb_bool swap_made;			/* flag to indicate swap made */
    long tmp;				/* for swapping */

    do {
		swap_made = LVB_FALSE;
		for (i = 0; i < matrix->nbranches; i++) {
			obj_no = objnos[i];
			if ((obj_no != UNSET) && (obj_no != i)) {
				tmp_1 = barray[obj_no];
				wherever_was_now_say(matrix, barray, obj_no, impossible_1);
				tmp_2 = barray[i];
				wherever_was_now_say(matrix, barray, i, impossible_2);
				if (tmp_1.parent == i) tmp_1.parent = impossible_2;
				if (tmp_1.left == i) tmp_1.left = impossible_2;
				if (tmp_1.right == i) tmp_1.right = impossible_2;
				if (tmp_2.parent == obj_no) tmp_2.parent = impossible_1;
				if (tmp_2.left == obj_no) tmp_2.left = impossible_1;
				if (tmp_2.right == obj_no) tmp_2.right = impossible_1;
				barray[i] = tmp_1;
				barray[obj_no] = tmp_2;
				wherever_was_now_say(matrix, barray, impossible_1, i);
				wherever_was_now_say(matrix, barray, impossible_2, obj_no);
				tmp = objnos[i];
				objnos[i] = objnos[obj_no];
				objnos[obj_no] = tmp;
				swap_made = LVB_TRUE;
			}
		}
    } while (swap_made == LVB_TRUE);

    /* patch up assignment of sset memory to prevent trouble in treecopy() */
    for (i = 0; i < matrix->nbranches; i++){
    	barray[i].sset = (Lvb_bit_lentgh *) (ss0_start + i * matrix->bytes);
    }

    for (i = 0; i < n_lines; i++) {
		if (barray[i].parent == UNSET){
			lvb_assert(root == UNSET);
			root = i;
		}
    }

    if (root != 0) {
    	lvb_reroot(matrix, barray, root, 0, LVB_TRUE);
    }

    /* check things didn't go haywire */
    for (i = 0; i < n_lines; i++) {
    	lvb_assert(objnos[i] == i);
    }
    for (i = n_lines; i < matrix->nbranches; i++) {
    	lvb_assert(objnos[i] == UNSET);
    }

} /* end tree_make_canonical() */

Branch *treealloc(Dataptr restrict matrix, Lvb_bool b_with_sset)
/* Return array of nbranches branches with scalars all UNSET, and all
 * statesets allocated for m characters but marked "dirty". Crash
 * verbosely if impossible. Memory is allocated once only, as a contiguous
 * block for the branch data structures followed by all their statesets.
 * So, to deallocate the tree, call the standard library function free()
 * ONCE ONLY, passing it the address of the first branch struct. If this
 * allocation approach is changed, be sure to change treecopy() too. */
{
    Branch *barray;						/* tree */
    unsigned char *barray_uchar_star;	/* tree as unsigned char */
    unsigned char *ss0_start;			/* start of first stateset */
    long i;								/* loop counter */

    lvb_assert(matrix->nbranches >= MIN_BRANCHES);
    lvb_assert(matrix->nbranches <= MAX_BRANCHES);

    if (b_with_sset) barray = alloc(matrix->tree_bytes, "tree with statesets");
    else{ /* don't need to do anything else */
    	barray = alloc(matrix->tree_bytes_without_sset, "tree without statesets");
    	return barray;
    }

    barray_uchar_star = (unsigned char *) barray;
    ss0_start = barray_uchar_star + matrix->nbranches * sizeof(Branch);

    /* crash if state set memory is misaligned for uint32_t */
    lvb_assert(((intptr_t) ss0_start % NIBBLE_WIDTH) == 0);
    lvb_assert((matrix->bytes % NIBBLE_WIDTH) == 0);

    for (i = 0; i < matrix->nbranches; i++){
    	barray[i].sset = (Lvb_bit_lentgh *) (ss0_start + i * matrix->bytes);
    	*barray[i].sset = 0U;  /* make durty */
    }

    /*make_dirty_tree(matrix, barray);  */
    return barray;

} /* end treealloc() */

static Lvb_bool *randtopology(Dataptr matrix, Branch *const barray, const long nobjs)
/* fill barray with tree of random topology, where barray[0] is root;
 * return static array where element i is LVB_TRUE if barray[i] is a
 * leaf or LVB_FALSE if it is not; this array will be overwritten in
 * subsequent calls */
{
    long i;		/* loop counter */
    long leaves = 0;	/* number of leaves */
    int nextfree = 0;	/* next unused element of barray */
    int togrow;	/* random candidate for sprouting */
    static Lvb_bool isleaf[MAX_BRANCHES];	/* return value */

    lvb_assert(nobjs == matrix->n);

    /* clear the leaf mask */
    for (i = 0; i < matrix->nbranches; i++)
    	isleaf[i] = LVB_FALSE;

    /* start with initial tree of 3 leaves */
    barray[0].parent = UNSET;
    isleaf[nextfree++] = LVB_TRUE;
    barray[0].left = nextfree;
    barray[nextfree].parent = 0;
    isleaf[nextfree++] = LVB_TRUE;
    barray[0].right = nextfree;
    barray[nextfree].parent = 0;
    isleaf[nextfree++] = LVB_TRUE;
    leaves = 3;

    /* sprout! */
    while(leaves < nobjs)
    {
		do	/* select a random leaf other than the root */
		{
			togrow = 1 + randpint(nextfree - 2);
		} while (isleaf[togrow] == LVB_FALSE);
		/* left child */
		barray[togrow].left = nextfree;
		barray[nextfree].parent = togrow;
		isleaf[nextfree++] = LVB_TRUE;
		/* right child */
		barray[togrow].right = nextfree;
		barray[nextfree].parent = togrow;
		isleaf[nextfree++] = LVB_TRUE;
		/* other updates */
		isleaf[togrow] = LVB_FALSE;
		leaves++;
    }

    return isleaf;

} /* end randtopology() */

static int *randleaf(Dataptr restrict matrix, Branch *const barray, const Lvb_bool *const leafmask)
/* randomly assign objects numbered 0 to objs - 1 to leaves of tree in
 * barray; leaves in barray must be indicated by corresponding
 * LVB_TRUEs in leafmask; returns static array of object numbers, where
 * elements 0..nbranches give the object associated with branches 0..nbranches
 * respectively (UNSET for internal branches); this array will be overwritten
 * in subsequent calls */
{
    long assigned = 0;				/* for safety: should == objs at end */
    int candidate;			/* random object */
    long i;							/* loop counter */
    static Lvb_bool *p_used; 		/* element i LVB_TRUE if object i has leaf */
    static int *p_objnos = NULL;	/* object associated with branches */

    if (p_objnos == NULL) p_objnos = alloc(matrix->nbranches * sizeof(int), "Alloc objnos");
    if (p_used == NULL) p_used = alloc(matrix->n * sizeof(Lvb_bool), "Alloc p_used");

    /* clear 'used' array */
    // for (i = 0; i < matrix->n; i++) p_used[i] = LVB_FALSE;
    memset(p_used, 0, matrix->n * sizeof(Lvb_bool));

    /* clear object nos array, defaulting to internal branch, i.e. UNSET */
    //for (i = 0; i < matrix->nbranches; i++) p_objnos[i] = UNSET;
    memset(p_objnos, 0xff, matrix->nbranches * sizeof(int));

    /* assign an object to every leaf */
    for (i = 0; i < matrix->nbranches; i++){
		if (leafmask[i] == LVB_TRUE)	/* leaf, requires object */
		{
			do	/* get a new object number */
			{
				candidate = randpint(matrix->n - 1);
			} while(p_used[candidate] == LVB_TRUE);
			/* assign object to leaf */
			p_objnos[i] = candidate;
			p_used[candidate] = LVB_TRUE;
			assigned++;
		}
    }

    lvb_assert(assigned == matrix->n);
    return p_objnos;

} /* end randleaf() */

void treeswap(Branch **const tree1, long *const root1,
 Branch **const tree2, long *const root2)
/* swap trees pointed to by tree1 and tree2; also swap records of their
 * roots, as pointed to by root1 and root2 */
{
    long tmproot;	/* temporary value-holder for swapping roots */
    Branch *tmptree;	/* temporary value-holder for swapping trees */

    /* swap pointers to branch arrays */
    tmptree = *tree1;
    *tree1 = *tree2;
    *tree2 = tmptree;

    /* swap roots */
    tmproot = *root1;
    *root1 = *root2;
    *root2 = tmproot;

} /* end treeswap() */

void treedump(Dataptr matrix, FILE *const stream, const Branch *const tree, Lvb_bool b_with_sset)
/* send tree as table of integers to file pointed to by stream */
{
    long i;				/* loop counter */
    long j;				/* loop counter */

    fprintf(stream, "Branch\tParent\tLeft\tRight\tChanges\tDirty\tSset_arr\tSsets\n");
    for (i = 0; i < matrix->nbranches; i++) {
    	fprintf(stream, "%ld\t%u\t%u\t%u\t%u", i, tree[i].parent, tree[i].left, tree[i].right, tree[i].changes);
    	if (tree[i].sset[0] == 0U) fprintf(stream, "\tyes");
    	else fprintf(stream, "\tno");
    	if (b_with_sset){
			fprintf(stream, "\t%p", (void *) tree[i].sset);
			for (j = 0; j < matrix->nwords; j++){
				fprintf(stream, "\t0%o", (unsigned) tree[i].sset[j]);
			}
    	}
    	else{
    		fprintf(stream, "\tversion without sset");
    	}
    	fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
    if(ferror(stream) != 0)
	cr_uxe(stream, "dumping tree");

} /* end treedump() */

void treedump_screen(Dataptr matrix, const Branch *const tree)
/* send tree as table of integers to file pointed to by stream */
{
    long i;				/* loop counter */

    printf("Branch\tParent\tLeft\tRight\tChanges\tDirty\n");
    for (i = 0; i < matrix->nbranches; i++) {
    	printf("%ld\t%u\t%u\t%u\t%u", i, tree[i].parent, tree[i].left, tree[i].right, tree[i].changes);
    	if (tree[i].sset == NULL) printf("\tNULL\n");
    	else if (tree[i].sset[0] == 0U) printf("\tyes\n");
    	else printf("\tno\n");
    }
    printf("\n");

} /* end treedump() */

static void cr_uxe(FILE *const stream, const char *const msg)
/* crash because of problem on file stream, with message consisting of
 * "FATAL ERROR: ", then "file error " if the error indicator for
 * stream is set or "unexpected end of file " if not, then "when ",
 * followed by msg */
{
    if (ferror(stream))
	crash("file error when %s", msg);
    else
	crash("unexpected end of file when %s", msg);

} /* end cr_uxe */

void lvb_treeprint (Dataptr matrix, FILE *const stream, const Branch *const barray, const long root)
/* print tree in barray (of root root) in bracketed text form to stream stream,
 * in unrooted form */
{
    ur_print(matrix, stream, barray, root);
} /* end lvb_treeprint() */

static void ur_print(Dataptr matrix, FILE *const stream, const Branch *const barray, const long root)
/* send tree in barray, of root root, to file pointed to by stream in
 * unrooted form */
{
    long obj;					/* current object */
    static Lvb_bool doneabsroot = LVB_FALSE;	/* have output root */
    static Lvb_bool usecomma;			/* output clade sep. */
    char *tmp_title;				/* temporary string */

    obj = root;

    if (doneabsroot == LVB_FALSE)	/* print whole tree */
    {
		/* start tree */
		tmp_title = alloc(strlen(matrix->rowtitle[obj]) + 1, "temp. title");
		strcpy(tmp_title, matrix->rowtitle[obj]);
		while(tmp_title[strlen(tmp_title) - 1] == ' '){
			tmp_title[strlen(tmp_title) - 1] = '\0';
		}
		fprintf(stream, "(%s", tmp_title);
		free(tmp_title);    /* VERY LOCAL dynamic heap memory */
		usecomma = LVB_TRUE;
		doneabsroot = LVB_TRUE;

		ur_print(matrix, stream, barray, barray[root].left);
		ur_print(matrix, stream, barray, barray[root].right);

		/* end tree */
		fprintf(stream, ");\n");
		if (ferror(stream))
			crash("file error when writing unrooted tree");

		/* clean up for next call */
		usecomma = LVB_FALSE;
		doneabsroot = LVB_FALSE;
    }
    else	/* print remainder of tree */
    {
		if (usecomma == LVB_TRUE) fprintf(stream, "%s", CLADESEP);
		if (root < matrix->n)	/* leaf */
		{
			tmp_title = alloc(strlen(matrix->rowtitle[obj]) + 1, "temp. title");
			strcpy(tmp_title, matrix->rowtitle[obj]);
			while(tmp_title[strlen(tmp_title) - 1] == ' '){
				tmp_title[strlen(tmp_title) - 1] = '\0';
			}
			fprintf(stream, "%s", tmp_title);
			free(tmp_title);	/* VERY LOCAL dynamic heap memory */
			usecomma = LVB_TRUE;
		}
		else
		{
			fprintf(stream, "(");
			usecomma = LVB_FALSE;
			ur_print(matrix, stream, barray, barray[root].left);
			ur_print(matrix, stream, barray, barray[root].right);
			fputc(')', stream);
			usecomma = LVB_TRUE;
		}
    }

} /* end ur_print() */



long treecmp(Dataptr matrix, Objset *sset_1, const Branch *const tree_2, Lvb_bool b_First)
/* return 0 if the topology of tree_1 (of root root_1) is the same as
 * that of tree_2 (of root root_2), or non-zero if different */
{
//	b_First = LVB_TRUE;
    if (b_First == LVB_TRUE) makesets(matrix, tree_2, 0 /* always root zero */);
    return setstcmp(matrix, sset_1, sset_2 /* this one is the static */);

} /* end treecmp() */

static long setstcmp(Dataptr matrix, Objset *const oset_1, Objset *const oset_2 /* this one is the static */)
/* return 0 if the same sets of objects are in oset_1 and oset_2,
 * and non-zero otherwise */
{
    long i;		/* loop counter */

    /* compare the set arrays */
/*    printf("########################\n######################\n");
    for (i = 0; i < matrix->nsets; i++){
    	for (int j = 0; j < oset_1[i].cnt; j++){
    	    printf("Set:%ld   cnt:%d   original:%d    comapre:%d\n", i, j, oset_1[i].set[j], oset_2[i].set[j]);
    	}
    }*/
    for (i = 0; i < matrix->nsets; i++){
    	if (oset_1[i].cnt != oset_2[i].cnt) return 1;
    	if (memcmp(oset_1[i].set, oset_2[i].set, sizeof(int) * oset_1[i].cnt) != 0) return 1;
    }
    return 0;
} /* end setstcmp() */


long setstcmp_with_sset2(Dataptr matrix, Objset *const oset_1)
/* return 0 if the same sets of objects are in oset_1 and oset_2,
 * and non-zero otherwise */
{
    long i;		/* loop counter */

    for (i = 0; i < matrix->nsets; i++){
    	if (oset_1[i].cnt != sset_2[i].cnt) return 1;
    	if (memcmp(oset_1[i].set, sset_2[i].set, sizeof(int) * oset_1[i].cnt) != 0) return 1;
    }
    return 0;
} /* end setstcmp() */


void dump_stack_to_screen(Dataptr matrix, Treestack *sp){
	for (int i = 0; i < sp->next; i++){
		printf("Stack number: %d\n", i);
		dump_objset_to_screen(matrix, sp->stack[i].p_sset);
	}
}


void dump_objset_to_screen(Dataptr matrix, Objset *oset_1){
	printf("################\n##################\n");
	for (int i = 0; i < matrix->nsets; i++){
		printf("%d    %d    ", i, oset_1[i].cnt);
		for (int x = 0; x < oset_1[i].cnt; x++) printf("%d   ", oset_1[i].set[x]);
		printf("\n");
	}
	printf("\n");
}

void dump_objset_to_screen_sset_2(Dataptr matrix){
	dump_objset_to_screen(matrix, sset_2);
}

void sort_array(int *p_array, int n_left, int n_rigth){
	int l_hold, r_hold;
	int l_pivot = *(p_array + ((n_left + n_rigth) / 2));
	int l_temp;

	l_hold = n_left;		//i=l;
	r_hold = n_rigth;     //j=r;
	do {
		while (p_array[l_hold] < l_pivot) l_hold++;
		while (p_array[r_hold] > l_pivot) r_hold--;

		if (l_hold <= r_hold){
			l_temp = p_array[l_hold];
			p_array[l_hold] = p_array[r_hold];
			p_array[r_hold] = l_temp;
			l_hold++;
			r_hold--;
		}
	} while (l_hold < r_hold);
	if (n_left < r_hold) sort_array(p_array, n_left, r_hold);
	if (l_hold < n_rigth) sort_array(p_array, l_hold, n_rigth);
}


static void sort(Dataptr matrix, Objset *const oset_2, const long nels)
/* sort the nels object sets in oset so that each is in order, and sort oset so
 * that the sets themselves are in order of size and content */
{
    /* first sort each set member list */
	omp_set_dynamic(0);	  /* disable dinamic threathing */
	#pragma omp parallel for num_threads(matrix->n_threads_getplen)
    	for (long i = 0; i < nels; i++)
    		sort_array(oset_2[i].set, 0, oset_2[i].cnt - 1);

   	/* now sort the arrays of sets by size and content */
   	qsort(oset_2, (size_t) nels, sizeof(Objset), osetcmp);
} /* end sort() */


static int osetcmp(const void *oset1, const void *oset2)
/* comparison function for object sets (type Objset):
 * return negative if *oset1 is a smaller set of objects than *oset2 or
 * is the same size but with a list of elements that compares lower;
 * return positive if *ostet1 is bigger or the same size but with a
 * list of elements that compares higher; return 0 if they are the
 * same; N.B. the object numbers must be in numerical order within the
 * sets */
{
    long i;						/* loop cntr */
    const Objset loset_1 = *((const Objset *) oset1);	/* typed */
    const Objset loset_2 = *((const Objset *) oset2);	/* typed */

    /* sets of different size differ */
    if (loset_1.cnt < loset_2.cnt) return -1;
    else if (loset_1.cnt > loset_2.cnt) return +1;

    /* if we reach here, sets are equal size, so we see if sets'
     * contents differ */
    for (i = 0; i < loset_1.cnt; i++){
    	if (loset_1.set[i] < loset_2.set[i]) return -1;
    	else if (loset_1.set[i] > loset_2.set[i]) return +1;
    }

    /* if we reach here, really the sets are the same */
    return 0;

} /* end osetcmp() */

void makesets(Dataptr matrix, const Branch *const tree_2, const int root)
/* fill static sset_1 and static sset_2 with arrays of object sets for
 * tree_1 and tree_2 (of root_1 and root_2 respectively), and return
 * the extent of each array;
 * the trees must have the same object in the root branch;
 * arrays will be overwritten on subsequent calls */
{
    if (sset_2[0].set == NULL){	/* first call, allocate memory  to the static sset_2*/
		ssarralloc(matrix, sset_2);
    }

    fillsets(matrix, sset_2, tree_2, root);
    sort(matrix, sset_2, matrix->nsets);
} /* end makesets() */

static void ssarralloc(Dataptr matrix, Objset *nobjset_2)
/* Fill nobjset[0..nsets-1] with pointers each pointing to newly
 * allocated space for setsize objects; assumes nobjset points to the
 * first element of an array with at least nsets elements. */
{
    long i; 	/* loop counter */
    for (i = 0; i < matrix->nsets; i++){
    	nobjset_2[i].set = alloc(matrix->mssz * sizeof(int), "object set object arrays");
    	nobjset_2[i].cnt = UNSET;
    }

} /* end ssarralloc() */

static void fillsets(Dataptr matrix, Objset *const sstruct, const Branch *const tree, const int root)
/* fill object sets in sstruct with all sets of objects in tree tree,
 * descended from but not including root and not including sets of one
 * object */
{
    static int i = UNSET;	/* current set being filled */

    if (i == UNSET)	/* not a recursive call */
    {
		i = 0;

		/* avoid generating sets for true root and leaves */
		if (tree[root].left >= matrix->n)	/* interior */
			fillsets(matrix, sstruct, tree, tree[root].left);
		if (tree[root].right >= matrix->n)	/* interior */
			fillsets(matrix, sstruct, tree, tree[root].right);

		i = UNSET;	/* clean up for next non-recursive call */
		return;
    }
	if (tree[root].left != UNSET)	/* not leaf */
	{
		getobjs(matrix, tree, root, sstruct[i].set, &sstruct[i].cnt);
		i++;
		fillsets(matrix, sstruct, tree, tree[root].left);
		fillsets(matrix, sstruct, tree, tree[root].right);
		return;
	}

} /* end fillsets */

static void getobjs(Dataptr matrix, const Branch *const barray, const int root, int *const objarr, int *const cnt)
/* fill objarr (which must be large enough) with numbers of all objects
 * in the tree in barray in the clade starting at branch root;
 * fill the number pointed to by cnt with the number of objects found
 * (i.e. the number of elements written to objarr) */
{
    *cnt = 0;
    realgetobjs(matrix, barray, root, objarr, cnt);

} /* end getobjs() */

static void realgetobjs(Dataptr matrix, const Branch *const barray, const int root, int *const objarr, int *const cnt)
/* fill objarr (which must be large enough) with numbers of all objects
 * in the tree in barray in the clade starting at branch root;
 * fill the number pointed to by cnt, which must initially be zero,
 * with the number of objects found (i.e. the number of elements
 * written to objarr); this function should not be called from anywhere
 * except getobjs(), which is a safer interface */
{
    if (root < matrix->n) {
		objarr[*cnt] = root;
		++(*cnt);
    }
    else {
		if (barray[root].left != UNSET)
			realgetobjs(matrix, barray, barray[root].left, objarr, cnt);
		if (barray[root].right != UNSET)
			realgetobjs(matrix, barray, barray[root].right, objarr, cnt);
    }

} /* end realgetobjs() */


void ss_init(Dataptr matrix, Branch *tree, Lvb_bit_lentgh **enc_mat)
/* copy m states from enc_mat to the stateset arrays for the leaves in tree,
 * including padding at the end; the nth entry in enc_mat is assumed to be the
 * encoded state sets for object number n in the tree; non-leaf branches in the
 * tree are marked "dirty"; the root branch struct is marked "clean" since it
 * is also a terminal */
{
    long i;			/* loop counter */
    for (i = 0; i < matrix->n; i++) memcpy(tree[i].sset, enc_mat[i], matrix->bytes);
    for (i = matrix->n; i < matrix->nbranches; i++) tree[i].sset[0] = 0U;

} /* end ss_init() */

