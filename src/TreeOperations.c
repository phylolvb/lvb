/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== TreeOperations.c - tree operations ========== */

#include "TreeOperations.h"

void nodeclear(TREESTACK_TREE_NODES *const BranchArray, const long brnch)
/* Initialize all scalars in branch brnch to UNSET or zero as appropriate,
 * and mark it "dirty" */
{
	BranchArray[brnch].left = UNSET;
	BranchArray[brnch].right = UNSET;
	BranchArray[brnch].parent = UNSET;
	BranchArray[brnch].changes = UNSET;
	BranchArray[brnch].sitestate[0] = 0U; /* "dirty" */

} /* end nodeclear() */

long tree_bytes(Dataptr restrict MSA)
/* return bytes required for contiguous allocation of a tree for the data
 * accessible by MSA, if branches and their statesets are allocated
 * as one contiguous array */
{
	/* return in bytes */
	return (MSA->numberofpossiblebranches * sizeof(TREESTACK_TREE_NODES)) + (MSA->numberofpossiblebranches * MSA->bytes);
} /* end tree_bytes() */

long tree_bytes_without_sitestate(Dataptr restrict MSA)
/* return bytes required for contiguous allocation of a tree for the data
 * accessible by MSA, */
{
	/* return in bytes */
	return (MSA->numberofpossiblebranches * sizeof(TREESTACK_TREE_NODES));
} /* end tree_bytes() */

void treeclear(Dataptr restrict MSA, TREESTACK_TREE_NODES *const BranchArray)
/* clear all branches in array BranchArray, on the assumption that its size fits
 * the data MSA; mark all branches dirty */
{
	long i; /* loop counter */
	for (i = 0; i < MSA->numberofpossiblebranches; i++)
		nodeclear(BranchArray, i);

} /* end treeclear() */

static void make_dirty_below(Dataptr restrict MSA, TREESTACK_TREE_NODES *tree, long dirty_node)
/* mark nodes "dirty" from branch dirty_node, which must not be the root,
 * down to (but not including) the root branch of the tree tree; the true
 * root lies outside the LVB tree data structure so cannot be marked
 * dirty, but will always be dirty after any rearrangement */
{

	lvb_assert(dirty_node >= MSA->n); /* not leaf/root */
	lvb_assert(tree[dirty_node].parent != UNSET);
	do
	{
		tree[dirty_node].sitestate[0] = 0U; /* " make dirty" */
		dirty_node = tree[dirty_node].parent;
	} while (tree[dirty_node].parent != UNSET);

} /* end make_dirty_below() */

void mutate_deterministic(Dataptr restrict MSA, TREESTACK_TREE_NODES *const desttree,
						  const TREESTACK_TREE_NODES *const sourcetree, long root, long p, Lvb_bool left)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a change in topology, the change being caused by nearest
 * neighbour interchange (NNI) rearrangement at branch p, involving the
 * right node if right is LVB_TRUE, otherwise the left node; N.B. code
 * is largely copied from mutate_nni() */
{
	TREESTACK_TREE_NODES *tree;
	long u, v, a, b, c;

	lvb_assert(p != root);
	lvb_assert(p >= MSA->n);

	/* for ease of reading, make alias of desttree, tree */
	tree = desttree;
	treecopy(MSA, tree, sourcetree, LVB_TRUE);

	u = p;
	v = tree[u].parent;
	a = tree[u].left;
	b = tree[u].right;

	if (tree[v].left == u)
		c = tree[v].right;
	else
		c = tree[v].left;

	if (left != LVB_TRUE)
	{
		if (tree[v].left == u)
			tree[v].right = b;
		else
			tree[v].left = b;
		tree[u].left = a;
		tree[u].right = c;
		tree[a].parent = tree[c].parent = u;
		tree[b].parent = v;
	}
	else
	{
		if (tree[v].left == u)
			tree[v].right = a;
		else
			tree[v].left = a;
		tree[u].left = b;
		tree[u].right = c;
		tree[b].parent = tree[c].parent = u;
		tree[a].parent = v;
	}

	make_dirty_below(MSA, tree, u);

} /* end mutate_nni() */

void mutate_nni(Dataptr restrict MSA, TREESTACK_TREE_NODES *const desttree, const TREESTACK_TREE_NODES *const sourcetree, long root)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a random change in topology, the change being caused by nearest
 * neighbour interchange (NNI) rearrangement; N.B. the code is mostly
 * the same in mutate_deterministic() */
{
	TREESTACK_TREE_NODES *tree;
	long u, v, a, b, c;

	/* for ease of reading, make alias of desttree, tree */
	tree = desttree;
	treecopy(MSA, tree, sourcetree, LVB_TRUE);

	/* get a random internal branch */
	u = randpint(MSA->numberofpossiblebranches - MSA->n - 1) + MSA->n;
	v = tree[u].parent;
	a = tree[u].left;
	b = tree[u].right;

	if (tree[v].left == u)
		c = tree[v].right;
	else
		c = tree[v].left;

	if (uni() < 0.5)
	{
		if (tree[v].left == u)
			tree[v].right = b;
		else
			tree[v].left = b;
		tree[u].left = a;
		tree[u].right = c;
		tree[a].parent = tree[c].parent = u;
		tree[b].parent = v;
	}
	else
	{
		if (tree[v].left == u)
			tree[v].right = a;
		else
			tree[v].left = a;
		tree[u].left = b;
		tree[u].right = c;
		tree[b].parent = tree[c].parent = u;
		tree[a].parent = v;
	}

	make_dirty_below(MSA, tree, u);

} /* end mutate_nni() */

static Lvb_bool is_descendant(TREESTACK_TREE_NODES *tree, long root, long ancestor, long candidate)
/* return LVB_TRUE if candidate is among the branches in the clade descending
from ancestor, LVB_FALSE otherwise */
{
	Lvb_bool val = LVB_FALSE;		   /* return value */
	long par = tree[candidate].parent; /* current parent */
	long newpar;					   /* next parent */

	while (par != UNSET)
	{
		if (par == ancestor)
		{
			val = LVB_TRUE;
			par = UNSET; /* exit the loop */
		}
		else
		{
			newpar = tree[par].parent;
			par = newpar;
		}
	}
	return val;

} /* end is_descendant() */

void mutate_spr(Dataptr restrict MSA, TREESTACK_TREE_NODES *const desttree, const TREESTACK_TREE_NODES *const sourcetree, long root)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a random change in topology, the change being caused by subtree
 * pruning and regrafting (SPR) rearrangement */
{
	long src;					/* branch to move */
	long dest;					/* destination of branch to move */
	long dest_parent;			/* parent of destination branch */
	long src_parent;			/* parent of branch to move */
	long excess_br;				/* branch temporarily excised */
	long orig_child = UNSET;	/* original child of destination */
	long parents_par;			/* parent of parent of br. to move */
	long src_sister;			/* sister of branch to move */
	TREESTACK_TREE_NODES *tree; /* destination tree */

	/* for ease of reading, make alias of desttree, tree */
	tree = desttree;
	treecopy(MSA, tree, sourcetree, LVB_TRUE);

	/* get random branch but not root and not root's immediate descendant */
	do
	{
		src = randpint(MSA->numberofpossiblebranches - 1);
	} while ((src == root) || (src == tree[root].left) || (src == tree[root].right));

	src_parent = tree[src].parent;
	lvb_assert(src_parent != UNSET);
	src_sister = getsister(tree, src);
	lvb_assert(src_sister != UNSET);

	/* get destination that is not source or its parent, sister or descendant
	 * or the root */
	do
	{
		dest = randpint(MSA->numberofpossiblebranches - 1);
	} while ((dest == src) || (dest == src_parent) || (dest == src_sister) || (dest == root) || is_descendant(tree, root, src, dest));

	/* excise source branch, leaving a damaged data structure */
	if (tree[src_parent].left == src)
	{
		tree[src_parent].left = UNSET;
	}
	else if (tree[src_parent].right == src)
	{
		tree[src_parent].right = UNSET;
	}
	else
	{
		cr_bpnc(tree, src);
	}
	tree[src].parent = UNSET;

	/* fix data structure by "freeing" the excess branch */
	parents_par = tree[src_parent].parent;
	lvb_assert(parents_par != UNSET);
	if (tree[parents_par].left == src_parent)
	{
		tree[parents_par].left = src_sister;
	}
	else
	{
		tree[parents_par].right = src_sister;
	}
	tree[src_sister].parent = parents_par;

	excess_br = src_parent; /* for ease of human understanding */
	nodeclear(tree, excess_br);

	/* make space at destination, re-using the excess branch */
	dest_parent = tree[dest].parent;
	if (tree[dest_parent].left == dest)
	{
		orig_child = tree[dest_parent].left;
		tree[dest_parent].left = excess_br;
	}
	else if (tree[dest_parent].right == dest)
	{
		orig_child = tree[dest_parent].right;
		tree[dest_parent].right = excess_br;
	}
	else
	{
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
	make_dirty_below(MSA, tree, excess_br);
	if (parents_par != root)
	{
		make_dirty_below(MSA, tree, parents_par);
	}
} /* end mutate_spr() */

void mutate_tbr(Dataptr restrict MSA, TREESTACK_TREE_NODES *const desttree, const TREESTACK_TREE_NODES *const sourcetree, long root)
/* make a copy of the tree sourcetree (of root root) in desttree,
 * with a random change in topology, the change being caused by subtree
 * pruning and regrafting (SPR) rearrangement */
{
	long src;					/* branch to move */
	long dest;					/* destination of branch to move */
	long dest_parent;			/* parent of destination branch */
	long src_parent;			/* parent of branch to move */
	long excess_br;				/* branch temporarily excised */
	long orig_child = UNSET;	/* original child of destination */
	long parents_par;			/* parent of parent of br. to move */
	long src_sister;			/* sister of branch to move */
	TREESTACK_TREE_NODES *tree; /* destination tree */

	long oldroot;
	long current;		 /* current branch */
	long parnt;			 /* parent of current branch */
	long sister = UNSET; /* sister of current branch */
	long oldsister = UNSET;
	long tempsister = UNSET;
	long previous = UNSET; /* previous branch */
	long newroot;
	static int *oldparent = NULL; /* element i was old static parent of i */

	/* for ease of reading, make alias of desttree, tree */
	tree = desttree;
	treecopy(MSA, tree, sourcetree, LVB_TRUE);

	/* get random branch but not root and not root's immediate descendant */
	do
	{
		src = randpint(MSA->numberofpossiblebranches - 1);
	} while ((src == root) || (src == tree[root].left) || (src == tree[root].right));

	src_parent = tree[src].parent;
	lvb_assert(src_parent != UNSET);
	src_sister = getsister(tree, src);
	lvb_assert(src_sister != UNSET);

	/* get destination that is not source or its parent, sister or descendant
	 * or the root */
	do
	{
		dest = randpint(MSA->numberofpossiblebranches - 1);
	} while ((dest == src) || (dest == src_parent) || (dest == src_sister) || (dest == root) || is_descendant(tree, root, src, dest));

	/* excise source branch, leaving a damaged data structure */
	if (tree[src_parent].left == src)
	{
		tree[src_parent].left = UNSET;
	}
	else if (tree[src_parent].right == src)
	{
		tree[src_parent].right = UNSET;
	}
	else
	{
		cr_bpnc(tree, src);
	}
	tree[src].parent = UNSET;

	/* fix data structure by "freeing" the excess branch */
	parents_par = tree[src_parent].parent;
	lvb_assert(parents_par != UNSET);
	if (tree[parents_par].left == src_parent)
	{
		tree[parents_par].left = src_sister;
	}
	else
	{
		tree[parents_par].right = src_sister;
	}
	tree[src_sister].parent = parents_par;

	excess_br = src_parent; /* for ease of human understanding */
	nodeclear(tree, excess_br);

	/* make space at destination, re-using the excess branch */
	dest_parent = tree[dest].parent;
	if (tree[dest_parent].left == dest)
	{
		orig_child = tree[dest_parent].left;
		tree[dest_parent].left = excess_br;
	}
	else if (tree[dest_parent].right == dest)
	{
		orig_child = tree[dest_parent].right;
		tree[dest_parent].right = excess_br;
	}
	else
	{
		crash("destination %ld is not a child of it's parent %ld\n", dest, dest_parent);
	}
	tree[excess_br].parent = dest_parent;
	tree[excess_br].left = dest;
	lvb_assert(orig_child != UNSET);
	tree[orig_child].parent = excess_br;

	if (oldparent == NULL)
		oldparent = (int *)alloc(MSA->numberofpossiblebranches * sizeof(int), "old parent alloc");

	int size = count(tree, src);
	int *arr = NULL;
	if (arr == NULL)
		arr = (int *)malloc(size * sizeof(*arr));

	int *mid_nodes = NULL;
	if (mid_nodes == NULL)
		mid_nodes = (int *)malloc(size * sizeof(*mid_nodes));
	int i = 0;

	/*reroot source branch (only if size of subtree > than 2) */
	if (size > 2)
	{
		oldroot = src;
		for (current = 0; current < MSA->numberofpossiblebranches; current++)
			oldparent[current] = tree[current].parent;

		addtoarray(tree, src, arr, 0);

		do
		{
			newroot = arr[randpint(size - 1)];
		} while (newroot == tree[oldroot].left || newroot == tree[oldroot].right);

		/* update the newroot */
		parnt = tree[newroot].parent;
		if (tree[parnt].left == newroot)
			sister = tree[parnt].right;
		else if (tree[parnt].right == newroot)
			sister = tree[parnt].left;

		tree[parnt].parent = previous;
		tree[parnt].left = oldparent[parnt];
		tree[parnt].right = newroot;

		oldsister = sister;
		previous = parnt;
		current = tree[parnt].left;
		lvb_assert(current != UNSET);

		/* loop for changing nodes between the newroot and oldroot */

		while (current != oldroot)
		{
			mid_nodes[i] = current;
			i++;

			parnt = oldparent[current];

			if (tree[current].left == previous)
				tempsister = tree[current].right;
			else if (tree[current].right == previous)
				tempsister = tree[current].left;
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
		if (tree[current].left == tree[oldroot].parent)
			tree[oldroot].left = oldsister;
		else if (tree[current].right == tree[oldroot].parent)
			tree[oldroot].right = oldsister;
		tree[oldsister].parent = current;
		src = tree[newroot].parent;
	}

	/* add source branch to this new location */
	tree[excess_br].right = src;
	tree[src].parent = excess_br;

	/* ensure recalculation of lengths where necessary */
	if (size > 2)
	{
		make_dirty_below(MSA, tree, oldroot);
		if (i > 0)
		{
			int j;
			for (j = i - 1; j >= 0; j--)
			{
				make_dirty_below(MSA, tree, mid_nodes[j]);
			}
		}
		make_dirty_below(MSA, tree, src);
	}
	make_dirty_below(MSA, tree, excess_br);
	if (parents_par != root)
	{
		make_dirty_below(MSA, tree, parents_par);
	}

	free(arr);
	free(mid_nodes);
} /* end mutate_tbr() */

/* Count total number of nodes of the tbr subtree */
int count(TREESTACK_TREE_NODES *const tree, int current)
{

	if (current == UNSET)
		return 0;
	if ((tree[current].left == UNSET) && (tree[current].right == UNSET))
	{
		return 1;
	}
	else
	{
		return count(tree, tree[current].left) + count(tree, tree[current].right);
	}
}

/* Generate an array of nodes from the tbr subtree */
int addtoarray(TREESTACK_TREE_NODES *const tree, int current, int arr[], int i)
{
	if (current == UNSET)
		return 0;
	if (tree[current].left == UNSET && tree[current].right == UNSET)
	{
		arr[i] = current;
		i++;
	}
	if (tree[current].left != UNSET)
		i = addtoarray(tree, tree[current].left, arr, i);
	if (tree[current].right != UNSET)
		i = addtoarray(tree, tree[current].right, arr, i);
	return i;
}

long lvb_reroot(Dataptr restrict MSA, TREESTACK_TREE_NODES *const BranchArray, const long oldroot, const long newroot, Lvb_bool b_with_sitestate)
/* Change the root of the tree in BranchArray from oldroot to newroot, which
 * must not be the same. Mark all internal nodes (everything but the leaves
 * and root) as "dirty". Return oldroot. */
{
	long current;				   /* current branch */
	long parnt;					   /* parent of current branch */
	long sister = UNSET;		   /* sister of current branch */
	long previous;				   /* previous branch */
	static long *oldparent = NULL; /* element i was old static parent of i */

	/* check new root is a leaf but not the current root */
	lvb_assert(newroot < MSA->n);
	lvb_assert(newroot != oldroot);
	if (oldparent == NULL)
		oldparent = (long int *)alloc(MSA->numberofpossiblebranches * sizeof(long int), "old parent alloc");

	/* create record of parents as they are now */
	for (current = 0; current < MSA->numberofpossiblebranches; current++)
		oldparent[current] = BranchArray[current].parent;

	current = newroot;
	previous = UNSET;
	while (current != oldroot)
	{
		lvb_assert(current != UNSET);
		parnt = oldparent[current]; /* original parent */
		if (current == BranchArray[parnt].left)
			sister = BranchArray[parnt].right;
		else if (current == BranchArray[parnt].right)
			sister = BranchArray[parnt].left;
		else /* error in tree structure */
			crash("internal error in function lvb_reroot(): current\n"
				  "branch %ld has old parent %ld, but old parent does not\n"
				  "have it as a child",
				  current, parnt);
		BranchArray[current].parent = previous; /* now chld of prev. */

		/* make former parent the new left child, and former sister the
		 * new right child of the current branch */
		BranchArray[current].left = parnt;
		BranchArray[current].right = sister;
		BranchArray[parnt].parent = current;
		BranchArray[sister].parent = current;

		/* move towards original root, i.e. to original parent of
		 * current branch */
		previous = current;
		current = parnt;
	}

	/* former root is now a normal leaf, without descendants */
	BranchArray[oldroot].left = UNSET;
	BranchArray[oldroot].right = UNSET;

	if (b_with_sitestate)
	{
		for (current = MSA->n; current < MSA->numberofpossiblebranches; current++)
			BranchArray[current].sitestate[0] = 0U;
	}
	return oldroot;
} /* end lvb_reroot() */

long arbreroot(Dataptr MSA, TREESTACK_TREE_NODES *const tree, const long oldroot)
/* Change tree's root arbitrarily, to a leaf other than oldroot.
 * Mark all nodes other than the leaves and root "dirty".
 * Return the number of the new root. */
{
	long newroot; /* new root */

	/* find a leaf that is not the current root */
	long ugg_minus_1 = MSA->n - 1;
	do
	{
		newroot = randpint(ugg_minus_1);
	} while (newroot == oldroot);

	lvb_reroot(MSA, tree, oldroot, newroot, LVB_TRUE);
	return newroot;

} /* end arbreroot() */

long getsister(const TREESTACK_TREE_NODES *const BranchArray, const long branch)
/* return number of sister of branch branch in tree in BranchArray, or UNSET if
 * branch has none */
{
	long parnt; /* parent of current branch */

	parnt = BranchArray[branch].parent;
	if (parnt == UNSET)
		return UNSET;
	if (branch == BranchArray[parnt].left)
		return BranchArray[parnt].right;
	else if (branch == BranchArray[parnt].right)
		return BranchArray[parnt].left;
	else /* error in tree structure */
	{
		cr_bpnc(BranchArray, branch);
		return 0; /* NEVER reached but it shuts up compilers */
	}

} /* end getsister() */

long childadd(TREESTACK_TREE_NODES *const tree, const long destination, const long newchild)
/* replace unset child of destination with newchild, and return
 * destination */
{
	if (tree[destination].right == UNSET)
		tree[destination].right = newchild;
	else if (tree[destination].left == UNSET)
		tree[destination].left = newchild;
	else /* error: destination already has 2 children */
		cr_chaf(tree, destination, newchild);

	tree[newchild].parent = destination;
	return destination;

} /* end childadd() */

void cr_chaf(const TREESTACK_TREE_NODES *const BranchArray, const long destination, const long newchild)
/* crash because we want to add branch newchild to the children of
 * branch destination in tree in BranchArray, but it already has two so
 * there is no room */
{
	crash("internal error in tree array %p: cannot make branch %ld a\n"
		  "child of branch %ld since this already has 2 children (left is\n"
		  "branch %ld, right is branch %ld)",
		  (const void *)BranchArray,
		  newchild, destination, BranchArray[destination].left,
		  BranchArray[destination].right);

} /* end cr_chaf() */

void cr_bpnc(const TREESTACK_TREE_NODES *const BranchArray, const long branch)
/* crash because branch branch in tree in BranchArray is not connected to
 * its parent, i.e. it is not a child of the branch it claims as
 * parent, according to that 'parent's' record of its own children */
{
	const long parnt = BranchArray[branch].parent; /* parent of branch */

	crash("internal error in tree array %p: branch record %ld says\n"
		  "it has parent %ld, but branch record %ld has left child %ld\n"
		  "and right child %ld",
		  (const void *)BranchArray, branch, parnt,
		  parnt, BranchArray[parnt].left, BranchArray[parnt].right);

} /* end cr_bpnc() */

TREESTACK_TREE_NODES *mvBranch(long numberofpossiblebranches, TREESTACK_TREE_NODES *const dest, const TREESTACK_TREE_NODES *const src)
{
	long i;
	Lvb_bit_length *tmp_sitestate;
	for (i = 0; i < numberofpossiblebranches; i++)
	{
		tmp_sitestate = dest[i].sitestate;
		dest[i] = src[i];
		dest[i].sitestate = tmp_sitestate;
	}
	return dest;
}

void treecopy(Dataptr restrict MSA, TREESTACK_TREE_NODES *const dest, const TREESTACK_TREE_NODES *const src, Lvb_bool b_with_sitestate)
/* copy tree from src to dest; dest must be totally distinct from source
 * in memory, and have enough space; the approach used below may fail if
 * treealloc() is changed */
{
	long i;							   /* loop counter */
	Lvb_bit_length *tmp_sitestate;	   /* temporary variable used in copy */
	unsigned char *src_statesets_all;  /* start of source's statesets */
	unsigned char *dest_statesets_all; /* start of dest's statesets */

	if (b_with_sitestate)
	{
		/* scalars */
		for (i = 0; i < MSA->numberofpossiblebranches; i++)
		{
			tmp_sitestate = dest[i].sitestate;
			dest[i] = src[i];
			dest[i].sitestate = tmp_sitestate; /* keep dest's stateset arrs for dest */
		}

		/* stateset arrays */
		src_statesets_all = ((unsigned char *)src) + MSA->numberofpossiblebranches * sizeof(TREESTACK_TREE_NODES);
		dest_statesets_all = ((unsigned char *)dest) + MSA->numberofpossiblebranches * sizeof(TREESTACK_TREE_NODES);
		memcpy(dest_statesets_all, src_statesets_all, MSA->numberofpossiblebranches * MSA->bytes);
	}
	else
	{
		/* only the scalars */
		for (i = 0; i < MSA->numberofpossiblebranches; i++)
		{
			tmp_sitestate = dest[i].sitestate;
			dest[i] = src[i];
			dest[i].sitestate = tmp_sitestate;
		}
	}

} /* end treecopy() */

void copy_sitestate(Dataptr restrict MSA, Objset *p_sitestate_1)
{

	long to_copy;
	for (long i = 0; i < MSA->nsets; i++)
	{
		to_copy = sitestate_2[i].cnt * sizeof(long);
		if (p_sitestate_1[i].set == NULL)
		{ // need to alloc memory
			p_sitestate_1[i].set = (long *)alloc(to_copy, "object set object arrays");
		}
		else if (p_sitestate_1[i].cnt != sitestate_2[i].cnt)
		{
			p_sitestate_1[i].set = (long *)realloc(p_sitestate_1[i].set, to_copy);
			if (p_sitestate_1[i].set == NULL)
			{
				crash("out of memory: cannot increase allocation for best sitestate %ld elements", to_copy);
			}
		}
		memcpy(p_sitestate_1[i].set, sitestate_2[i].set, to_copy);
		p_sitestate_1[i].cnt = sitestate_2[i].cnt;
	}
}

void PullRandomTree(Dataptr MSA, TREESTACK_TREE_NODES *const BranchArray)
/* fill BranchArray with a random tree, where BranchArray[0] is the root; all branches
 * in this random tree are marked as "dirty" */
{
	Lvb_bool *leafmask;		   /* LVB_TRUE where branch in array is a leaf */
	long *currentbranchobject; /* element i is obj associated with branch i */

	treeclear(MSA, BranchArray);
	leafmask = GenerateRandomTopology(MSA, BranchArray, MSA->n);
	currentbranchobject = randleaf(MSA, BranchArray, leafmask, MSA->n);
	tree_make_canonical(MSA, BranchArray, currentbranchobject);

} /* end PullRandomTree() */

static void wherever_was_now_say(Dataptr restrict MSA, TREESTACK_TREE_NODES *const BranchArray, long was, long now)
{
	long branchno; /* loop counter */

	for (branchno = 0; branchno < MSA->numberofpossiblebranches; branchno++)
	{
		if (BranchArray[branchno].parent == was)
			BranchArray[branchno].parent = now;
		if (BranchArray[branchno].left == was)
			BranchArray[branchno].left = now;
		if (BranchArray[branchno].right == was)
			BranchArray[branchno].right = now;
	}

} /* end wherever_was_now_say() */

void tree_make_canonical(Dataptr MSA, TREESTACK_TREE_NODES *const BranchArray, long *currentbranchobject)
/* ensure that objects 0, 1, 2, ... n-1 are associated with branches 0, 1, 2,
 * ... n-1, respectively; currentbranchobject indicates for each branch the currently
 * assigned object or UNSET for internal branches */
{
	long i;		 /* loop counter */
	long obj_no; /* current object number */
	long numberofpossiblebranches = MSA->numberofpossiblebranches;
	long n_lines = MSA->n;
	long impossible_1 = numberofpossiblebranches;						  /* an out-of-range branch index */
	long impossible_2 = numberofpossiblebranches + 1;					  /* an out-of-range branch index */
	long root = UNSET;													  /* root branch index */
	TREESTACK_TREE_NODES tmp_1, tmp_2;									  /* temporary branches for swapping */
	unsigned char *ss0_start = (unsigned char *)BranchArray[0].sitestate; /* start of state set memory */
	Lvb_bool swap_made;													  /* flag to indicate swap made */
	long tmp;															  /* for swapping */

	do
	{
		swap_made = LVB_FALSE;
		for (i = 0; i < numberofpossiblebranches; i++)
		{
			obj_no = currentbranchobject[i];
			if ((obj_no != UNSET) && (obj_no != i))
			{
				tmp_1 = BranchArray[obj_no];
				wherever_was_now_say(MSA, BranchArray, obj_no, impossible_1);
				tmp_2 = BranchArray[i];
				wherever_was_now_say(MSA, BranchArray, i, impossible_2);
				if (tmp_1.parent == i)
					tmp_1.parent = impossible_2;
				if (tmp_1.left == i)
					tmp_1.left = impossible_2;
				if (tmp_1.right == i)
					tmp_1.right = impossible_2;
				if (tmp_2.parent == obj_no)
					tmp_2.parent = impossible_1;
				if (tmp_2.left == obj_no)
					tmp_2.left = impossible_1;
				if (tmp_2.right == obj_no)
					tmp_2.right = impossible_1;
				BranchArray[i] = tmp_1;
				BranchArray[obj_no] = tmp_2;
				wherever_was_now_say(MSA, BranchArray, impossible_1, i);
				wherever_was_now_say(MSA, BranchArray, impossible_2, obj_no);
				tmp = currentbranchobject[i];
				currentbranchobject[i] = currentbranchobject[obj_no];
				currentbranchobject[obj_no] = tmp;
				swap_made = LVB_TRUE;
			}
		}
	} while (swap_made == LVB_TRUE);

	/* patch up assignment of sitestate memory to prevent trouble in treecopy() */
	for (i = 0; i < numberofpossiblebranches; i++)
	{
		BranchArray[i].sitestate = (Lvb_bit_length *)(ss0_start + i * MSA->bytes);
	}

	for (i = 0; i < n_lines; i++)
	{
		if (BranchArray[i].parent == UNSET)
		{
			lvb_assert(root == UNSET);
			root = i;
		}
	}

	if (root != 0)
	{
		lvb_reroot(MSA, BranchArray, root, 0, LVB_TRUE);
	}

	/* check things didn't go haywire */
	for (i = 0; i < n_lines; i++)
	{
		lvb_assert(currentbranchobject[i] == i);
	}
	for (i = n_lines; i < numberofpossiblebranches; i++)
	{
		lvb_assert(currentbranchobject[i] == UNSET);
	}

} /* end tree_make_canonical() */

TREESTACK_TREE_NODES *treealloc(Dataptr restrict MSA, Lvb_bool b_with_sitestate)
/* Return array of numberofpossiblebranches branches with scalars all UNSET, and all
 * statesets allocated for m characters but marked "dirty". Crash
 * verbosely if impossible. Memory is allocated once only, as a contiguous
 * block for the branch data structures followed by all their statesets.
 * So, to deallocate the tree, call the standard library function free()
 * ONCE ONLY, passing it the address of the first branch struct. If this
 * allocation approach is changed, be sure to change treecopy() too. */
{
	TREESTACK_TREE_NODES *BranchArray;			/* tree */
	unsigned char *CurrentTreeArray_uchar_star; /* tree as unsigned char */
	unsigned char *ss0_start;					/* start of first stateset */
	long i;										/* loop counter */

	lvb_assert(MSA->numberofpossiblebranches >= MIN_BRANCHES);
	lvb_assert(MSA->numberofpossiblebranches <= MAX_BRANCHES);

	if (b_with_sitestate)
		BranchArray = (TREESTACK_TREE_NODES *)alloc(MSA->tree_bytes, "tree with statesets");
	else
	{ /* don't need to do anything else */
		BranchArray = (TREESTACK_TREE_NODES *)alloc(MSA->tree_bytes_without_sitestate, "tree without statesets");
		return BranchArray;
	}

	CurrentTreeArray_uchar_star = (unsigned char *)BranchArray;
	ss0_start = CurrentTreeArray_uchar_star + MSA->numberofpossiblebranches * sizeof(TREESTACK_TREE_NODES);

	/* crash if state set memory is misaligned for uint32_t */
	lvb_assert(((intptr_t)ss0_start % NIBBLE_WIDTH) == 0);
	lvb_assert((MSA->bytes % NIBBLE_WIDTH) == 0);

	for (i = 0; i < MSA->numberofpossiblebranches; i++)
	{
		BranchArray[i].sitestate = (Lvb_bit_length *)(ss0_start + i * MSA->bytes);
		*BranchArray[i].sitestate = 0U; /* make durty */
	}

	/*make_dirty_tree(MSA, BranchArray);  */
	return BranchArray;

} /* end treealloc() */

Lvb_bool *GenerateRandomTopology(Dataptr MSA, TREESTACK_TREE_NODES *const BranchArray, const long nobjs)
/* fill BranchArray with tree of random topology, where BranchArray[0] is root;
 * return static array where element i is LVB_TRUE if BranchArray[i] is a
 * leaf or LVB_FALSE if it is not; this array will be overwritten in
 * subsequent calls */
{
	long i;								  /* loop counter */
	long leaves = 0;					  /* number of leaves */
	long nextfree = 0;					  /* next unused element of BranchArray */
	long togrow;						  /* random candidate for sprouting */
	static Lvb_bool isleaf[MAX_BRANCHES]; /* return value */

	lvb_assert(nobjs == MSA->n);

	/* clear the leaf mask */
	for (i = 0; i < MSA->numberofpossiblebranches; i++)
		isleaf[i] = LVB_FALSE;

	/* start with initial tree of 3 leaves */
	BranchArray[0].parent = UNSET;
	isleaf[nextfree++] = LVB_TRUE;
	BranchArray[0].left = nextfree;
	BranchArray[nextfree].parent = 0;
	isleaf[nextfree++] = LVB_TRUE;
	BranchArray[0].right = nextfree;
	BranchArray[nextfree].parent = 0;
	isleaf[nextfree++] = LVB_TRUE;
	leaves = 3;

	/* sprout! */
	while (leaves < nobjs)
	{
		do /* select a random leaf other than the root */
		{
			togrow = 1 + randpint(nextfree - 2);
		} while (isleaf[togrow] == LVB_FALSE);

		/* left child */
		BranchArray[togrow].left = nextfree;
		BranchArray[nextfree].parent = togrow;
		isleaf[nextfree++] = LVB_TRUE;

		/* right child */
		BranchArray[togrow].right = nextfree;
		BranchArray[nextfree].parent = togrow;
		isleaf[nextfree++] = LVB_TRUE;

		/* other updates */
		isleaf[togrow] = LVB_FALSE;
		leaves++;
	}

	return isleaf;

} /* end GenerateRandomTopology() */

long *randleaf(Dataptr MSA, TREESTACK_TREE_NODES *const BranchArray, const Lvb_bool *const leafmask, const long objs)
/* randomly assign objects numbered 0 to objs - 1 to leaves of tree in
 * BranchArray; leaves in BranchArray must be indicated by corresponding
 * LVB_TRUEs in leafmask; returns static array of object numbers, where
 * elements 0..numberofpossiblebranches give the object associated with branches 0..numberofpossiblebranches
 * respectively (UNSET for internal branches); this array will be overwritten
 * in subsequent calls */
{
	long assigned = 0;							   /* for safety: should == objs at end */
	long candidate;								   /* random object */
	long i;										   /* loop counter */
	static Lvb_bool used[MAX_N];				   /* element i LVB_TRUE if object
													* i has leaf */
	static long currentbranchobject[MAX_BRANCHES]; /* object associated with branches */

	lvb_assert(objs < MAX_N);
	lvb_assert(objs == MSA->n);

	/* clear 'used' array */
	for (i = 0; i < objs; i++)
		used[i] = LVB_FALSE;

	/* clear object nos array, defaulting to internal branch, i.e. UNSET */
	for (i = 0; i < MSA->numberofpossiblebranches; i++)
		currentbranchobject[i] = UNSET;

	/* assign an object to every leaf */
	for (i = 0; i < MSA->numberofpossiblebranches; i++)
	{
		if (leafmask[i] == LVB_TRUE) /* leaf, requires object */
		{
			do /* get a new object number */
			{
				candidate = randpint(objs - 1);
			} while (used[candidate] == LVB_TRUE);

			/* assign object to leaf */
			currentbranchobject[i] = candidate;
			used[candidate] = LVB_TRUE;
			assigned++;
		}
	}

	lvb_assert(assigned == objs);

	return currentbranchobject;

} /* end randleaf() */

void SwapTrees(TREESTACK_TREE_NODES **const tree1, long *const root1,
			   TREESTACK_TREE_NODES **const tree2, long *const root2)
/* swap trees pointed to by tree1 and tree2; also swap records of their
 * roots, as pointed to by root1 and root2 */
{
	long tmproot;				   /* temporary value-holder for swapping roots */
	TREESTACK_TREE_NODES *tmptree; /* temporary value-holder for swapping trees */

	/* swap pointers to branch arrays */
	tmptree = *tree1;
	*tree1 = *tree2;
	*tree2 = tmptree;

	/* swap roots */
	tmproot = *root1;
	*root1 = *root2;
	*root2 = tmproot;

} /* end SwapTrees() */

void treedump(Dataptr MSA, FILE *const stream, const TREESTACK_TREE_NODES *const tree, Lvb_bool b_with_sitestate)
/* send tree as table of integers to file pointed to by stream */
{
	long i; /* loop counter */
	long j; /* loop counter */

	fprintf(stream, "TREESTACK_TREE_NODES\tParent\tLeft\tRight\tChanges\tDirty\tSset_arr\tSsets\n");
	for (i = 0; i < MSA->numberofpossiblebranches; i++)
	{
		fprintf(stream, "%ld\t%ld\t%ld\t%ld\t%ld", i, tree[i].parent, tree[i].left, tree[i].right, tree[i].changes);
		if (tree[i].sitestate[0] == 0U)
			fprintf(stream, "\tyes");
		else
			fprintf(stream, "\tno");
		if (b_with_sitestate)
		{
			fprintf(stream, "\t%p", (void *)tree[i].sitestate);
			for (j = 0; j < MSA->nwords; j++)
			{
				fprintf(stream, "\t0%o", (unsigned)tree[i].sitestate[j]);
			}
		}
		else
		{
			fprintf(stream, "\tversion without sitestate");
		}
		fprintf(stream, "\n");
	}
	fprintf(stream, "\n");
	if (ferror(stream) != 0)
		cr_uxe(stream, "dumping tree");

} /* end treedump() */

void treedump_screen(Dataptr MSA, const TREESTACK_TREE_NODES *const tree)
/* send tree as table of integers to file pointed to by stream */
{
	long i; /* loop counter */

	printf("TREESTACK_TREE_NODES\tParent\tLeft\tRight\tChanges\tDirty\n");
	for (i = 0; i < MSA->numberofpossiblebranches; i++)
	{
		printf("%ld\t%ld\t%ld\t%ld\t%ld", i, tree[i].parent, tree[i].left, tree[i].right, tree[i].changes);
		if (tree[i].sitestate == NULL)
			printf("\tNULL\n");
		else if (tree[i].sitestate[0] == 0U)
			printf("\tyes\n");
		else
			printf("\tno\n");
	}
	printf("\n");

} /* end treedump() */

void cr_uxe(FILE *const stream, const char *const msg)
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

void lvb_treeprint(Dataptr MSA, FILE *const stream, const TREESTACK_TREE_NODES *const BranchArray, const long root)
/* print tree in BranchArray (of root root) in bracketed text form to stream stream,
 * in unrooted form */
{
	ur_print(MSA, stream, BranchArray, root);
} /* end lvb_treeprint() */

void ur_print(Dataptr MSA, FILE *const stream, const TREESTACK_TREE_NODES *const BranchArray, const long root)
/* send tree in BranchArray, of root root, to file pointed to by stream in
 * unrooted form */
{
	long obj;								 /* current object */
	static Lvb_bool doneabsroot = LVB_FALSE; /* have output root */
	static Lvb_bool usecomma;				 /* output clade sep. */
	char *tmp_title;						 /* temporary string */

	obj = root;

	if (doneabsroot == LVB_FALSE) /* print whole tree */
	{
		/* start tree */
		tmp_title = (char *)alloc(strlen(MSA->rowtitle[obj]) + 1, "temp. title");
		strcpy(tmp_title, MSA->rowtitle[obj]);

		while (tmp_title[strlen(tmp_title) - 1] == ' ')
		{
			tmp_title[strlen(tmp_title) - 1] = '\0';
		}
		fprintf(stream, "(%s", tmp_title);
		free(tmp_title); /* VERY LOCAL dynamic heap memory */
		usecomma = LVB_TRUE;
		doneabsroot = LVB_TRUE;

		ur_print(MSA, stream, BranchArray, BranchArray[root].left);
		ur_print(MSA, stream, BranchArray, BranchArray[root].right);
		/* end tree */
		fprintf(stream, ");\n");
		if (ferror(stream))
			crash("file error when writing unrooted tree");

		/* clean up for next call */
		usecomma = LVB_FALSE;
		doneabsroot = LVB_FALSE;
	}
	else /* print remainder of tree */
	{
		if (usecomma == LVB_TRUE)
			fprintf(stream, "%s", CLADESEP);
		if (root < MSA->n) /* leaf */
		{
			tmp_title = (char *)alloc(strlen(MSA->rowtitle[obj]) + 1, "temp. title");
			strcpy(tmp_title, MSA->rowtitle[obj]);
			while (tmp_title[strlen(tmp_title) - 1] == ' ')
			{
				tmp_title[strlen(tmp_title) - 1] = '\0';
			}
			fprintf(stream, "%s", tmp_title);
			free(tmp_title); /* VERY LOCAL dynamic heap memory */
			usecomma = LVB_TRUE;
		}
		else
		{
			fprintf(stream, "(");
			usecomma = LVB_FALSE;
			ur_print(MSA, stream, BranchArray, BranchArray[root].left);
			ur_print(MSA, stream, BranchArray, BranchArray[root].right);
			fputc(')', stream);
			usecomma = LVB_TRUE;
		}
	}

} /* end ur_print() */

long TopologyComparison(Dataptr MSA, Objset *sitestate_1, const TREESTACK_TREE_NODES *const tree_2, Lvb_bool b_First)
/* return 0 if the topology of tree_1 (of root root_1) is the same as
 * that of tree_2 (of root root_2), or non-zero if different */
{
	//	b_First = LVB_TRUE;
	if (b_First == LVB_TRUE)
	{
		makesets(MSA, tree_2, 0 /* always root zero */);
	}
	return setstcmp(MSA, sitestate_1, sitestate_2, b_First /* this one is the static */);
} /* end TopologyComparison() */

long setstcmp(Dataptr MSA, Objset *const oset_1, Objset *const oset_2, Lvb_bool b_First) /* this one is the static */
/* return 0 if the same sets of objects are in oset_1 and oset_2,
 * and non-zero otherwise */
{
	long i; /* loop counter */

	/* compare the set arrays */
	for (i = 0; i < MSA->nsets; i++)
	{
		if (oset_1[i].cnt != oset_2[i].cnt)
			return 1;
		if (memcmp(oset_1[i].set, oset_2[i].set, sizeof(long) * oset_1[i].cnt) != 0)
			return 1;
	}
	return 0;
} /* end setstcmp() */

long setstcmp_with_sitestate2(Dataptr MSA, Objset *const oset_1)
/* return 0 if the same sets of objects are in oset_1 and oset_2,
 * and non-zero otherwise */
{
	long i; /* loop counter */

	for (i = 0; i < MSA->nsets; i++)
	{
		if (oset_1[i].cnt != sitestate_2[i].cnt)
			return 1;
		if (memcmp(oset_1[i].set, sitestate_2[i].set, sizeof(long) * oset_1[i].cnt) != 0)
			return 1;
	}
	return 0;
} /* end setstcmp() */

void dump_stack_to_screen(Dataptr MSA, TREESTACK *sp)
{
	for (int i = 0; i < sp->next; i++)
	{
		printf("Stack number: %d\n", i);
		dump_objset_to_screen(MSA, sp->stack[i].p_sitestate);
	}
}

void dump_objset_to_screen(Dataptr MSA, Objset *oset_1)
{
	printf("################\n##################\n");
	for (int i = 0; i < MSA->nsets; i++)
	{
		printf("%d    %ld    ", i, oset_1[i].cnt);
		for (int x = 0; x < oset_1[i].cnt; x++)
			printf("%ld   ", oset_1[i].set[x]);
		printf("\n");
	}
	printf("\n");
}

void dump_objset_to_file(Dataptr MSA, Objset *oset_1)
{
	FILE *objset = fopen("PrintObjectset", "w");
	FILE *allobjset = fopen("PrintAllObjectset", "a+");
	for (int i = 0; i < MSA->nsets; i++)
	{
		fprintf(objset, "%d    %ld    ", i, oset_1[i].cnt);
		fprintf(allobjset, "%d    %ld    ", i, oset_1[i].cnt);
		for (int x = 0; x < oset_1[i].cnt; x++)
		{
			fprintf(objset, "%ld   ", oset_1[i].set[x]);
			fprintf(allobjset, "%ld   ", oset_1[i].set[x]);
		}
		fprintf(objset, "\n");
		fprintf(allobjset, "\n");
	}
	fclose(objset);
	fclose(allobjset);
}

void dump_objset_to_screen_sitestate_2(Dataptr MSA)
{
	dump_objset_to_screen(MSA, sitestate_2);
}

void sort_array(long *p_array, int n_left, int n_rigth)
{
	int l_hold, r_hold;
	long l_pivot = *(p_array + ((n_left + n_rigth) / 2));
	long l_temp;

	l_hold = n_left;  // i=l;
	r_hold = n_rigth; // j=r;
	do
	{
		while (p_array[l_hold] < l_pivot)
			l_hold++;
		while (p_array[r_hold] > l_pivot)
			r_hold--;

		if (l_hold <= r_hold)
		{
			l_temp = p_array[l_hold];
			p_array[l_hold] = p_array[r_hold];
			p_array[r_hold] = l_temp;
			l_hold++;
			r_hold--;
		}
	} while (l_hold < r_hold);
	if (n_left < r_hold)
		sort_array(p_array, n_left, r_hold);
	if (l_hold < n_rigth)
		sort_array(p_array, l_hold, n_rigth);
}

void Sort(Dataptr MSA, Objset *const oset_2, const long nels)
/* sort the nels object sets in oset so that each is in order, and sort oset so
 * that the sets themselves are in order of size and content */
{
	/* first sort each set member list */
	omp_set_dynamic(0); /* disable dinamic threathing */
#pragma omp parallel for num_threads(MSA->n_threads_getplen)
	for (long i = 0; i < nels; i++)
		sort_array(oset_2[i].set, 0, oset_2[i].cnt - 1);

	/* now sort the arrays of sets by size and content */
	qsort(oset_2, (size_t)nels, sizeof(Objset), osetcmp);
} /* end sort() */

int osetcmp(const void *oset1, const void *oset2)
/* comparison function for object sets (type Objset):
 * return negative if *oset1 is a smaller set of objects than *oset2 or
 * is the same size but with a list of elements that compares lower;
 * return positive if *ostet1 is bigger or the same size but with a
 * list of elements that compares higher; return 0 if they are the
 * same; N.B. the object numbers must be in numerical order within the
 * sets */
{
	long i;											 /* loop cntr */
	const Objset loset_1 = *((const Objset *)oset1); /* typed */
	const Objset loset_2 = *((const Objset *)oset2); /* typed */

	/* sets of different size differ */
	if (loset_1.cnt < loset_2.cnt)
		return -1;
	else if (loset_1.cnt > loset_2.cnt)
		return +1;

	/* if we reach here, sets are equal size, so we see if sets'
	 * contents differ */
	for (i = 0; i < loset_1.cnt; i++)
	{
		if (loset_1.set[i] < loset_2.set[i])
			return -1;
		else if (loset_1.set[i] > loset_2.set[i])
			return +1;
	}

	/* if we reach here, really the sets are the same */
	return 0;

} /* end osetcmp() */

void makesets(Dataptr MSA, const TREESTACK_TREE_NODES *const tree_2, const long root)
/* fill static sitestate_1 and static sitestate_2 with arrays of object sets for
 * tree_1 and tree_2 (of root_1 and root_2 respectively), and return
 * the extent of each array;
 * the trees must have the same object in the root branch;
 * arrays will be overwritten on subsequent calls */
{
	if (sitestate_2[0].set == NULL)
	{ /* first call, allocate memory  to the static sitestate_2*/
		ssarralloc(MSA, sitestate_2);
	}

	fillsets(MSA, sitestate_2, tree_2, root);
	Sort(MSA, sitestate_2, MSA->nsets);

} /* end makesets() */

void ssarralloc(Dataptr MSA, Objset *nobjset_2)
/* Fill nobjset[0..nsets-1] with pointers each pointing to newly
 * allocated space for setsize objects; assumes nobjset points to the
 * first element of an array with at least nsets elements.
 * Now allocates contiguous memory
 */
{
	long i;		  /* loop counter */
	long space;	  /* space for object set contents (bytes) */
	long *memory; /* memory for the object set contents */

	/* if long int arithmetic will overflow, crash instead */
	lvb_assert(log_wrapper(MSA->mssz) + log_wrapper(MSA->nsets) + log_wrapper(sizeof(long)) < log_wrapper(LONG_MAX));

	space = MSA->mssz * MSA->nsets * sizeof(long);
	memory = (long *)alloc(space, "object set arrays");
	for (i = 0; i < MSA->nsets; i++)
	{
		nobjset_2[i].set = &(memory[MSA->mssz * i]);
		nobjset_2[i].cnt = UNSET;
	}

} /* end ssarralloc() */

void fillsets(Dataptr MSA, Objset *const sstruct, const TREESTACK_TREE_NODES *const tree, const long root)
/* fill object sets in sstruct with all sets of objects in tree tree,
 * descended from but not including root and not including sets of one
 * object */
{
	static long i = UNSET; /* current set being filled */

	if (i == UNSET) /* not a recursive call */
	{
		i = 0;

		/* avoid generating sets for true root and leaves */
		if (tree[root].left >= MSA->n) /* interior */
			fillsets(MSA, sstruct, tree, tree[root].left);
		if (tree[root].right >= MSA->n) /* interior */
			fillsets(MSA, sstruct, tree, tree[root].right);

		i = UNSET; /* clean up for next non-recursive call */
		return;
	}
	if (tree[root].left != UNSET) /* not leaf */
	{
		getobjs(MSA, tree, root, sstruct[i].set, &sstruct[i].cnt);
		i++;
		fillsets(MSA, sstruct, tree, tree[root].left);
		fillsets(MSA, sstruct, tree, tree[root].right);
		return;
	}

} /* end fillsets */

void getobjs(Dataptr MSA, const TREESTACK_TREE_NODES *const BranchArray, const long root,
			 long *const objarr, long *const cnt)
/* fill objarr (which must be large enough) with numbers of all objects
 * in the tree in BranchArray in the clade starting at branch root;
 * fill the number pointed to by cnt with the number of objects found
 * (i.e. the number of elements written to objarr) */
{
	*cnt = 0;
	realgetobjs(MSA, BranchArray, root, objarr, cnt);

} /* end getobjs() */

void realgetobjs(Dataptr MSA, const TREESTACK_TREE_NODES *const BranchArray, const long root,
				 long *const objarr, long *const cnt)
/* fill objarr (which must be large enough) with numbers of all objects
 * in the tree in BranchArray in the clade starting at branch root;
 * fill the number pointed to by cnt, which must initially be zero,
 * with the number of objects found (i.e. the number of elements
 * written to objarr); this function should not be called from anywhere
 * except getobjs(), which is a safer interface */
{
	if (root < MSA->n)
	{
		objarr[*cnt] = root;
		++(*cnt);
	}
	else
	{
		if (BranchArray[root].left != UNSET)
			realgetobjs(MSA, BranchArray, BranchArray[root].left, objarr, cnt);
		if (BranchArray[root].right != UNSET)
			realgetobjs(MSA, BranchArray, BranchArray[root].right, objarr, cnt);
	}

} /* end realgetobjs() */

void ss_init(Dataptr MSA, TREESTACK_TREE_NODES *tree, Lvb_bit_length **enc_mat)
/* copy m states from enc_mat to the stateset arrays for the leaves in tree,
 * including padding at the end; the nth entry in enc_mat is assumed to be the
 * encoded state sets for object number n in the tree; non-leaf branches in the
 * tree are marked "dirty"; the root branch struct is marked "clean" since it
 * is also a terminal */
{
	long i; /* loop counter */
	for (i = 0; i < MSA->n; i++)
		memcpy(tree[i].sitestate, enc_mat[i], MSA->bytes);
	for (i = MSA->n; i < MSA->numberofpossiblebranches; i++)
		tree[i].sitestate[0] = 0U;

} /* end ss_init() */

std::string MakeHashSet(Dataptr MSA, const TREESTACK_TREE_NODES *const tree_2, const long root)
/* fill static sitestate_2 with arrays of object sets for
 * tree_2, and return sitestate_2 string
 * arrays and strings overwritten on subsequent calls */
{
	if (sitestate_2[0].set == NULL)
	{ /* first call, allocate memory  to the static sitestate_2*/
		ssarralloc(MSA, sitestate_2);
	}

	fillsets(MSA, sitestate_2, tree_2, root);
	Sort(MSA, sitestate_2, MSA->nsets);
	// dump_objset_to_file(MSA, sitestate_2);

	return ConvertSiteSetToString(MSA, sitestate_2);
} /* end MakeHashSet() */

/* convert sset to string  */
std::string ConvertSiteSetToString(Dataptr MSA, Objset *oset_1)
{
	std::ostringstream os;
	for (int i = 0; i < MSA->nsets; i++)
	{
		os << i << "    " << oset_1[i].cnt << "    ";
		for (int x = 0; x < oset_1[i].cnt; x++)
			os << oset_1[i].set[x] << "   ";
		os << std::endl;
	}

	std::string sitesetstr(os.str());
	return sitesetstr;
}
