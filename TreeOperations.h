/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
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

/* ========== TreeOperations.c - Interface for TreeOperations.c ========== */

#ifndef LVB_TREEOPERATIONS_H_
#define LVB_TREEOPERATIONS_H_

#include "LVB.h"

#define CLADESEP ","	/* clade separator for trees */

#ifdef LVB_MAPREDUCE

void map_pushSets(int itask, KeyValue *kv, void *ptr);
static Objset sitestate_1[MAX_N - 3] = { { NULL, 0 } };

#endif

/* object sets for tree comparison */
static Objset sitestate_2[MAX_N - 3] = { { NULL, 0 } };

static void cr_bpnc(const TREESTACK_TREE_BRANCH *const BranchArray, const long branch);
static void cr_chaf(const TREESTACK_TREE_BRANCH *const BranchArray, const long destination, const long newchild);
static void cr_uxe(FILE *const stream, const char *const msg);
static long getsister(const TREESTACK_TREE_BRANCH *const BranchArray, const long branch);
static int osetcmp(const void *oset1, const void *oset2);
static void tree_make_canonical(Dataptr restrict, TREESTACK_TREE_BRANCH *const BranchArray, long *currentbranchobject);
static void fillsets(Dataptr, Objset *const sstruct, const TREESTACK_TREE_BRANCH *const tree, const long root);
static void getobjs(Dataptr, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, long *const objarr, long *const cnt);
static long *randleaf(Dataptr, TREESTACK_TREE_BRANCH *const BranchArray, const Lvb_bool *const leafmask, const long objs);
static void realgetobjs(Dataptr, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, long *const objarr, long *const cnt);
static Lvb_bool *GenerateRandomTopology(Dataptr, TREESTACK_TREE_BRANCH *const BranchArray, const long nobjs);
static void ur_print(Dataptr restrict, FILE *const stream, const TREESTACK_TREE_BRANCH *const BranchArray, const long root);
static long setstcmp(Dataptr restrict, Objset *const oset_1, Objset *const oset_2, Lvb_bool b_First);
static void Sort(Dataptr MSA, Objset *const oset_2, const long nels);
static void ssarralloc(Dataptr restrict MSA, Objset *nobjset_2);

void dump_objset_to_screen(Dataptr MSA, Objset *oset_1);
void dump_objset_to_file(Dataptr MSA, Objset *oset_1);
void dump_objset_to_screen_sitestate_2(Dataptr MSA);

#endif  /* LVB_TREEOPERATIONS_H_ */
