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
#include "Hash.h"

#define CLADESEP ","	/* clade separator for trees */

/* object sets for tree comparison */
Objset sitestate_2[MAX_N - 3] = { { NULL, 0 } };

void cr_bpnc(const TREESTACK_TREE_NODES *const BranchArray, const long branch);
void cr_chaf(const TREESTACK_TREE_NODES *const BranchArray, const long destination, const long newchild);
void cr_uxe(FILE *const stream, const char *const msg);
long getsister(const TREESTACK_TREE_NODES *const BranchArray, const long branch);
int osetcmp(const void *oset1, const void *oset2);
void tree_make_canonical(Dataptr restrict, TREESTACK_TREE_NODES *const BranchArray, long *currentbranchobject);
void fillsets(Dataptr, Objset *const sstruct, const TREESTACK_TREE_NODES *const tree, const long root);
void getobjs(Dataptr, const TREESTACK_TREE_NODES *const BranchArray, const long root, long *const objarr, long *const cnt);
long *randleaf(Dataptr, TREESTACK_TREE_NODES *const BranchArray, const Lvb_bool *const leafmask, const long objs);
void realgetobjs(Dataptr, const TREESTACK_TREE_NODES *const BranchArray, const long root, long *const objarr, long *const cnt);
Lvb_bool *GenerateRandomTopology(Dataptr, TREESTACK_TREE_NODES *const BranchArray, const long nobjs);
void ur_print(Dataptr restrict, FILE *const stream, const TREESTACK_TREE_NODES *const BranchArray, const long root);
long setstcmp(Dataptr restrict, Objset *const oset_1, Objset *const oset_2, Lvb_bool b_First);
void Sort(Dataptr MSA, Objset *const oset_2, const long nels);
void ssarralloc(Dataptr restrict MSA, Objset *nobjset_2);

void dump_objset_to_screen(Dataptr MSA, Objset *oset_1);
void dump_objset_to_file(Dataptr MSA, Objset *oset_1);
void dump_objset_to_screen_sitestate_2(Dataptr MSA);

#endif  /* LVB_TREEOPERATIONS_H_ */
