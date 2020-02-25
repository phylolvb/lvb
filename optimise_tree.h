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

#include "lvb.h"

#define CLADESEP ","	/* clade separator for trees */

#ifndef NP_Implementation
/* maximum number of object sets per tree */
#define MAX_SSET_SIZE (MAX_N - 3)
#endif

static void cr_bpnc(const Branch *const barray, const long branch);
static void cr_chaf(const Branch *const barray, const long destination, const long newchild);
static void cr_uxe(FILE *const stream, const char *const msg);
static int osetcmp(const void *oset1, const void *oset2);
static void fillsets(Dataptr, Objset *const sstruct, const Branch *const tree, const long root);
static long getsister(const Branch *const barray, const long branch);
static Lvb_bool *randtopology(Dataptr restrict, Branch *const barray, const long nobjs);
static void getobjs(Dataptr, const Branch *const barray, const long root, long *const objarr, long *const cnt);
static void realgetobjs(Dataptr restrict, const Branch *const barray, const long root, long *const objarr, long *const cnt);
static long *randleaf(Dataptr restrict, Branch *const barray, const Lvb_bool *const leafmask, const long objs);
static long setstcmp(Dataptr restrict, Objset *const oset_1, Objset *const oset_2, Lvb_bool b_First);
static void ur_print(Dataptr restrict, FILE *const stream, const Branch *const barray, const long root);
static int objnocmp(const void *o1, const void *o2);


#ifndef NP_Implementation
static void makesets(Dataptr matrix, const Branch *const tree_1, const Branch *const tree_2, const long root, Lvb_bool b_First);


static void sort(Dataptr matrix, Objset *const oset_1, Objset *const oset_2, const long nels, Lvb_bool b_First);
static void ssarralloc(Dataptr restrict matrix, Objset *nobjset_1, Objset *nobjset_2);
static void tree_make_canonical(Dataptr restrict, Branch *const barray, long *objnos); 

#else
void makesets(Dataptr restrict, const Branch *const tree_2, const int root);
static void sort(Dataptr matrix, Objset *const oset_1, Objset *const oset_2, const long nels, Lvb_bool b_First);
static void ssarralloc(Dataptr restrict matrix, Objset *nobjset_2);
static void tree_make_canonical(Dataptr restrict, Branch *const barray, long *objnos);

#endif

#ifndef NP_Implementation
#ifdef MAP_REDUCE_SINGLE
	void map_pushSets(int itask, KeyValue *kv, void *ptr);
#endif

/* object sets for tree 1 in comparison */
static Objset sset_1[MAX_N - 3] = { { NULL, 0 } };
#endif

/* object sets for tree 2 in comparison */
static Objset sset_2[MAX_N - 3] = { { NULL, 0 } };

// ----------------------------------------------------------------------------