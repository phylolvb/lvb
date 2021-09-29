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

/* ========== Solve.h - interface for Solve.c ========== */

#ifndef LVB_SOLVE_H
#define LVB_SOLVE_H

#include "Clock.h"
#include "DataOperations.h"
#include "Hash.h"
#include "LVB.h"
#include "Print.h"
#include "Verbose.h"

static TREESTACK treestack;	/* overall best tree stack */
static TREESTACK stack_treevo;

#ifdef LVB_MAPREDUCE
long GetSoln(Dataptr restrict MSA, Parameters rcstruct, long *iter_p, Lvb_bool log_progress,
				MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer);
long CompareMapReduceTreesAnneal(Dataptr, TREESTACK *, const TREESTACK_TREE_NODES *const, long, 
                        int *, int, MISC *, MapReduce *, MapReduce *);
long CompareMapReduceTreesGetSoln(Dataptr, TREESTACK *, const TREESTACK_TREE_NODES *const, long, 
                        int *, MISC *, MapReduce *, MapReduce *);
#else
long GetSoln(Dataptr restrict, Parameters, long *, Lvb_bool);

#endif

#endif