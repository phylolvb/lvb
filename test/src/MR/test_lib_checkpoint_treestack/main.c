/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Chang Sik Kim,
Maximilian Strobl and Martyn Winn
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


/* test of checkpoint_treestack() and restore_treestack() */

#define RAND_TREES 705		/* compare treestacks of this many trees */
#define SEED 49203002		/* arbitrary */

#ifdef CHECKPOINT_INTERVAL
	#undef CHECKPOINT_INTERVAL
#endif
#define CHECKPOINT_INTERVAL 4	/* checkpoint every ... trees */


int main(int argc, char **argv)
{
    FILE *fp;				/* checkpoint file */
    long i;				/* loop counter */
    Dataptr matrix;			/* data matrix */
    DataSeqPtr matrix_seq_data;
    Params rcstruct;			/* configurable parameters */
    long root1;				/* root of tree 1 */
    long root2;				/* root of tree 2 */
    long success_cnt = 0;		/* successful tree comparisons */
    Branch *tree1;			/* first tree to compare */
    Branch *tree2;			/* second tree to compare */
    int my_id;				/* MPI process ID */
    Treestack *s_no_checkpoint;		/* tree stack without checkpointing */
    Treestack *s_with_checkpoint;	/* tree stack with checkpointing */

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (my_id == 0) {

		lvb_initialize();
		getparam(&rcstruct, argc, argv);
		matrix = (Dataptr) alloc(sizeof(DataStructure), "alloc data structure");
		matrix_seq_data = (DataSeqPtr) alloc(sizeof(DataSeqStructure), "alloc data structure");
		phylip_dna_matrin("infile", FORMAT_PHYLIP, matrix, matrix_seq_data);
		matchange(matrix, matrix_seq_data, rcstruct);
		tree1 = treealloc(matrix, LVB_TRUE);
		tree2 = treealloc(matrix, LVB_TRUE);
		s_with_checkpoint = treestack_new();
		s_no_checkpoint = treestack_new();

		/* fill a treestack without checkpointing */
		rinit(SEED);
		for (i = 0; i < RAND_TREES; i++){
			randtree(matrix, tree1);
			root1 = arbreroot(matrix, tree1, 0);
			treestack_push(matrix, s_no_checkpoint, tree1, root1, LVB_FALSE);
		}

		/* fill a treestack with frequent checkpointing */
		rinit(SEED);
		for (i = 0; i < RAND_TREES; i++)
		{
			randtree(matrix, tree2);
			root2 = arbreroot(matrix, tree2, 0);
			treestack_push(matrix, s_with_checkpoint, tree2, root2, LVB_FALSE);
			if ((i % CHECKPOINT_INTERVAL) == 0) {
				fp = fopen("treestack_checkpoint", "wb");
				checkpoint_treestack(fp, s_with_checkpoint, matrix, LVB_FALSE);
				lvb_assert(fclose(fp) == 0);
				treestack_free(s_with_checkpoint);
				fp = fopen("treestack_checkpoint", "rb");
				restore_treestack(fp, s_with_checkpoint, matrix, LVB_FALSE);
				lvb_assert(fclose(fp) == 0);
			}
			remove("treestack_checkpoint");
		}

		/* compare the two treestacks */
		for (i = 0; i < RAND_TREES; i++)
		{
			treestack_pop(matrix, tree1, &root1, s_no_checkpoint, LVB_FALSE);
			treestack_pop(matrix, tree2, &root2, s_with_checkpoint, LVB_FALSE);
			if (treecmp(matrix, tree1, tree2, root1, LVB_TRUE) == 0) success_cnt++;
		}

		if (success_cnt == RAND_TREES) {
			printf("test passed\n");
		}
		else {
			printf("FATAL ERROR\n");
		}
		free(tree1);
		free(tree2);
    }
    MPI_Finalize();
}



