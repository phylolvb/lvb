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


#include "LVB.h"


#define RAND_TREES 705		/* compare treestacks of this many trees */
#define SEED 40291989		/* arbitrary */

/* test of checkPoint() */
int main(int argc, char **argv)
{
    FILE *fp;						/* checkpoint file */
    long i;							/* loop counter */
    int my_id, is_process_finished, n_number_blocks;
    char filename[] = "test_file";
    Parameters rcstruct, rcstruct_2;	/* configurable parameters */
    Dataptr MSA;					/* data MSA */
    Dataptr MSA;
    long root1;						/* root of tree 1 */
    TREESTACK_TREE_BRANCH *tree1;					/* first tree to compare */
    TREESTACK *tree_checkpoint;		/* tree stack without checkpointing */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (my_id == 0) {
    	rinit(SEED);

    	/* create structure, it can be in different order, only the two first int need to be there */
    	/*		1: int is_process_finished;		*/
    	/*		2: int number of blocks;		*/
    	/*		3: uni struture					*/
    	/*		4: Parameters struture				*/
    	/*		5: tree stack struture			*/

    	getparam(&rcstruct, argc, argv);
    	MSA = (Dataptr) alloc(sizeof(DataStructure), "alloc data structure");
    	MSA = (Dataptr) alloc(sizeof(DataStructure), "alloc data structure");
    	phylip_dna_matrin("infile", FORMAT_PHYLIP, MSA);
    	matchange(MSA, rcstruct);
    	tree_checkpoint = CreateNewTreestack();

		/* fill a treestack without checkpointing */
    	tree1 = treealloc(MSA, LVB_TRUE);
		rinit(SEED);
		for (i = 0; i < RAND_TREES; i++){
			randtree(MSA, tree1);
			root1 = arbreroot(MSA, tree1, 0);
			CompareTreeToTreestack(MSA, tree_checkpoint, tree1, root1, LVB_FALSE);
		}

    	is_process_finished = CHECK_POINT_PROCESS_FINISHED;

    	/* start write */
    	fp = fopen(filename, "wb");
    	lvb_assert(fp != NULL);
    	fwrite(&is_process_finished, sizeof(is_process_finished), 1, fp);
    	n_number_blocks = 3;
    	fwrite(&n_number_blocks, sizeof(n_number_blocks), 1, fp);
    	uni();		/* call one time */
    	checkpoint_uni(fp);
    	checkpoint_params(fp, &rcstruct);
    	checkpoint_treestack(fp, tree_checkpoint, MSA, LVB_FALSE);
		FreeTreestackMemory(tree_checkpoint);
    	lvb_assert(fclose(fp) == 0);
    	lvb_assert(test_consistency_state_file(filename, 1) == LVB_TRUE);
    	lvb_assert(is_process_ended(filename) == LVB_TRUE);

    	/* point the file to the being of the position */
		fp = fopen(filename, "rb");
    	lvb_assert(point_file_pointer_to_block(fp, STATE_BLOCK_PARAMETERS /* block parameters*/) == LVB_TRUE);

    	/* read parameters block and compare with the one that was saved before */
    	restore_params(fp, &rcstruct_2);
    	lvb_assert(compare_params(&rcstruct, &rcstruct_2, LVB_TRUE) == LVB_TRUE);

    	if (rcstruct.seed == 0) rcstruct.seed = 1;
    	else rcstruct.seed >>= 1;
    	lvb_assert(compare_params(&rcstruct, &rcstruct_2, LVB_TRUE) == LVB_FALSE);

    	/* restore uni data */
    	lvb_assert(point_file_pointer_to_block(fp, STATE_BLOCK_UNI /* block parameters*/) == LVB_TRUE);
    	restore_uni(fp);

    	/* a blcok ID that doesn't exist */
    	lvb_assert(point_file_pointer_to_block(fp, STATE_BLOCK_ONLY_TESTS_FAIL /* block parameters*/) == LVB_FALSE);
    	lvb_assert(fclose(fp) == 0);

    	free(tree1);
    	free(MSA);
    	free(MSA);
    	remove(filename);
    	printf("test passed\n");
    }
    MPI_Finalize();
}
