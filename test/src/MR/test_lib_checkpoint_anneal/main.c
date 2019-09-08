/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
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


#include "Lvb.h"
#include "Store_states.h"

/* test of checkpoint_uni() and restore_uni() */

#define SEED 40291989		/* arbitrary */

int main(int argc, char **argv)
{
    FILE *fp;			/* checkpoint file */
    int my_id;			/* MPI process ID */
    Dataptr matrix;					/* data matrix */
    DataSeqPtr matrix_seq_data;
    Params rcstruct;

    Branch *p_current_tree, *p_current_tree_2;	/* first tree to compare */
    Branch *p_proposed_tree, *p_proposed_tree_2;	/* first tree to compare */
    long accepted = 54, accepted_2;		/* changes accepted */
	Lvb_bool dect = LVB_FALSE, dect_2;		/* should decrease temperature */
	double deltah = 0.23423, deltah_2;		/* change in energy (1 - C.I.) */
	long deltalen = 24, deltalen_2;		/* change in length with new tree */
	long failedcnt = 0, failedcnt_2; 	/* "failed count" for temperatures */
	long iter = 0, iter_2;		/* iteration of mutate/evaluate loop */
	long current_iter = 23423, current_iter_2;
	long len = 34, len_2;			/* length of current tree */
	long lenbest = 55, lenbest_2;		/* bet length found so far */
	long lendash = 53, lendash_2;		/* length of proposed new tree */
	double ln_t = 6.78, ln_t_2;		/* ln(current temperature) */
	long t_n = 543, t_n_2;		/* ordinal number of current temperature */
	double t0 = 0.002345, t0_2;		/* ordinal number of current temperature */
	double pacc = 24.1, pacc_2;		/* prob. of accepting new config. */
	long proposed = 345, proposed_2;		/* trees proposed */
	double r_lenmin, r_lenmin_2;		/* minimum length for any tree */
	long rootdash = 0, rootdash_2;		/* root of new configuration */
	double t = 0.454545, t_2;				/* current temperature */
	double grad_geom = 0.99, grad_geom_2;			/* "gradient" of the geometric schedule */
	double grad_linear = 3.64 * LVB_EPS, grad_linear_2; 	/* gradient of the linear schedule */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (my_id == 0) {
    	rinit(SEED);
    	matrix = (Dataptr) alloc(sizeof(DataStructure), "alloc data structure");
    	matrix_seq_data = (DataSeqPtr) alloc(sizeof(DataSeqStructure), "alloc data structure");
    	getparam(&rcstruct, argc, argv);
    	phylip_dna_matrin("infile", FORMAT_PHYLIP, matrix, matrix_seq_data);
    	matchange(matrix, matrix_seq_data, rcstruct);

    	fp = fopen("uni_anneal", "wb");
    	/* fill a treestack without checkpointing */
    	Lvb_bool b_with_sset_current_tree = LVB_TRUE;
    	Lvb_bool b_with_sset_proposed_tree = LVB_TRUE;
    	p_current_tree = treealloc(matrix, b_with_sset_current_tree);
    	p_proposed_tree = treealloc(matrix, b_with_sset_proposed_tree);
    	p_current_tree_2 = treealloc(matrix, b_with_sset_current_tree);
    	p_proposed_tree_2 = treealloc(matrix, b_with_sset_proposed_tree);
    	rinit(SEED);
    	randtree(matrix, p_current_tree);// rootdash = arbreroot(matrix, p_current_tree, rootdash);
    	randtree(matrix, p_proposed_tree);// rootdash = arbreroot(matrix, p_proposed_tree, rootdash);
    	r_lenmin = (double) matrix->min_len_tree;
    	checkpoint_anneal(fp, matrix, accepted, dect, deltah, deltalen, failedcnt, iter, current_iter, len, lenbest,
    						lendash, ln_t, t_n, t0, pacc, proposed, r_lenmin, rootdash, t, grad_geom, grad_linear,
    						p_current_tree, b_with_sset_current_tree, p_proposed_tree, b_with_sset_proposed_tree);
    	lvb_assert(fclose(fp) == 0);
		fp = fopen("uni_anneal", "rb");
    	restore_anneal(fp, matrix, &accepted_2, &dect_2, &deltah_2, &deltalen_2, &failedcnt_2, &iter_2, &current_iter_2, &len_2, &lenbest_2,
    			    &lendash_2, &ln_t_2, &t_n_2, &t0_2, &pacc_2, &proposed_2, &r_lenmin_2, &rootdash_2, &t_2, &grad_geom_2,
    				&grad_linear_2, p_current_tree_2, b_with_sset_current_tree, p_proposed_tree_2, b_with_sset_proposed_tree);
    	lvb_assert(fclose(fp) == 0);
    	remove("uni_anneal");

    	if (accepted_2 == accepted && dect == dect_2 && deltah == deltah_2 && deltalen == deltalen_2 &&
				failedcnt == failedcnt_2 && iter == iter_2 && current_iter == current_iter_2 && len == len_2 &&
				lenbest == lenbest_2 && lendash == lendash_2 && ln_t == ln_t_2 &&
				t_n == t_n_2 && t0 == t0_2 && pacc == pacc_2 && proposed == proposed_2 && r_lenmin == r_lenmin_2 &&
				rootdash == rootdash_2 && t == t_2 && grad_geom == grad_geom_2 &&
				grad_linear == grad_linear_2) {

    		treedump_screen(matrix, p_current_tree);
    		treedump_screen(matrix, p_current_tree_2);
    		printf("test passed\n");
    	}
    	else {
    		printf("FATAL ERROR\n");
    	}
    	free(p_current_tree);
    	free(p_proposed_tree);
    	free(p_current_tree_2);
    	free(p_proposed_tree_2);
    	free(matrix);
    	free(matrix_seq_data);
    }
    MPI_Finalize();
}
