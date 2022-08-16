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

/* ========== Solve.c - solving functions ========== */

#include "Solve.h"
#include "Para_lib.h"
#include "mpi.h"
#include "Clock.h"

using namespace std;
using namespace std::chrono;

static void lenlog(FILE *lengthfp, TREESTACK *treestack_ptr, long iteration, long length, double temperature)
/* write a message to file pointer lengthfp; iteration gives current iteration;
 * crash verbosely on write error */
{
    fprintf(lengthfp, "  %-15.8f%-15ld%-16ld%-15ld\n", temperature, iteration, treestack_ptr->next, length);
    if (ferror(lengthfp)){
    	crash("file error when logging search progress");
    }

} /* end lenlog() */

#ifdef LVB_MAPREDUCE  
long deterministic_hillclimb(Dataptr MSA, TREESTACK *treestack_ptr, const TREESTACK_TREE_NODES *const inittree,
	Parameters rcstruct, long root, FILE * const lenfp, long *current_iter, Lvb_bool log_progress, 
	MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer)
#else 
long deterministic_hillclimb(Dataptr MSA, TREESTACK *treestack_ptr, const TREESTACK_TREE_NODES *const inittree,
		Parameters rcstruct, long root, FILE * const lenfp, long *current_iter, Lvb_bool log_progress)
#endif
/* perform a deterministic hill-climbing optimization on the tree in inittree,
 * using NNI on all internal branches until no changes are accepted; return the
 * length of the best tree found; current_iter should give the iteration number
 * at the start of this call and will be used in any statistics sent to lenfp,
 * and will be updated on return */
{
    long i;							/* loop counter */
    long j;							/* loop counter */
    long number_of_internal_branches = 0;				/* count of internal branches */
    long current_tree_length;						/* current length */
    long proposed_tree_length;					/* length of proposed new config */
    long proposed_tree_root = root;			/* root of proposed new config */
    long tree_length_change;					/* change in length */
    Lvb_bool newtree;				/* accepted a new configuration */
    TREESTACK_TREE_NODES *p_current_tree;			/* current configuration */
    TREESTACK_TREE_NODES *p_proposed_tree;		/* proposed new configuration */
    unsigned int *branch_numbers_arr;				/* array of internal branch numbers */
    Lvb_bool leftright[] = { LVB_FALSE, LVB_TRUE }; /* to loop through left and right */
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */

	#ifdef LVB_MAPREDUCE  
		int *total_count;
	    int check_cmp;
	#endif

    /* "local" dynamic heap memory */
    p_current_tree = treealloc(MSA, LVB_TRUE);
    p_proposed_tree = treealloc(MSA, LVB_TRUE);
    branch_numbers_arr = (unsigned int *) alloc(MSA->numberofpossiblebranches * sizeof(unsigned int), "old parent alloc");

    treecopy(MSA, p_current_tree, inittree, LVB_TRUE);      /* current configuration */
	alloc_memory_to_getplen(MSA, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	current_tree_length = getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);

    /* identify internal branches */
    for (i = MSA->n; i < MSA->numberofpossiblebranches; i++) branch_numbers_arr[number_of_internal_branches++] = i;
    lvb_assert(number_of_internal_branches == MSA->numberofpossiblebranches - MSA->n);

    do {
		newtree = LVB_FALSE;
		for (i = 0; i < number_of_internal_branches; i++) {
			for (j = 0; j < 2; j++) {
				mutate_deterministic(MSA, p_proposed_tree, p_current_tree, root, branch_numbers_arr[i], leftright[j]);
				proposed_tree_length = getplen(MSA, p_proposed_tree, rcstruct, proposed_tree_root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
				lvb_assert (proposed_tree_length >= 1L);
				tree_length_change = proposed_tree_length - current_tree_length;

				#ifdef LVB_MAPREDUCE  
				MPI_Bcast(&tree_length_change, 1, MPI_LONG, 0,    MPI_COMM_WORLD);
				MPI_Bcast(&proposed_tree_length,  1, MPI_LONG, 0,    MPI_COMM_WORLD);

					if (tree_length_change <= 0) {
						if (tree_length_change < 0)  /* very best so far */
						{
							ClearTreestack(treestack_ptr);
							PushCurrentTreeToStack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE);
							misc->ID = treestack_ptr->next;
							misc->SB = 1;
							tree_setpush(MSA, p_proposed_tree, proposed_tree_root, mrTreeStack, misc);
							current_tree_length = proposed_tree_length;
						} else {

						  misc->SB = 0;
						  tree_setpush(MSA, p_proposed_tree, proposed_tree_root, mrBuffer, misc);
						  mrBuffer->add(mrTreeStack);
						  mrBuffer->collate(NULL);

						  misc->count = (int *) alloc( (treestack_ptr->next+1) * sizeof(int), "int array for tree comp using MR");
						  total_count = (int *) alloc( (treestack_ptr->next+1) * sizeof(int), "int array for tree comp using MR");
						  for(int i=0; i<=treestack_ptr->next; i++) misc->count[i] = 0;
						  mrBuffer->reduce(reduce_count, misc);

						  for(int i=0; i<=treestack_ptr->next; i++) total_count[i] = 0;
						  MPI_Reduce( misc->count, total_count, treestack_ptr->next+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

						  check_cmp = 1;
						  if (misc->rank == 0) { /* sum to root process */
							for(int i=1; i<=treestack_ptr->next; i++) {
								if (misc->nsets == total_count[i]) {
									check_cmp = 0; /* different */
									break;
								}
							}
						  }

						  MPI_Barrier(MPI_COMM_WORLD);
						  MPI_Bcast(&check_cmp, 1, MPI_INT, 0,    MPI_COMM_WORLD);
						  if (check_cmp == 0) {
								misc->SB = 1;
								tree_setpush(MSA, p_proposed_tree, proposed_tree_root, mrBuffer, misc);
								mrTreeStack->add(mrBuffer);
							PushCurrentTreeToStack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE);
								misc->ID = treestack_ptr->next;

								newtree = LVB_TRUE;
								SwapTrees(&p_current_tree, &root, &p_proposed_tree, &proposed_tree_root);
							}
						  free(misc->count);
						  free(total_count);
						}

					#else 
					if (tree_length_change <= 0) {
					if (tree_length_change < 0)  /* very best so far */
					{
						ClearTreestack(treestack_ptr);
						current_tree_length = proposed_tree_length;
					}
						#ifdef LVB_HASH
							if (CompareHashTreeToHashstack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1) 
						#else
							if (CompareTreeToTreestack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1) 
						#endif
					{
						newtree = LVB_TRUE;
						SwapTrees(&p_current_tree, &root, &p_proposed_tree, &proposed_tree_root);
					}
					
					#endif
				}
				if ((log_progress == LVB_TRUE) && ((*current_iter % STAT_LOG_INTERVAL) == 0)) {
					lenlog(lenfp, treestack_ptr, *current_iter, current_tree_length, 0);
				}
				*current_iter += 1;
			}
		}
    } while (newtree == LVB_TRUE);

    /* free "local" dynamic heap memory */
    free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
    free(p_current_tree);
    free(p_proposed_tree);
    free(branch_numbers_arr);

    return current_tree_length;
} /* end deterministic_hillclimb */

#ifdef LVB_MAPREDUCE  
long Anneal(Dataptr MSA, TREESTACK* treestack_ptr, TREESTACK* treevo, const TREESTACK_TREE_NODES* const inittree, Parameters rcstruct,
	long root, const double t0, const long maxaccept, const long maxpropose,
	const long maxfail, FILE* const lenfp, long* current_iter,
	Lvb_bool log_progress, MISC* misc, MapReduce* mrTreeStack, MapReduce* mrBuffer)

#else
long Anneal(Dataptr MSA, TREESTACK* treestack_ptr, TREESTACK* treevo, const TREESTACK_TREE_NODES* const inittree, Parameters rcstruct,
	long root, const double t0, const long maxaccept, const long maxpropose,
	const long maxfail, FILE* const lenfp, long* current_iter,
	Lvb_bool log_progress)

#endif
	/* seek parsimonious tree from initial tree in inittree (of root root)
	 * with initial temperature t0, and subsequent temperatures obtained by
	 * multiplying the current temperature by (t1 / t0) ** n * t0 where n is
	 * the ordinal number of this temperature, after at least maxaccept changes
	 * have been accepted or maxpropose changes have been proposed, whichever is
	 * sooner;
	 * return the length of the best tree(s) found after maxfail consecutive
	 * temperatures have led to no new accepted solution;
	 * lenfp is for output of current tree length and associated details;
	 * *current_iter should give the iteration number at the start of this call and
	 * will be used in any statistics sent to lenfp, and will be updated on
	 * return */
{
	long accepted = 0;		/* changes accepted */
	Lvb_bool dect;		/* should decrease temperature */
	double deltah;		/* change in energy (1 - C.I.) */
	long tree_length_change;		/* change in length with new tree */
	long failedcnt = 0; 	/* "failed count" for temperatures */
	long iter = 0;		/* iteration of mutate/evaluate loop */
	long current_tree_length;			/* length of current tree */
	long best_tree_length;		/* best length found so far */
	long proposed_tree_length;		/* length of proposed new tree */
	double ln_t;		/* ln(current temperature) */
	long t_n = 0;		/* ordinal number of current temperature */
	double pacc;		/* prob. of accepting new config. */
	long proposed = 0;		/* trees proposed */
	double tree_minimum_length;		/* minimum length for any tree */
	long proposed_tree_root;		/* root of new configuration */
	double t = t0;		/* current temperature */
	double grad_geom = 0.99;		/* "gradient" of the geometric schedule */
	double grad_linear = 10 * LVB_EPS; /* gradient of the linear schedule */
	/* double grad_linear = 10 * LVB_EPS; */ /* gradient of the linear schedule */
	TREESTACK_TREE_NODES* p_current_tree;			/* current configuration */
	TREESTACK_TREE_NODES* p_proposed_tree;		/* proposed new configuration */
	long* p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
	long* p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
	int* p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */
#ifdef LVB_MAPREDUCE  
	int* total_count;
	int check_cmp;
#endif

	/* variables that could calculate immediately */
	const double log_wrapper_LVB_EPS = log_wrapper(LVB_EPS);
	const double log_wrapper_grad_geom = log_wrapper(grad_geom);
	const double log_wrapper_t0 = log_wrapper(t0);
	/* REND variables that could calculate immediately */

	long w_changes_prop = 0;
	long w_changes_acc = 0;
	p_proposed_tree = treealloc(MSA, LVB_TRUE);
	p_current_tree = treealloc(MSA, LVB_TRUE);

	treecopy(MSA, p_current_tree, inittree, LVB_TRUE);	/* current configuration */

	alloc_memory_to_getplen(MSA, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	current_tree_length = getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
	dect = LVB_FALSE;		/* made LVB_TRUE as necessary at end of loop */

	lvb_assert(((float)t >= (float)LVB_EPS) && (t <= 1.0) && (grad_geom >= LVB_EPS) && (grad_linear >= LVB_EPS));

	best_tree_length = current_tree_length;

#ifdef LVB_MAPREDUCE  
	MPI_Bcast(&best_tree_length, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	PushCurrentTreeToStack(MSA, treestack_ptr, inittree, root, LVB_FALSE);
	misc->ID = treestack_ptr->next;
	misc->SB = 1;
	tree_setpush(MSA, inittree, root, mrTreeStack, misc);

#elif LVB_HASH
	CompareHashTreeToHashstack(MSA, treestack_ptr, inittree, root, LVB_FALSE);	/* init. tree initially best */
#else
	CompareTreeToTreestack(MSA, treestack_ptr, inittree, root, LVB_FALSE);	/* init. tree initially best */
#endif

	double trops_counter[3] = { 1,1,1 };
	double trops_probs[3] = { 0,0,0 };
	long trops_total = trops_counter[0] + trops_counter[1] + trops_counter[2];
	long trops_id = 0;




	/*Writing output to table.tsv*/
	FILE* pFile;
	char change[10] = "";
	if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {
		if (rcstruct.verbose == LVB_TRUE)
		{
			pFile = fopen("changeAccepted.tsv", "w");
			fprintf(pFile, "Iteration\tAlgorithm\tAccepted\tLength\tTemperature\tTreestack Size\n");
		}
	}
	tree_minimum_length = (double)MSA->min_len_tree;


#ifndef parallel
	//Control parameter for parallel_selection==1(cluster SA)
	long p = 1000;//for cluster SA
	double r = 0.05;
	long iter_cluster = 0;

	if (rcstruct.parallel_selection != 1)
		p = -1;//Clustering loop for generating best initial treee won't launch for other parallel selections

	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);



	int pad_size = sizeof(Lvb_bit_length*);
	int nItems = 2;
	int          	blocklengths_2[2] = { 4, pad_size };//content of sitestate doesn't matter
	MPI_Datatype 	types_2[2] = { MPI_LONG, MPI_BYTE };
	MPI_Datatype 	MPI_BRANCH;
	MPI_Aint     	displacements_2[2];
	displacements_2[0] = offsetof(TREESTACK_TREE_NODES, parent);
	displacements_2[1] = offsetof(TREESTACK_TREE_NODES, sitestate);
	MPI_Type_create_struct(nItems, blocklengths_2, displacements_2, types_2, &MPI_BRANCH);
	MPI_Type_commit(&MPI_BRANCH);

	if (log_progress == LVB_TRUE&&rcstruct.parallel_selection==1)
	{
		printf("\n--------------------------------------------------------");
		printf("\nCluster SA launches with %d processes, p=%ld, r=%lf, propertis of best partial tree will show every p iterations", nprocs, p, r);
		printf("\n  p:  Rearrangement:   Length:  From rank:\n");
		printf("--------------------------------------------------------\n");
	}


#endif	
	while (1) {
#ifndef parallel
		

		if (p >= 0)
		{
			iter_cluster++;

			if (iter_cluster == (p+1))
			{
				

				int from_rank;
				Bcast_best_partial_tree(MSA, current_tree_length, rank, nprocs, p_current_tree, &root, &from_rank, MPI_BRANCH);//Bcast best tree with root

				current_tree_length = getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);

				if(rank==0)
					printf("  %-15ld%-15ld%-15ld%-10d\n", p, *current_iter, current_tree_length, from_rank);

				iter_cluster = 1;
				p -= r * (double)iter;
			}

			if (p < 0)//reset annealing by best initial tree
			{
				accepted = 0;		/* changes accepted */
				failedcnt = 0; 	/* "failed count" for temperatures */
				iter = 0;		/* iteration of mutate/evaluate loop */
				*current_iter = 0;
				current_tree_length = getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);			/* length of current tree */
				t_n = 0;		/* ordinal number of current temperature */
				proposed = 0;		/* trees proposed */
				t = StartingTemperature(MSA, p_current_tree, rcstruct, root, LVB_FALSE);		/* current temperature, no need to printf */

				w_changes_prop = 0;
				w_changes_acc = 0;
				for (int i = 0; i < 3; i++)
				{
					trops_counter[i] = 1;
					trops_probs[i] = 0;
				}
				trops_total = trops_counter[0] + trops_counter[1] + trops_counter[2];
				trops_id = 0;

				if (rank == 0)
					printf("\n\n Generate the best initial tree, current_tree_length=%ld, root=%ld, new start temperature: %f\n\n", current_tree_length, root, t);

				//Canonical  LVB starts, but only Rank 0 print its progress
				if ((log_progress == LVB_TRUE) && (*current_iter == 0)) 
				{

					printf("\n--------------------------------------------------------");
					fprintf(lenfp, "\n  Temperature:   Rearrangement: TreeStack size: Length:\n");
					printf("--------------------------------------------------------\n");
				}

			}
			//ClearTreestack(treestack_ptr);
			//CompareTreeToTreestack(MSA, treestack_ptr, p_current_tree, root, LVB_FALSE);
		}


		int changeAcc = 0;
		*current_iter += 1;
		/* occasionally re-root, to prevent influence from root position */
		if ((*current_iter % REROOT_INTERVAL) == 0) {
			root = arbreroot(MSA, p_current_tree, root);

			if (p < 0)
			if ((log_progress == LVB_TRUE) && ((*current_iter % STAT_LOG_INTERVAL) == 0)) {
				lenlog(lenfp, treestack_ptr, *current_iter, current_tree_length, t);
			}
		}
#endif
		lvb_assert(t > DBL_EPSILON);

		/* mutation: alternate between the two mutation functions */
		proposed_tree_root = root;
		if (rcstruct.algorithm_selection == 2)
		{
			trops_total = trops_counter[0] + trops_counter[1] + trops_counter[2];
			trops_probs[0] = trops_counter[0] / trops_total;
			trops_probs[1] = trops_counter[1] / trops_total;
			trops_probs[2] = trops_counter[2] / trops_total;
		}

		if (rcstruct.algorithm_selection >= 1)
		{
			double random_val = uni();
			if (random_val < trops_probs[0]) {
				mutate_nni(MSA, p_proposed_tree, p_current_tree, root);	/* local change */
				strcpy(change, "NNI");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 0;
			}
			else if (random_val < trops_probs[0] + trops_probs[1]) {
				mutate_spr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
				strcpy(change, "SPR");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 1;
			}
			else {
				mutate_tbr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
				strcpy(change, "TBR");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 2;
			}
		}
		else
		{
			if (iter & 0x01) {
				mutate_spr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
				strcpy(change, "SPR");
			}
			else {
				mutate_nni(MSA, p_proposed_tree, p_current_tree, root);	/* local change */
				strcpy(change, "NNI");
			}
		}

		proposed_tree_length = getplen(MSA, p_proposed_tree, rcstruct, proposed_tree_root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
		lvb_assert(proposed_tree_length >= 1L);
		tree_length_change = proposed_tree_length - current_tree_length;
		deltah = (tree_minimum_length / (double)current_tree_length) - (tree_minimum_length / (double)proposed_tree_length);
		if (deltah > 1.0) deltah = 1.0; /* MinimumTreeLength() problem with ambiguous sites */

#ifdef LVB_MAPREDUCE  
		MPI_Bcast(&tree_length_change, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&deltah, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&proposed_tree_length, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

		if (tree_length_change <= 0)	/* accept the change */
		{
#ifdef LVB_MAPREDUCE  
			if (proposed_tree_length <= best_tree_length)	/* store tree if new */
			{
				if (proposed_tree_length < best_tree_length) {
					ClearTreestack(treestack_ptr);
					mrTreeStack->map(mrTreeStack, map_clean, NULL);
				}
				auto start_timer = high_resolution_clock::now();

				if (CompareMapReduceTrees(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, total_count,
					check_cmp, misc, mrTreeStack, mrBuffer) == 1) {
					accepted++;
					MPI_Bcast(&accepted, 1, MPI_LONG, 0, MPI_COMM_WORLD);
				}
				auto stop = high_resolution_clock::now();

				auto duration = duration_cast<microseconds>(stop - start_timer);

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE* timefunction = fopen("FunctionTimesMR", "a+");
					fprintf(timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
				}
			}
#else
			if (proposed_tree_length <= best_tree_length)	/* store tree if new */
			{
				if (proposed_tree_length < best_tree_length) {
					ClearTreestack(treestack_ptr);	/* discard old bests */
				}
#ifdef LVB_HASH
				auto start_timer = high_resolution_clock::now();
				if (CompareHashTreeToHashstack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1)
#else
				auto start_timer = high_resolution_clock::now();
				if (CompareTreeToTreestack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1)
#endif
				{
					accepted++;
				}
				auto stop = high_resolution_clock::now();

				auto duration = duration_cast<microseconds>(stop - start_timer);

#ifdef LVB_HASH

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE* timefunction = fopen("FunctionTimesHASH", "a+");
					fprintf(timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
				}
#else

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE* timefunction = fopen("FunctionTimesNP", "a+");
					fprintf(timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
				}
#endif
			}
#endif
			/* update current tree and its stats */
			current_tree_length = proposed_tree_length;
			SwapTrees(&p_current_tree, &root, &p_proposed_tree, &proposed_tree_root);

			/* very best so far */
			if (proposed_tree_length < best_tree_length) {
				best_tree_length = proposed_tree_length;
#ifdef LVB_MAPREDUCE  
				MPI_Bcast(&best_tree_length, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
			}
			if (rcstruct.algorithm_selection == 1)
				changeAcc = 1;

		}
		else	/* poss. accept change for the worse */
		{
			if (rcstruct.algorithm_selection == 2)
				w_changes_prop++;
			/* Mathematically,
			 *     Pacc = e ** (-1/T * deltaH)
			 *     therefore ln Pacc = -1/T * deltaH
			 *
			 * Computationally, if Pacc is going to be small, we
			 * can assume Pacc is 0 without actually working it
			 * out.
			 * i.e.,
			 *     if ln Pacc < ln eps, let Pacc = 0
			 * substituting,
			 *     if -deltaH / T < ln eps, let Pacc = 0
			 * rearranging,
			 *     if -deltaH < T * ln eps, let Pacc = 0
			 * This lets us work out whether Pacc will be very
			 * close to zero without dividing anything by T. This
			 * should prevent overflow. Since T is no less
			 * than eps and ln eps is going to have greater
			 * magnitude than eps, underflow when calculating
			 * T * ln eps is not possible. */
			if (-deltah < t * log_wrapper_LVB_EPS)
			{
				pacc = 0.0;
				/* Call uni() even though its not required. It
				 * would have been called in LVB 1.0A, so this
				 * helps make results identical to results with
				 * that version. */
				(void)uni();
			}
			else	/* possibly accept the change */
			{
				pacc = exp_wrapper(-deltah / t);
				if (uni() < pacc)	/* do accept the change */
				{
					SwapTrees(&p_current_tree, &root, &p_proposed_tree, &proposed_tree_root);
					if (rcstruct.algorithm_selection == 2)
						w_changes_acc++;
					current_tree_length = proposed_tree_length;
					if (rcstruct.algorithm_selection == 1)
						changeAcc = 1;
				}
			}
		}
		proposed++;
#ifdef LVB_MAPREDUCE  
		MPI_Bcast(&proposed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif


		/* decide whether to reduce temperature */
		if (accepted >= maxaccept) {	/* enough new trees */

			failedcnt = 0;  /* this temperature a 'success' */
			dect = LVB_TRUE;
		}
		else if (proposed >= maxpropose) {	/* enough proposals */
			failedcnt++;
#ifdef LVB_MAPREDUCE  
			int check_stop = 0;
			if (misc->rank == 0 && failedcnt >= maxfail && t < FROZEN_T) check_stop = 1;
			MPI_Bcast(&check_stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (check_stop == 1) {
				MPI_Bcast(&failedcnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
#endif
			if (failedcnt >= maxfail && t < FROZEN_T)	/* system frozen */
			{
				/* Preliminary experiments yielded that the freezing
				 * criterion used in previous versions of LVB is not
				 * suitable for the new version and regularly results
				 * in premature termination of the search. An easy fix
				 * is to only apply the freezing criterion in the
				 * temperature ranges of previous versions of LVB
				 * (LVB_EPS < t < 10^-4). Future work will look at optimising
				 * maxpropose, maxaccept and maxfail directly. */
				break; /* end of cooling */
			}
			else {	/* system not frozen, so further decrease temp. */
				dect = LVB_TRUE;
			}
		}

		if (dect == LVB_TRUE)
		{
			t_n++;	/* originally n is 0 */

			if (rcstruct.cooling_schedule == 0)  /* Geometric cooling */
			{
				/* Ensure t doesn't go out of bound */
				ln_t = ((double)t_n) * log_wrapper_grad_geom + log_wrapper_t0;
				if (ln_t < log_wrapper_LVB_EPS) t = LVB_EPS;
				else t = pow_wrapper(grad_geom, (double)t_n) * t0; /* decrease the temperature */
				if (rcstruct.algorithm_selection == 1)
				{
					trops_probs[2] = t / t0;
					trops_probs[1] = (1 - trops_probs[2]) / 2;
					trops_probs[0] = trops_probs[1];
				}
			}
			else /* Linear cooling */
			{
				t = t0 - grad_linear * t_n;
				/* Make sure t doesn't go out of bounce */
				if (t < DBL_EPSILON || t <= LVB_EPS) t = LVB_EPS;
			}
			proposed = 0;
#ifdef LVB_MAPREDUCE  
			MPI_Bcast(&proposed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			accepted = 0;
			MPI_Bcast(&accepted, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			dect = LVB_FALSE;
#else
			accepted = 0;
			dect = LVB_FALSE;
#endif
			if (rcstruct.algorithm_selection == 2)
			{
				w_changes_prop = 0;
				w_changes_acc = 0;
			}
		}

		iter++;

		if (rcstruct.n_number_max_trees > 0 && treestack_ptr->next >= rcstruct.n_number_max_trees) {
			break;
		}

		if (rcstruct.algorithm_selection == 2)
		{
			if (changeAcc == 1) {
				trops_counter[trops_id]++;
			}
			else {
				for (int i = 0; i < 3; i++) {
					if (trops_id != i)
						trops_counter[i] = trops_counter[i] + 0.5;
				}
			}
		}
		if (rcstruct.verbose == LVB_TRUE)
			fprintf(pFile, "%ld\t%s\t%d\t%ld\t%lf\t%ld\n", iter, change, changeAcc, current_tree_length, t * 10000, treestack_ptr->next);



#ifdef LVB_MAPREDUCE  
		MPI_Barrier(MPI_COMM_WORLD);

	}
	print_sets(MSA, treestack_ptr, misc);
#else 
}
#endif

	/* free "local" dynamic heap memory */
	if (rcstruct.verbose == LVB_TRUE)
		fclose(pFile);
	free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	free(p_current_tree);
	free(p_proposed_tree);

	return best_tree_length;

} /* end Anneal() */


//pointer of rcstruct, instead of rcstruct itself
long Slave_Anneal(Dataptr MSA, TREESTACK* treestack_ptr, TREESTACK* treevo, const TREESTACK_TREE_NODES* const inittree, Parameters* p_rcstruct,
	long root, const double t0, const long maxaccept, const long maxpropose,
	const long maxfail, FILE* const lenfp, long* current_iter,
	Lvb_bool log_progress, int* p_n_state_progress, int* p_n_number_tried_seed, int myMPIid,
	Slave_record *record_slave, double * critical_t, MPI_Datatype mpi_recv_data, MPI_Datatype mpi_data_from_master)
	/* seek parsimonious tree from initial tree in inittree (of root root)
	 * with initial temperature t0, and subsequent temperatures obtained by
	 * multiplying the current temperature by (t1 / t0) ** n * t0 where n is
	 * the ordinal number of this temperature, after at least maxaccept changes
	 * have been accepted or maxpropose changes have been proposed, whichever is
	 * sooner;
	 * return the length of the best tree(s) found after maxfail consecutive
	 * temperatures have led to no new accepted solution;
	 * lenfp is for output of current tree length and associated details;
	 * *current_iter should give the iteration number at the start of this call and
	 * will be used in any statistics sent to lenfp, and will be updated on
	 * return */
{
	long accepted = 0;		/* changes accepted */
	Lvb_bool dect;		/* should decrease temperature */
	double deltah;		/* change in energy (1 - C.I.) */
	long tree_length_change;		/* change in length with new tree */
	long failedcnt = 0; 	/* "failed count" for temperatures */
	long iter = 0;		/* iteration of mutate/evaluate loop */
	long current_tree_length;			/* length of current tree */
	long best_tree_length;		/* best length found so far */
	long proposed_tree_length;		/* length of proposed new tree */
	double ln_t;		/* ln(current temperature) */
	long t_n = 0;		/* ordinal number of current temperature */
	double pacc;		/* prob. of accepting new config. */
	long proposed = 0;		/* trees proposed */
	double tree_minimum_length;		/* minimum length for any tree */
	long proposed_tree_root;		/* root of new configuration */
	double t = t0;		/* current temperature */
	double grad_geom = 0.99;		/* "gradient" of the geometric schedule */
	double grad_linear = 10 * LVB_EPS; /* gradient of the linear schedule */
	/* double grad_linear = 10 * LVB_EPS; */ /* gradient of the linear schedule */
	TREESTACK_TREE_NODES* p_current_tree;			/* current configuration */
	TREESTACK_TREE_NODES* p_proposed_tree;		/* proposed new configuration */
	long* p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
	long* p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
	int* p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */

#ifndef parallel
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int nFlag;
	MPI_Status mpi_status;
	MPI_Request request_handle_send, request_message_from_master;


	/* REND variables that could calculate immediately */
	SendInfoToMaster* p_data_info_to_master;
	p_data_info_to_master = (SendInfoToMaster*)malloc(sizeof(SendInfoToMaster));
	RecvInfoFromMaster* p_data_info_from_master;
	p_data_info_from_master = (RecvInfoFromMaster*)malloc(sizeof(RecvInfoFromMaster));



	int n = (maxpropose + INTERVAL_FOR_SPECIFIC_HEAT - 1) / INTERVAL_FOR_SPECIFIC_HEAT;
	double HI[n];
	int num_of_HI = 0;
	double Critical_temp = t0;//temperature at which specific heat is the maximum
	double tmp_sh=0, peak_sh=0;//specific heat

	p_data_info_to_master->start_temperature = t0;
	*p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_REPEAT; /* we consider always is necessary to repeat */ //来自函数形参

	Parameters rcstruct = *p_rcstruct;

	int first_commu_to_master = 1;
	p_data_info_from_master->n_is_to_continue = MPI_IS_TO_CONTINUE_ANNEAL;


	clock_t anneal_start, anneal_end;
	clock_t comm_start, comm_end;
	double t1 = 0, t2 = 0;//t1 is for whole anneal, t2 is just for communication

	Slave_record* tmp;
	tmp = (Slave_record*)malloc(sizeof(Slave_record));
	tmp->no_seed = *p_n_number_tried_seed;

	anneal_start = clock();
#endif


#ifdef LVB_MAPREDUCE  
	int* total_count;
	int check_cmp;
#endif

	/* variables that could calculate immediately */
	const double log_wrapper_LVB_EPS = log_wrapper(LVB_EPS);
	const double log_wrapper_grad_geom = log_wrapper(grad_geom);
	const double log_wrapper_t0 = log_wrapper(t0);
	/* REND variables that could calculate immediately */

	long w_changes_prop = 0;
	long w_changes_acc = 0;
	p_proposed_tree = treealloc(MSA, LVB_TRUE);
	p_current_tree = treealloc(MSA, LVB_TRUE);

	treecopy(MSA, p_current_tree, inittree, LVB_TRUE);	/* current configuration */

	alloc_memory_to_getplen(MSA, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	current_tree_length = getplen(MSA, p_current_tree, *p_rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
	dect = LVB_FALSE;		/* made LVB_TRUE as necessary at end of loop */

	lvb_assert(((float)t >= (float)LVB_EPS) && (t <= 1.0) && (grad_geom >= LVB_EPS) && (grad_linear >= LVB_EPS));

	best_tree_length = current_tree_length;

#ifdef LVB_MAPREDUCE  
	MPI_Bcast(&best_tree_length, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	PushCurrentTreeToStack(MSA, treestack_ptr, inittree, root, LVB_FALSE);
	misc->ID = treestack_ptr->next;
	misc->SB = 1;
	tree_setpush(MSA, inittree, root, mrTreeStack, misc);

#elif LVB_HASH
	CompareHashTreeToHashstack(MSA, treestack_ptr, inittree, root, LVB_FALSE);	/* init. tree initially best */
#else
	CompareTreeToTreestack(MSA, treestack_ptr, inittree, root, LVB_FALSE);	/* init. tree initially best */
#endif

	double trops_counter[3] = { 1,1,1 };
	double trops_probs[3] = { 0,0,0 };
	long trops_total = trops_counter[0] + trops_counter[1] + trops_counter[2];
	long trops_id = 0;


	if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {

		printf("\n--------------------------------------------------------");
		fprintf(lenfp, "\n  Temperature:   Rearrangement: TreeStack size: Length:\n");
		printf("--------------------------------------------------------\n");
	}

	/*Writing output to table.tsv*/
	FILE* pFile;
	char change[10] = "";
	if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {
		if (p_rcstruct->verbose == LVB_TRUE)
		{
			pFile = fopen("changeAccepted.tsv", "w");
			fprintf(pFile, "Iteration\tAlgorithm\tAccepted\tLength\tTemperature\tTreestack Size\n");
		}
	}
	tree_minimum_length = (double)MSA->min_len_tree;






	while (1) {
#ifndef parallel
		//添加HI: cost for specific heat
		if ((proposed % INTERVAL_FOR_SPECIFIC_HEAT) == 0)
		{
			HI[num_of_HI] = 1.0-(tree_minimum_length / (double)current_tree_length);
			num_of_HI++;
		}
#endif
		int changeAcc = 0;
		*current_iter += 1;
		/* occasionally re-root, to prevent influence from root position */
		if ((*current_iter % REROOT_INTERVAL) == 0) {
			root = arbreroot(MSA, p_current_tree, root);
			if ((log_progress == LVB_TRUE) && ((*current_iter % STAT_LOG_INTERVAL) == 0)) {
				lenlog(lenfp, treestack_ptr, *current_iter, current_tree_length, t);
			}

#ifndef parallel
			//Only Dynamic multi-instances needs to decide if being killed
			if ((p_rcstruct->parallel_selection == PARALLEL_DYNAMIC_AMI_KILL || p_rcstruct->parallel_selection == PARALLEL_DYNAMIC_AMI_KILL_PHASE_TRANSITION)
				&& ((*current_iter % STAT_LOG_INTERVAL) == 0))
			{
				//record communication time
				comm_start = clock();
				
				if(first_commu_to_master!=1)
				{
					MPI_Wait(&request_handle_send, MPI_STATUS_IGNORE);
				}

				if(first_commu_to_master!=1)
				{
					MPI_Wait(&request_message_from_master, MPI_STATUS_IGNORE);
				}

				//上一轮的temp，算出来不用kill;之前出错是因为第一轮没收到from_master故根本没有送与发，此时wait自然会sigment 11错误
				if (p_data_info_from_master->n_is_to_continue == MPI_IS_TO_CONTINUE_ANNEAL)
				{
					p_data_info_to_master->n_iterations = *current_iter;
					p_data_info_to_master->n_seed = p_rcstruct->seed;
					p_data_info_to_master->l_length = best_tree_length;
					p_data_info_to_master->stack_size = treestack_ptr->size;
					p_data_info_to_master->stack_next = treestack_ptr->next;
					p_data_info_to_master->temperature = t;
					p_data_info_to_master->n_finish_message = MPI_IS_TO_CONTINUE;//it means this message is from reaching STAT_LOG_INTERVAL
					p_data_info_to_master->Critical_temp = Critical_temp;
					p_data_info_to_master->peak_sh = peak_sh;

					//printf("Process:%d   send temperature:%.3g   iterations:%ld\n", myMPIid, t, *current_iter); 
					MPI_Issend(p_data_info_to_master, 1, mpi_recv_data, MPI_MAIN_PROCESS, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, &request_handle_send);
					MPI_Irecv(p_data_info_from_master, 1, mpi_data_from_master, MPI_MAIN_PROCESS, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, &request_message_from_master);
					//printf("After Process:%d   send temperature:%.3g   iterations:%ld\n", myMPIid, t, *current_iter);

					//end of this communication
					comm_end = clock();
					t2 += ((double)(comm_end - comm_start)) / CLOCKS_PER_SEC;

				}
				else if (p_data_info_from_master->n_is_to_continue == MPI_IS_TO_RESTART_ANNEAL || p_data_info_from_master->n_is_to_continue == MPI_IS_NOT_TO_RESTART)
				{
					*p_n_state_progress = MESSAGE_ANNEAL_KILLED;

					//end of this communication
					comm_end = clock();
					t2 += ((double)(comm_end - comm_start)) / CLOCKS_PER_SEC;
					break;
				}

				first_commu_to_master = 0;
			}
#endif

		}



		
		lvb_assert(t > DBL_EPSILON);

		/* mutation: alternate between the two mutation functions */
		proposed_tree_root = root;
		if (rcstruct.algorithm_selection == 2)
		{
			trops_total = trops_counter[0] + trops_counter[1] + trops_counter[2];
			trops_probs[0] = trops_counter[0] / trops_total;
			trops_probs[1] = trops_counter[1] / trops_total;
			trops_probs[2] = trops_counter[2] / trops_total;
		}

		if (rcstruct.algorithm_selection >= 1)
		{
			double random_val = uni();
			if (random_val < trops_probs[0]) {
				mutate_nni(MSA, p_proposed_tree, p_current_tree, root);	/* local change */
				strcpy(change, "NNI");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 0;
			}
			else if (random_val < trops_probs[0] + trops_probs[1]) {
				mutate_spr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
				strcpy(change, "SPR");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 1;
			}
			else {
				mutate_tbr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
				strcpy(change, "TBR");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 2;
			}
		}
		else
		{
			if (iter & 0x01) {
				mutate_spr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
				strcpy(change, "SPR");
			}
			else {
				mutate_nni(MSA, p_proposed_tree, p_current_tree, root);	/* local change */
				strcpy(change, "NNI");
			}
		}

		proposed_tree_length = getplen(MSA, p_proposed_tree, rcstruct, proposed_tree_root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
		lvb_assert(proposed_tree_length >= 1L);
		tree_length_change = proposed_tree_length - current_tree_length;
		deltah = (tree_minimum_length / (double)current_tree_length) - (tree_minimum_length / (double)proposed_tree_length);
		if (deltah > 1.0) deltah = 1.0; /* MinimumTreeLength() problem with ambiguous sites */

#ifdef LVB_MAPREDUCE  
		MPI_Bcast(&tree_length_change, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&deltah, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&proposed_tree_length, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

		if (tree_length_change <= 0)	/* accept the change */
		{
#ifdef LVB_MAPREDUCE  
			if (proposed_tree_length <= best_tree_length)	/* store tree if new */
			{
				if (proposed_tree_length < best_tree_length) {
					ClearTreestack(treestack_ptr);
					mrTreeStack->map(mrTreeStack, map_clean, NULL);
				}
				auto start_timer = high_resolution_clock::now();

				if (CompareMapReduceTrees(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, total_count,
					check_cmp, misc, mrTreeStack, mrBuffer) == 1) {
					accepted++;
					MPI_Bcast(&accepted, 1, MPI_LONG, 0, MPI_COMM_WORLD);
				}
				auto stop = high_resolution_clock::now();

				auto duration = duration_cast<microseconds>(stop - start_timer);

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE* timefunction = fopen("FunctionTimesMR", "a+");
					fprintf(timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
				}
			}
#else
			if (proposed_tree_length <= best_tree_length)	/* store tree if new */
			{

				if (proposed_tree_length < best_tree_length) {
					ClearTreestack(treestack_ptr);	/* discard old bests */
				}

#ifdef LVB_HASH
				auto start_timer = high_resolution_clock::now();
				if (CompareHashTreeToHashstack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1)
#else
				auto start_timer = high_resolution_clock::now();
				if (CompareTreeToTreestack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1)
#endif
				{
					accepted++;
				}
				auto stop = high_resolution_clock::now();

				auto duration = duration_cast<microseconds>(stop - start_timer);

#ifdef LVB_HASH

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE* timefunction = fopen("FunctionTimesHASH", "a+");
					fprintf(timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
				}
#else

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE* timefunction = fopen("FunctionTimesNP", "a+");
					fprintf(timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
				}
#endif
			}
#endif
			/* update current tree and its stats */
			current_tree_length = proposed_tree_length;
			SwapTrees(&p_current_tree, &root, &p_proposed_tree, &proposed_tree_root);

			/* very best so far */
			if (proposed_tree_length < best_tree_length) {
				best_tree_length = proposed_tree_length;
#ifdef LVB_MAPREDUCE  
				MPI_Bcast(&best_tree_length, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
			}
			if (rcstruct.algorithm_selection == 1)
				changeAcc = 1;

		}
		else	/* poss. accept change for the worse */
		{
			if (rcstruct.algorithm_selection == 2)
				w_changes_prop++;
			/* Mathematically,
			 *     Pacc = e ** (-1/T * deltaH)
			 *     therefore ln Pacc = -1/T * deltaH
			 *
			 * Computationally, if Pacc is going to be small, we
			 * can assume Pacc is 0 without actually working it
			 * out.
			 * i.e.,
			 *     if ln Pacc < ln eps, let Pacc = 0
			 * substituting,
			 *     if -deltaH / T < ln eps, let Pacc = 0
			 * rearranging,
			 *     if -deltaH < T * ln eps, let Pacc = 0
			 * This lets us work out whether Pacc will be very
			 * close to zero without dividing anything by T. This
			 * should prevent overflow. Since T is no less
			 * than eps and ln eps is going to have greater
			 * magnitude than eps, underflow when calculating
			 * T * ln eps is not possible. */
			if (-deltah < t * log_wrapper_LVB_EPS)
			{
				pacc = 0.0;
				/* Call uni() even though its not required. It
				 * would have been called in LVB 1.0A, so this
				 * helps make results identical to results with
				 * that version. */
				(void)uni();
			}
			else	/* possibly accept the change */
			{
				pacc = exp_wrapper(-deltah / t);
				if (uni() < pacc)	/* do accept the change */
				{
					SwapTrees(&p_current_tree, &root, &p_proposed_tree, &proposed_tree_root);
					if (rcstruct.algorithm_selection == 2)
						w_changes_acc++;
					current_tree_length = proposed_tree_length;
					if (rcstruct.algorithm_selection == 1)
						changeAcc = 1;
				}
			}
		}
		proposed++;
#ifdef LVB_MAPREDUCE  
		MPI_Bcast(&proposed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif


		/* decide whether to reduce temperature */
		if (accepted >= maxaccept) {	/* enough new trees */

			failedcnt = 0;  /* this temperature a 'success' */
			dect = LVB_TRUE;
		}
		else if (proposed >= maxpropose) {	/* enough proposals */
			failedcnt++;
#ifdef LVB_MAPREDUCE  
			int check_stop = 0;
			if (misc->rank == 0 && failedcnt >= maxfail && t < FROZEN_T) check_stop = 1;
			MPI_Bcast(&check_stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (check_stop == 1) {
				MPI_Bcast(&failedcnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
#endif
			if (failedcnt >= maxfail && t < FROZEN_T)	/* system frozen */
			{
				/* Preliminary experiments yielded that the freezing
				 * criterion used in previous versions of LVB is not
				 * suitable for the new version and regularly results
				 * in premature termination of the search. An easy fix
				 * is to only apply the freezing criterion in the
				 * temperature ranges of previous versions of LVB
				 * (LVB_EPS < t < 10^-4). Future work will look at optimising
				 * maxpropose, maxaccept and maxfail directly. */
				break; /* end of cooling */
			}
			else {	/* system not frozen, so further decrease temp. */
				dect = LVB_TRUE;
			}
		}

		if (dect == LVB_TRUE)
		{
			
#ifndef parallel
			/*这里要计算specific  heat了*/
//添加HI: cost for specific heat
			tmp_sh = Calculate_specific_heat(HI, num_of_HI, t);
			if (tmp_sh > peak_sh)
			{
				peak_sh = tmp_sh;
				Critical_temp = t;
			}
			num_of_HI = 0;
#endif


			t_n++;	/* originally n is 0 */

			if (rcstruct.cooling_schedule == 0)  /* Geometric cooling */
			{
				/* Ensure t doesn't go out of bound */
				ln_t = ((double)t_n) * log_wrapper_grad_geom + log_wrapper_t0;
				if (ln_t < log_wrapper_LVB_EPS) t = LVB_EPS;
				else t = pow_wrapper(grad_geom, (double)t_n) * t0; /* decrease the temperature */
				if (rcstruct.algorithm_selection == 1)
				{
					trops_probs[2] = t / t0;
					trops_probs[1] = (1 - trops_probs[2]) / 2;
					trops_probs[0] = trops_probs[1];
				}
			}
			else /* Linear cooling */
			{
				t = t0 - grad_linear * t_n;
				/* Make sure t doesn't go out of bounce */
				if (t < DBL_EPSILON || t <= LVB_EPS) t = LVB_EPS;
			}
			proposed = 0;
#ifdef LVB_MAPREDUCE  
			MPI_Bcast(&proposed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			accepted = 0;
			MPI_Bcast(&accepted, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			dect = LVB_FALSE;
#else
			accepted = 0;
			dect = LVB_FALSE;
#endif
			if (rcstruct.algorithm_selection == 2)
			{
				w_changes_prop = 0;
				w_changes_acc = 0;
			}
		}

		iter++;

		if (rcstruct.n_number_max_trees > 0 && treestack_ptr->next >= rcstruct.n_number_max_trees) {
			break;
		}

		if (rcstruct.algorithm_selection == 2)
		{
			if (changeAcc == 1) {
				trops_counter[trops_id]++;
			}
			else {
				for (int i = 0; i < 3; i++) {
					if (trops_id != i)
						trops_counter[i] = trops_counter[i] + 0.5;
				}
			}
		}
		if (rcstruct.verbose == LVB_TRUE)
			fprintf(pFile, "%ld\t%s\t%d\t%ld\t%lf\t%ld\n", iter, change, changeAcc, current_tree_length, t * 10000, treestack_ptr->next);
#ifdef LVB_MAPREDUCE  
		MPI_Barrier(MPI_COMM_WORLD);

	}
	print_sets(MSA, treestack_ptr, misc);
#else 
		
}
#endif



#ifndef old
	p_data_info_to_master->n_iterations = *current_iter;
	p_data_info_to_master->n_seed = p_rcstruct->seed;
	p_data_info_to_master->l_length = best_tree_length;
	p_data_info_to_master->stack_size = treestack_ptr->size;
	p_data_info_to_master->stack_next = treestack_ptr->next;
	p_data_info_to_master->temperature = t;
	p_data_info_to_master->Critical_temp = Critical_temp;
	p_data_info_to_master->peak_sh = peak_sh;

	comm_start = clock();
	//printf("\nno seed: %d, critical temperature:%lf,peak specific heat:%lf\n ", *p_n_number_tried_seed,Critical_temp, peak_sh);

	Slave_wait_final_message(&request_message_from_master, &request_handle_send, p_n_state_progress, p_data_info_from_master, p_rcstruct,
		p_n_number_tried_seed, p_data_info_to_master, critical_t, mpi_recv_data, mpi_data_from_master, first_commu_to_master);
	comm_end = clock();
	t2 += ((double)(comm_end - comm_start)) / CLOCKS_PER_SEC;

	anneal_end = clock();
	t1 = ((double)(anneal_end - anneal_start)) / CLOCKS_PER_SEC;
	
	tmp->anneal_once = t1;
	tmp->comm_cost = t2;
	tmp->next = record_slave->next;
	record_slave->next = tmp;

	//printf("\n\nwarn Cal  cost %lf\n\n", ti);
	

	free(p_data_info_to_master);
	free(p_data_info_from_master);
#endif




	/* free "local" dynamic heap memory */
	if (p_rcstruct->verbose == LVB_TRUE)
		fclose(pFile);
	free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	free(p_current_tree);
	free(p_proposed_tree);
	return best_tree_length;

} /* end Anneal() */








#ifdef LVB_MAPREDUCE
long GetSoln(Dataptr restrict MSA, Parameters rcstruct, long *iter_p, Lvb_bool log_progress,
				MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer)
#else
long GetSoln(Dataptr restrict MSA, Parameters rcstruct, long *iter_p, Lvb_bool log_progress)

#endif
/* get and output solution(s) according to parameters in rcstruct;
 * return length of shortest tree(s) found */
{
    static char fnam[LVB_FNAMSIZE];	/* current file name */
    long fnamlen;			/* length of current file name */
    long i;				/* loop counter */
    double t0;		/* SA cooling cycle initial temp */
    long maxaccept = MAXACCEPT_SLOW;	/* SA cooling cycle maxaccept */
    long maxpropose = MAXPROPOSE_SLOW;	/* SA cooling cycle maxpropose */
    long maxfail = MAXFAIL_SLOW;	/* SA cooling cycly maxfail */
    long treec;				/* number of trees found */
    long treelength = LONG_MAX;		/* length of each tree found */
    long initroot;			/* initial tree's root */
    FILE *sumfp;			/* best length file */
    FILE *resfp;			/* results file */
    TREESTACK_TREE_NODES *tree;			/* initial tree */
    Lvb_bit_length **enc_mat;	/* encoded data mat. */
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */
	#ifdef LVB_MAPREDUCE
	int *total_count;
	int check_cmp;
	#endif
	
     /* NOTE: These variables and their values are "dummies" and are no longer
     * used in the current version of LVB. However, in order to keep the
     * formatting of the output compatible with that of previous versions of
     * LVB these variables will continue to be used and written to the summary
     * files.  */
    long cyc = 0;	/* current cycle number */
    long start = 0;	/* current random (re)start number */


#ifndef parallel

	//TREESTACK_TREES best_tree;//the best tree from multiple runs
	TREESTACK best_treestack;//ts stand for treestack 
	best_treestack = CreateNewTreestack();


	long best_treelength=LONG_MAX;// the best final length from multiple runs
	int rank,nprocs;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	

	//Create MPI_Datatype for communication
	MPI_Datatype MPI_BRANCH, MPI_SLAVEtoMASTER, MPI_MASTERtoSLAVE;
	Create_MPI_Datatype(&MPI_BRANCH, &MPI_SLAVEtoMASTER, &MPI_MASTERtoSLAVE);


	double critical_t = -1;
	    int n_state_progress;		/* has a state to see if is necessary to repeat or finish */
	    int n_number_tried_seed_next = rank; 	/* it has the number of tried seed, (the number of next seed this proc will try */
	    int n_number_tried_seed = rank;	 	/* it has the number of tried seed */


		//Variables below are for collecting experimental data
		Info_record *record;
		Slave_record* record_slave;
		record_slave = (Slave_record *)malloc(sizeof(Slave_record));//header
		record_slave->next = NULL;
		int num_anneal = 0;
		clock_t Start, End;


		if (rcstruct.parallel_selection != PARALLEL_CLUSTER)
			if (rcstruct.nruns < (nprocs - 1))
				crash("number of runs should not be less than number of processes-1\n");
	
	if(rank==0)
	{
		record = (Info_record*)malloc(rcstruct.nruns * sizeof(Info_record));

		//Overall time taken
		Start = clock();

		//初始先送给rank 1
		MPI_Send(&rcstruct.seed, 1, MPI_INT, 1, MPI_TAG_PARAMS, MPI_COMM_WORLD);
		record[0].start = clock();

		//Rank 0 send initial seed for Slaves
		for (i = 2; i < nprocs; i++) 
		{
			/* send different to each proc */
			
			rcstruct.seed = get_other_seed_to_run_a_process();//如果有8个process，7个slave，会导致第七个新的seed不会被送出去，而是master自己留着
			MPI_Send(&rcstruct.seed, 1, MPI_INT, i, MPI_TAG_PARAMS, MPI_COMM_WORLD); 

			//record start time for each run corresponding to different seed
			record[i - 1].start = clock();
		}

		//Master start monitoring Slave, and do not need run SA
		if (rcstruct.parallel_selection != PARALLEL_CLUSTER)
		{
			int no_first_critical_temp = -1;// the number of the first seed that launched with critial temperature, if no seed launch with critical temp, it equals to -1
			double critical_temp = 0.0;
			get_temperature_and_control_process_from_other_process(nprocs, rcstruct.nruns, record, &no_first_critical_temp, &critical_temp, MPI_SLAVEtoMASTER, MPI_MASTERtoSLAVE);
			End= clock();
	
			Master_recv_record_from_Slave(record, rcstruct.nruns);
			write_final_results(record, rcstruct, no_first_critical_temp, critical_temp, ((double)(End - Start)) / CLOCKS_PER_SEC);
			free(record);
			goto finish;
		}

		free(record);
	}
	else
	{
		MPI_Recv(&rcstruct.seed, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_TAG_PARAMS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	rinit(rcstruct.seed);//reset seed


#endif
	
        /* dynamic "local" heap memory */
    tree = treealloc(MSA, LVB_TRUE);

    /* Allocation of the initial encoded MSA is non-contiguous because
     * this MSA isn't used much, so any performance penalty won't matter. */
    enc_mat = (Lvb_bit_length **) malloc((MSA->n) * sizeof(Lvb_bit_length *));
    for (i = 0; i < MSA->n; i++)
		enc_mat[i] = (Lvb_bit_length *) alloc(MSA->bytes, "state sets");
    DNAToBinary(MSA, enc_mat);

    /* open and entitle statistics file shared by all cycles
     * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */

    if (rcstruct.verbose == LVB_TRUE) {
		sumfp = clnopen(SUMFNAM, "w");
		fprintf(sumfp, "StartNo\tCycleNo\tCycInit\tCycBest\tCycTrees\n");
    }
    else{
        sumfp = NULL;
    }


#ifndef parallel
    	    
//start  loop of annealing
do{

	//reset seed at the end of the loop
    //rinit(rcstruct.seed);
	*iter_p = 0;

    /* determine starting temperature */
    PullRandomTree(MSA, tree);	/* initialise required variables */
    ss_init(MSA, tree, enc_mat);
    initroot = 0;

	if (critical_t != -1&& 
		(rcstruct.parallel_selection == PARALLEL_DYNAMIC_AMI_PHASE_TRANSITION|| rcstruct.parallel_selection == PARALLEL_DYNAMIC_AMI_KILL_PHASE_TRANSITION))
		t0 = critical_t;
	else
		t0 = StartingTemperature(MSA, tree, rcstruct, initroot, log_progress);

	n_state_progress = 0;	/* there's no state in beginning */
	n_number_tried_seed = n_number_tried_seed_next;//seed_next was received from alst anneal
	
#endif

    PullRandomTree(MSA, tree);	/* begin from scratch */
    ss_init(MSA, tree, enc_mat);
    initroot = 0;

    if (rcstruct.verbose) PrintStartMessage(start, cyc);
    	CheckStandardOutput();

    /* start cycles's entry in sum file
     * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions.  */

    if(rcstruct.verbose == LVB_TRUE) {
        alloc_memory_to_getplen(MSA, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
		fprintf(sumfp, "%ld\t%ld\t%ld\t", start, cyc, getplen(MSA, tree, rcstruct, initroot, p_todo_arr, p_todo_arr_sum_changes, p_runs));
		free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
		PrintInitialTree(MSA, tree, start, cyc, initroot);
    }

	#ifdef LVB_MAPREDUCE
		MPI_Barrier(MPI_COMM_WORLD);
		/* find solution(s) */
		treelength = Anneal(MSA, &treestack, &stack_treevo, tree, rcstruct, initroot, t0, maxaccept,
				maxpropose, maxfail, stdout, iter_p, log_progress, misc, mrTreeStack, mrBuffer );

		PullTreefromTreestack(MSA, tree, &initroot, &treestack, LVB_FALSE);
		CompareTreeToTreestack(MSA, &treestack, tree, initroot, LVB_FALSE);

		CompareMapReduceTrees(MSA, &treestack, tree, initroot, total_count,
					check_cmp, misc, mrTreeStack, mrBuffer);

	#else


	    /* find solution(s) */
	if(rcstruct.parallel_selection == PARALLEL_CLUSTER)
		treelength = Anneal(MSA, &treestack, &stack_treevo, tree, rcstruct, initroot, t0, maxaccept,
			maxpropose, maxfail, stdout, iter_p, log_progress);
	else
		treelength = Slave_Anneal(MSA, &treestack, &stack_treevo, tree, &rcstruct, initroot, t0, maxaccept,
					maxpropose, maxfail, stdout, iter_p, log_progress,&n_state_progress,&n_number_tried_seed_next,
				rank, record_slave, &critical_t, MPI_SLAVEtoMASTER, MPI_MASTERtoSLAVE);
		
    PullTreefromTreestack(MSA, tree, &initroot, &treestack, LVB_FALSE);

	#ifdef LVB_HASH
		CompareHashTreeToHashstack(MSA, &treestack, tree, initroot, LVB_FALSE);
	#else
		CompareTreeToTreestack(MSA, &treestack, tree, initroot, LVB_FALSE);
	#endif

    /* treelength = deterministic_hillclimb(MSA, &treestack, tree, rcstruct, initroot, stdout,
				iter_p, log_progress); */

	#endif

	/* log this cycle's solution and its details
	 * NOTE: There are no cycles anymore in the current version
     * of LVB. The code bellow is purely to keep the output consistent
     * with that of previous versions. */

    if (rcstruct.verbose == LVB_TRUE){
		fnamlen = sprintf(fnam, "%s_start%ld_cycle%ld", RESFNAM, start, cyc);
		lvb_assert(fnamlen < LVB_FNAMSIZE);	/* really too late */
		resfp = clnopen(fnam, "w");
		treec = PrintTreestack(MSA, &treestack, resfp, LVB_FALSE);
		clnclose(resfp, fnam);
		fprintf(sumfp, "%ld\t%ld\n", treelength, treec);

		/* won't use length summary file until end of next cycle */
		fflush(sumfp);
		if (ferror(sumfp)){
			crash("write error on file %s", SUMFNAM);
		}
    }


    if (rcstruct.verbose == LVB_TRUE) // printf("Ending start %ld cycle %ld\n", start, cyc);
    CheckStandardOutput();

    if (rcstruct.verbose == LVB_TRUE) clnclose(sumfp, SUMFNAM);



#ifdef parallel

/* compare the best tree with last best*/
    //printf("\n\n rank:%d  seed: %d;treelength=%d",rank, seed[i],treelength);

    if(treelength<best_treelength)
    {
	    //swap best treestack and the working treestacki, reduce complexity of creating new tree
	    long temp_size=best_treestack.size;			/* number of trees currently allocated for */
	    long temp_next=best_treestack.next;			/* next unused element of stack */
            TREESTACK_TREES *temp_stack=best_treestack.stack;

	    best_treestack.size=treestack.size;
	    best_treestack.next=treestack.next;
	    best_treestack.stack=treestack.stack;

	    treestack.size=temp_size;
	    treestack.next=temp_next;
	    treestack.stack=temp_stack;

	    best_treelength=treelength; 
    }
	//clear working treestack for next iteration
	ClearTreestack(&treestack);
#endif

#ifndef parallel
	num_anneal++;

	if (rcstruct.parallel_selection == PARALLEL_CLUSTER)
		break;

	//Below will also swap stack and best stack
    int is_conti = Slave_after_anneal_once(MSA,tree,n_state_progress,initroot,&treestack,&best_treestack,rank,iter_p,&treelength,&best_treelength, n_number_tried_seed_next,rcstruct);
	
    if(!is_conti)
	    break;


}while(1);
#endif

#ifndef old




/* "local" dynamic heap memory */
free(tree);
for (i = 0; i < MSA->n; i++) free(enc_mat[i]);
free(enc_mat);

if (rcstruct.parallel_selection != PARALLEL_CLUSTER)
	Slave_send_record_to_Master(num_anneal, record_slave->next);

	
finish:
free(record_slave);
//find best treestack
//output best treestack

/*Swap the best as return*/
if (rcstruct.parallel_selection != PARALLEL_CLUSTER)
{
	FreeTreestackMemory(MSA, &treestack);
	treestack.size = best_treestack.size;
	treestack.next = best_treestack.next;
	treestack.stack = best_treestack.stack;
	treelength = best_treelength;
}
else//In Cluster SA, best_treestack has no meaning 
	FreeTreestackMemory(MSA, &best_treestack);


//for master in parallel option of multi-instance, treelength == LONG_MAX
	Root_get_best_treestack(MSA, &treelength, 0, rank, nprocs, &treestack, MPI_BRANCH);

	//send record_slave to master

	
#endif

    return treelength;

} /* end GetSoln() */
