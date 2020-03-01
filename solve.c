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

/* ========== solve.c - solving functions ========== */

#include "lvb.h"

	static void lenlog(FILE *lengthfp, Treestack *bstackp, int myMPIid, long iteration, long length, double temperature)

	/* write a message to file pointer lengthfp; iteration gives current iteration;
	 * crash verbosely on write error */
	{
		fprintf(lengthfp, "%-15.8g%-15ld%-16ld%-15ld\n", temperature, iteration, bstackp->next, length);

		if (ferror(lengthfp)){
			crash("file error when logging search progress");
		}

	} /* end lenlog() */

	long deterministic_hillclimb(Dataptr restrict matrix, Treestack *bstackp, const Branch *const inittree,
			Params rcstruct, long root, FILE * const lenfp, long *current_iter, int myMPIid, Lvb_bool log_progress
			#ifdef MAP_REDUCE_SINGLE 
			, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer
			#endif
			)
	/* perform a deterministic hill-climbing optimization on the tree in inittree,
	 * using NNI on all internal branches until no changes are accepted; return the
	 * length of the best tree found; current_iter should give the iteration number
	 * at the start of this call and will be used in any statistics sent to lenfp,
	 * and will be updated on return */
	{
	    long i;													/* loop counter */
	    long j;													/* loop counter */
	    long todo_cnt = 0;										/* count of internal branches */
	    long len;												/* current length */
	    long lendash;											/* length of proposed new config */
	    long rootdash = root;									/* root of proposed new config */
	    long deltalen;											/* change in length */
	    Lvb_bool newtree;										/* accepted a new configuration */
	    Branch *p_current_tree;									/* current configuration */
	    Branch *p_proposed_tree;								/* proposed new configuration */
		static Lvb_bool leftright[] = { LVB_FALSE, LVB_TRUE };  // to loop through left and right
	    static long todo[MAX_BRANCHES];							/* array of internal branch numbers */
		// unsigned int *todo;										// array of internal branch numbers	 
	    long *p_todo_arr; 										/* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
	    long *p_todo_arr_sum_changes; 							/*used in openMP, to sum the partial changes */
	    int *p_runs; 											/*used in openMP, 0 if not run yet, 1 if it was processed */

#ifdef MAP_REDUCE_SINGLE
	    int  *total_count;
	    int check_cmp;
#endif

	    /* "local" dynamic heap memory */
	    p_current_tree  = treealloc(matrix, LVB_TRUE);
	    p_proposed_tree = treealloc(matrix, LVB_TRUE);
	    treecopy(matrix, p_current_tree, inittree, LVB_TRUE);      /* current configuration */
	    alloc_memory_to_getplen(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	    len = getplen(matrix, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
		
		// todo = alloc(matrix->nbranches * sizeof(unsigned int), "old parent alloc");
		
	    /* identify internal branches */
	    for (i = matrix->n; i < matrix->nbranches; i++) todo[todo_cnt++] = i;
	    lvb_assert(todo_cnt == matrix->nbranches - matrix->n);
	    do {
			newtree = LVB_FALSE;
			for (i = 0; i < todo_cnt; i++) {
				for (j = 0; j < 2; j++) {
					mutate_deterministic(matrix, p_proposed_tree, p_current_tree, root, todo[i], leftright[j]);
					lendash = getplen(matrix, p_proposed_tree, rcstruct, rootdash, p_todo_arr, p_todo_arr_sum_changes, p_runs);

					lvb_assert (lendash >= 1L);
					deltalen = lendash - len;


#ifdef MAP_REDUCE_SINGLE
					MPI_Bcast(&deltalen, 1, MPI_LONG, 0,    MPI_COMM_WORLD);
					MPI_Bcast(&lendash,  1, MPI_LONG, 0,    MPI_COMM_WORLD);
#endif
						if (deltalen <= 0) {
							if (deltalen < 0)  /* very best so far */
						{
							treestack_clear(bstackp);
								#ifdef MAP_REDUCE_SINGLE
								treestack_push_only(matrix, bstackp, p_proposed_tree, rootdash, LVB_FALSE);
								misc->ID = bstackp->next;
								misc->SB = 1;
								tree_setpush(matrix, p_proposed_tree, rootdash, mrTreeStack, misc);
								#endif
							len = lendash;
						}
#ifdef MAP_REDUCE_SINGLE
else {

						  misc->SB = 0;
						  tree_setpush(matrix, p_proposed_tree, rootdash, mrBuffer, misc);
						  mrBuffer->add(mrTreeStack);
						  mrBuffer->collate(NULL);

						  misc->count = (int *) alloc( (bstackp->next+1) * sizeof(int), "int array for tree comp using MR");
						  total_count = (int *) alloc( (bstackp->next+1) * sizeof(int), "int array for tree comp using MR");
						  for(int i=0; i<=bstackp->next; i++) misc->count[i] = 0;
						  mrBuffer->reduce(reduce_count, misc);

						  for(int i=0; i<=bstackp->next; i++) total_count[i] = 0;
						  MPI_Reduce( misc->count, total_count, bstackp->next+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

						  check_cmp = 1;
						  if (misc->rank == 0) { /* sum to root process */
							for(int i=1; i<=bstackp->next; i++) {
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
								tree_setpush(matrix, p_proposed_tree, rootdash, mrBuffer, misc);
								mrTreeStack->add(mrBuffer);
								treestack_push_only(matrix, bstackp, p_proposed_tree, rootdash, LVB_FALSE);
								misc->ID = bstackp->next;
#else
						  if (treestack_push(matrix, bstackp, p_proposed_tree, rootdash, LVB_FALSE) == 1) {
#endif
						  newtree = LVB_TRUE;
						  treeswap(&p_current_tree, &root, &p_proposed_tree, &rootdash);
						  }
#ifdef MAP_REDUCE_SINGLE
						  free(misc->count);
						  free(total_count);
						  }
#endif
						  }
						  if ((log_progress == LVB_TRUE) && ((*current_iter % rcstruct.STAT_LOG_INTERVAL) == 0)) {
							  #ifdef MAP_REDUCE_SINGLE
							  if(misc->rank == 0)
							  #endif
						  lenlog(lenfp, bstackp, myMPIid, *current_iter, len, 0);
						  }
						*current_iter += 1;
				}
			}
	    } while (newtree == LVB_TRUE);	

	    /* free "local" dynamic heap memory */
	    free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
	    free(p_current_tree);
	    free(p_proposed_tree);

	    return len;
	}

long anneal(Dataptr restrict matrix, Treestack *bstackp, Treestack *treevo, const Branch *const inittree, Params rcstruct, Params *p_rcstruct,
		long root, const double t0, const long maxaccept, const long maxpropose,
		const long maxfail, FILE *const lenfp, long *current_iter, int myMPIid, Lvb_bool log_progress


#ifdef MAP_REDUCE_SINGLE
		, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer
#endif
)

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
	 * return;
	 * if this is a restart after checkpointing, restore_fnam must
	 * give the file name of the checkpoint file; if it is not a restart, restore_fnam
	 * should be NULL */
	{
	    long accepted = 0;		/* changes accepted */
	    Lvb_bool dect;		/* should decrease temperature */
	    double deltah;		/* change in energy (1 - C.I.) */
	    long deltalen;		/* change in length with new tree */
	    long failedcnt = 0; 	/* "failed count" for temperatures */
	    long iter = 0;		/* iteration of mutate/evaluate loop */
	    long len;			/* length of current tree */
	    long lenbest;		/* bet length found so far */
	    long lendash;		/* length of proposed new tree */
	    double ln_t;		/* ln(current temperature) */
	    long t_n = 0;		/* ordinal number of current temperature */
	    double pacc;		/* prob. of accepting new config. */
	    long proposed = 0;		/* trees proposed */
	    double r_lenmin;		/* minimum length for any tree */
	    long rootdash;		/* root of new configuration */
	    double t = t0;				/* current temperature */
	    Branch *p_current_tree;			/* current configuration */
	    Branch *p_proposed_tree;			/* proposed new configuration */
	    long *p_todo_arr; 				/* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
	    long *p_todo_arr_sum_changes; 		/*used in openMP, to sum the partial changes */
	    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */
		double grad_geom = GRAD_GEOM;			/* "gradient" of the geometric schedule */
	    double grad_linear = GRAD_LINEAR; 	/* gradient of the linear schedule */
		long lenmin;						// minimum length for any tree

#ifdef MAP_REDUCE_SINGLE
	    int *total_count;
	    int check_cmp;
#endif

	    /* variables that could calculate immediately */
	    const double log_wrapper_LVB_EPS = log_wrapper(LVB_EPS);
	    const double log_wrapper_grad_geom = log_wrapper(grad_geom);
	    const double log_wrapper_t0 = log_wrapper(t0);
    	/* REND variables that could calculate immediately */
	    
		alloc_memory_to_getplen(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);

		long w_changes_prop = 0;
	    long w_changes_acc = 0;
		p_proposed_tree = treealloc(matrix, LVB_TRUE);
    	p_current_tree = treealloc(matrix, LVB_TRUE);


		treecopy(matrix, p_current_tree, inittree, LVB_TRUE);	/* current configuration */
		alloc_memory_to_getplen(matrix, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
		len = getplen(matrix, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);

		dect = LVB_FALSE;		/* made LVB_TRUE as necessary at end of loop */
		lvb_assert( ((float) t >= (float) LVB_EPS) && (t <= 1.0) && (grad_geom >= LVB_EPS) && (grad_linear >= LVB_EPS));
		lenbest = len;
			treestack_push(matrix, bstackp, inittree, root, LVB_FALSE);	/* init. tree initially best */

		if(rcstruct.algorithm_selection == 2)
			treestack_push(matrix, treevo, inittree, root, LVB_FALSE);

		double trops_counter[3] = {1,1,1};
		double trops_probs[3] = {0,0,0};
		long trops_total = trops_counter[0]+trops_counter[1]+trops_counter[2];
		long trops_id = 0;

		#ifdef MAP_REDUCE_SINGLE
		MPI_Bcast(&lenbest,  1, MPI_LONG, 0, MPI_COMM_WORLD);
		misc->ID = bstackp->next;
		misc->SB = 1;
		tree_setpush(matrix, inittree, root, mrTreeStack, misc);
		#endif

		if ((log_progress == LVB_TRUE) && (*current_iter == 0))
		#ifdef MAP_REDUCE_SINGLE
		if(misc->rank == 0) 
		#endif
			fprintf(lenfp, "\nTemperature:   Rearrangement: TreeStack size: Length:\n");

		if (rcstruct.verbose == LVB_TRUE)
        	fprintf(lenfp, "\nTemperature:   Rearrangement: TreeStack size: Length:\n");

		/*Writing output to table.tsv*/
    			FILE * pFile;
				pFile = (FILE *) alloc(sizeof(FILE), "alloc FILE structure");
    			char change[10]="";
    			if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {
					if (rcstruct.verbose == LVB_TRUE)
				{
	  			pFile = fopen ("changeAccepted.tsv","w");
	   			fprintf (pFile, "Iteration\tAlgorithm\tAccepted\tLength\tTemperature\n");
				}
				}
				
				lenmin = getminlen(matrix);
				r_lenmin = (double) lenmin;
				
			while (1) {
				int changeAcc = 0;
				*current_iter += 1;
			
			/* occasionally re-root, to prevent influence from root position */
			if ((*current_iter % REROOT_INTERVAL) == 0) {
				root = arbreroot(matrix, p_current_tree, root);
				if ((log_progress == LVB_TRUE) && ((*current_iter % rcstruct.STAT_LOG_INTERVAL) == 0)) {
					#ifdef MAP_REDUCE_SINGLE
	        		if(misc->rank == 0)  
					#endif
					lenlog(lenfp, bstackp, myMPIid, *current_iter, len, t);
					}
			}
			lvb_assert(t > DBL_EPSILON);
			rootdash = root;

			if (rcstruct.algorithm_selection ==2)
			{
				trops_total = trops_counter[0]+trops_counter[1]+trops_counter[2];
				trops_probs[0]=trops_counter[0]/trops_total;
				trops_probs[1]=trops_counter[1]/trops_total;
				trops_probs[2]=trops_counter[2]/trops_total;
			}

		 	if (rcstruct.algorithm_selection >= 1)
		  	{
				double random_val = uni();
 				if (random_val < trops_probs[0]) {
        			mutate_nni(matrix, p_proposed_tree, p_current_tree, root);	/* local change */
				strcpy(change,"NNI");
				if (rcstruct.algorithm_selection == 2)
					trops_id = 0;
			}
    		else if (random_val < trops_probs[0] + trops_probs[1]) {
        		mutate_spr(matrix, p_proposed_tree, p_current_tree, root);	/* global change */
            strcpy(change,"SPR");
			if (rcstruct.algorithm_selection == 2)
				trops_id = 1;
        	}
    		else {
    	   		mutate_tbr(matrix, p_proposed_tree, p_current_tree, root);	/* global change */
            strcpy(change,"TBR");
			if (rcstruct.algorithm_selection == 2)
				trops_id = 2;
        	}
		  		}
			else
			{
				if (iter & 0x01) { mutate_spr(matrix, p_proposed_tree, p_current_tree, root);	/* global change */
			strcpy(change,"SPR"); }
				else { mutate_nni (matrix, p_proposed_tree, p_current_tree, root);	/* local change */
			strcpy(change,"NNI"); }
			}

			lendash = getplen(matrix, p_proposed_tree, *p_rcstruct, rootdash, p_todo_arr, p_todo_arr_sum_changes, p_runs);
			lvb_assert (lendash >= 1L);
			deltalen = lendash - len;
			deltah = (r_lenmin / (double) len) - (r_lenmin / (double) lendash);
			if (deltah > 1.0) deltah = 1.0; /* getminlen() problem with ambiguous sites */
			
			#ifdef MAP_REDUCE_SINGLE
			MPI_Bcast(&deltalen, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(&deltah,   1, MPI_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(&lendash,  1, MPI_LONG, 0, MPI_COMM_WORLD);
			#endif

			if (deltalen <= 0)	/* accept the change */
			{
				if (lendash <= lenbest)	/* store tree if new */
				{
			#ifdef MAP_REDUCE_SINGLE
			if (lendash < lenbest) {
						treestack_clear(bstackp);
						mrTreeStack->map( mrTreeStack, map_clean, NULL );

						treestack_push_only(matrix, bstackp, p_proposed_tree, rootdash, LVB_FALSE);
						misc->ID = bstackp->next;

					    misc->SB = 1;
						tree_setpush(matrix, p_proposed_tree, rootdash, mrTreeStack, misc);

						accepted = 1;
						MPI_Bcast(&accepted,  1, MPI_LONG, 0, MPI_COMM_WORLD);

						MPI_Barrier(MPI_COMM_WORLD);
					} else {
	//		if(misc->rank == 0) cerr << "checking treecmp! " << endl;

	                    misc->SB = 0;
						tree_setpush(matrix, p_proposed_tree, rootdash, mrBuffer, misc);
						mrBuffer->add(mrTreeStack);
						mrBuffer->collate(NULL);

						misc->count = (int *) alloc( (bstackp->next+1) * sizeof(int), "int array for tree comp using MR");
						total_count = (int *) alloc( (bstackp->next+1) * sizeof(int), "int array for tree comp using MR");

						for(int i=0; i<=misc->ID; i++) misc->count[i] = 0;
				//		cerr << "Reduce ********************* " << endl;
				//		cerr << "Reduce ********************* " << endl;
						mrBuffer->reduce(reduce_count, misc);
				//		cerr << "END Reduce ********************* " << endl;
				//		cerr << "END Reduce ********************* " << endl;
						for(int i=0; i<=misc->ID; i++) total_count[i] = 0;
						MPI_Reduce( misc->count, total_count, misc->ID+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

						check_cmp = 1;
						if (misc->rank == 0) {
							for(int i=1; i<=misc->ID; i++) {
							//	if (misc->nsets == total_count[i]) {
								if (total_count[0] == total_count[i]) {
									check_cmp = 0;
									break;
								}
							}
	//  	cerr << "Sanity Check = " << endl;
	//  	cerr << "Sanity Check = " << endl;
	//	for(int i=0; i<=misc->ID; i++) cerr << "Sanity Check = " << i << " " << check_cmp << " " << total_count[i] << " " << total_count[0] << endl;

						}

						MPI_Barrier(MPI_COMM_WORLD);
						MPI_Bcast(&check_cmp, 1, MPI_INT, 0,    MPI_COMM_WORLD);
						if (check_cmp == 1) {

							treestack_push_only(matrix, bstackp, p_proposed_tree, rootdash, LVB_FALSE);
	                                                misc->ID = bstackp->next;

							misc->SB = 1;
							tree_setpush(matrix, p_proposed_tree, rootdash, mrBuffer, misc);
							mrTreeStack->add(mrBuffer);
							accepted++;
							MPI_Bcast(&accepted,  1, MPI_LONG, 0, MPI_COMM_WORLD);
						}

						free(misc->count);
						free(total_count);

					}

				}
			#else
				/*printf("%ld\n", *current_iter);*/
					if (lendash < lenbest) treestack_clear(bstackp);	/* discard old bests */
					if (treestack_push(matrix, bstackp, p_proposed_tree, rootdash, LVB_FALSE) == 1){
						accepted++;
					}
				}
			#endif

			/* update current tree and its stats */
				len = lendash;
				treeswap(&p_current_tree, &root, &p_proposed_tree, &rootdash);

				/* very best so far */
				if (lendash < lenbest) lenbest = lendash;
				if (rcstruct.algorithm_selection ==1)
				changeAcc = 1;

				#ifdef MAP_REDUCE_SINGLE
				MPI_Bcast(&lenbest,  1, MPI_LONG, 0, MPI_COMM_WORLD);
				#endif

				}
				else	/* poss. accept change for the worse */
				{
					if (rcstruct.algorithm_selection == 2)
					w_changes_prop ++;
				if (-deltah < t * log_wrapper_LVB_EPS)
				{
					pacc = 0.0;
					/* Call uni() even though its not required. It
					 * would have been called in LVB 1.0A, so this
					 * helps make results identical to results with
					 * that version. */
					(void) uni();
				}
				else	/* possibly accept the change */
				{
					pacc = exp_wrapper(-deltah/t);
					if (uni() < pacc)	/* do accept the change */
					{
						treeswap(&p_current_tree, &root, &p_proposed_tree, &rootdash);
						if (rcstruct.algorithm_selection == 2)
						w_changes_acc++;
						len = lendash;
						if (rcstruct.algorithm_selection == 1)
						changeAcc = 1;
					}
				}
			}
			proposed++;

#ifdef MAP_REDUCE_SINGLE
			MPI_Bcast(&proposed,  1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
			/* decide whether to reduce temperature */
			if (accepted >= maxaccept){	/* enough new trees */
				failedcnt = 0;  /* this temperature a 'success' */
				dect = LVB_TRUE;
			}
			else if (proposed >= maxpropose){	/* enough proposals */
				failedcnt++;

#ifdef MAP_REDUCE_SINGLE
				int check_stop = 0;
				if (misc->rank == 0 && failedcnt >= maxfail && t < FROZEN_T) check_stop = 1;
				MPI_Bcast(&check_stop,  1, MPI_INT, 0, MPI_COMM_WORLD);
				if (check_stop == 1) {
					MPI_Bcast(&failedcnt,  1, MPI_LONG, 0, MPI_COMM_WORLD);
					MPI_Bcast(&t,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				}
#endif
				if (failedcnt >= maxfail && t < FROZEN_T)	/* system frozen */
				{
					break; /* end of cooling, break from while(1) */
				}
				else{	/* system not frozen, so further decrease temp. */
					dect = LVB_TRUE;
				}
			}

			if (dect == LVB_TRUE)
			{
				t_n++;	/* originally n is 0 */
				if (rcstruct.cooling_schedule == 0) // Geometric cooling
				{
					/* Ensure t doesn't go out of bound */
					ln_t = ((double) t_n) * log_wrapper_grad_geom + log_wrapper_t0;
					if (ln_t < log_wrapper_LVB_EPS) t = LVB_EPS;
					else t = pow_wrapper(grad_geom, (double) t_n) * t0; /* decrease the temperature */
				if (rcstruct.algorithm_selection == 1)
					{
    			    trops_probs[2] = t/t0;
 			        trops_probs[1] = (1 - trops_probs[2])/2;
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
				accepted = 0;
				dect = LVB_FALSE;
				if (rcstruct.algorithm_selection == 2)
				{ w_changes_prop = 0;
        		w_changes_acc = 0;
				}
			}
			iter++;

	if (rcstruct.algorithm_selection == 2)
	{
	if (changeAcc == 1) {
	    trops_counter[trops_id]++;
	}
	else {
	    for (int i=0; i < 3; i++) {
	     if (trops_id != i)
	        trops_counter[i] = trops_counter[i] + 0.5;
	    }
	}
	}
	#ifdef MAP_REDUCE_SINGLE
	MPI_Barrier(MPI_COMM_WORLD);
	#endif
	if (rcstruct.verbose == LVB_TRUE)
	fprintf (pFile, "%ld\t%s\t%d\t%ld\t%lf\t%f\n", iter, change, changeAcc, len, t*10000, (float) r_lenmin/len);
    }

	#ifdef MAP_REDUCE_SINGLE
	print_sets(matrix, bstackp, misc);
	#endif

    /* free "local" dynamic heap memory */
	if (rcstruct.verbose == LVB_TRUE)
	fclose(pFile);
	free_memory_to_getplen(&p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
    free(p_current_tree);
    free(p_proposed_tree);
    return lenbest;

} /* end anneal() */