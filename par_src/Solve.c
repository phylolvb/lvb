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
long Anneal(Dataptr MSA, TREESTACK *treestack_ptr, TREESTACK *treevo, const TREESTACK_TREE_NODES *const inittree, Parameters rcstruct,
	long root, const double t0, const long maxaccept, const long maxpropose,
	const long maxfail, FILE *const lenfp, long *current_iter,
	Lvb_bool log_progress, MISC *misc, MapReduce *mrTreeStack, MapReduce *mrBuffer)

#else
long Anneal(Dataptr MSA, TREESTACK *treestack_ptr, TREESTACK *treevo, const TREESTACK_TREE_NODES *const inittree, Parameters rcstruct,
	long root, const double t0, const long maxaccept, const long maxpropose,
	const long maxfail, FILE *const lenfp, long *current_iter,
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
    TREESTACK_TREE_NODES *p_current_tree;			/* current configuration */
    TREESTACK_TREE_NODES *p_proposed_tree;		/* proposed new configuration */
    long *p_todo_arr; /* [MAX_BRANCHES + 1];	 list of "dirty" branch nos */
    long *p_todo_arr_sum_changes; /*used in openMP, to sum the partial changes */
    int *p_runs; 				/*used in openMP, 0 if not run yet, 1 if it was processed */

#ifdef old

    MPI_Request request_handle_send = 0, request_message_from_master = 0;
    MPI_Status mpi_status;

    int				nItems = 3;
    int          	blocklengths[3] = {2, 1, 1};
    MPI_Datatype 	types[3] = {MPI_INT, MPI_LONG, MPI_DOUBLE};
    MPI_Datatype 	mpi_recv_data;
    MPI_Aint     	displacements[3];
    displacements[0] = offsetof(SendInfoToMaster, n_iterations);
    displacements[1] = offsetof(SendInfoToMaster, l_length);
    displacements[2] = offsetof(SendInfoToMaster, temperature);
    MPI_Type_create_struct(nItems, blocklengths, displacements, types, &mpi_recv_data);
    MPI_Type_commit(&mpi_recv_data);

    nItems = 1;
    int          	blocklengths_2[1] = {3};
    MPI_Datatype 	types_2[1] = {MPI_INT};
    MPI_Datatype 	mpi_data_from_master;
    MPI_Aint     	displacements_2[1];
    displacements_2[0] = offsetof(RecvInfoFromMaster, n_seed);
    MPI_Type_create_struct(nItems, blocklengths_2, displacements_2, types_2, &mpi_data_from_master);
    MPI_Type_commit(&mpi_data_from_master);

    /* REND variables that could calculate immediately */
    SendInfoToMaster * p_data_info_to_master;
    p_data_info_to_master = (SendInfoToMaster *) malloc(sizeof(SendInfoToMaster));
    RecvInfoFromMaster * p_data_info_from_master;
    p_data_info_from_master = (RecvInfoFromMaster *) malloc(sizeof(RecvInfoFromMaster));

    *p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_REPEAT; /* we consider always is necessary to repeat */ //来自函数形参

#endif


	#ifdef LVB_MAPREDUCE  
		int *total_count;
	    int check_cmp;
	#endif

    /* variables that could calculate immediately */
    const double log_wrapper_LVB_EPS = log_wrapper(LVB_EPS);
    const double log_wrapper_grad_geom = log_wrapper(grad_geom);
    const double log_wrapper_t0 =  log_wrapper(t0);
    /* REND variables that could calculate immediately */

	long w_changes_prop = 0;
    long w_changes_acc = 0;  
    p_proposed_tree = treealloc(MSA, LVB_TRUE);
    p_current_tree = treealloc(MSA, LVB_TRUE);

    treecopy(MSA, p_current_tree, inittree, LVB_TRUE);	/* current configuration */

    alloc_memory_to_getplen(MSA, &p_todo_arr, &p_todo_arr_sum_changes, &p_runs);
    current_tree_length = getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
    dect = LVB_FALSE;		/* made LVB_TRUE as necessary at end of loop */

    lvb_assert( ((float) t >= (float) LVB_EPS) && (t <= 1.0) && (grad_geom >= LVB_EPS) && (grad_linear >= LVB_EPS));

    best_tree_length = current_tree_length;

	#ifdef LVB_MAPREDUCE  
		MPI_Bcast(&best_tree_length,  1, MPI_LONG, 0, MPI_COMM_WORLD);
		PushCurrentTreeToStack(MSA, treestack_ptr, inittree, root, LVB_FALSE);
		misc->ID = treestack_ptr->next;
		misc->SB = 1;
		tree_setpush(MSA, inittree, root, mrTreeStack, misc);
	
	#elif LVB_HASH
		CompareHashTreeToHashstack(MSA, treestack_ptr, inittree, root, LVB_FALSE);	/* init. tree initially best */
	#else
		CompareTreeToTreestack(MSA, treestack_ptr, inittree, root, LVB_FALSE);	/* init. tree initially best */ 
	#endif

	double trops_counter[3] = {1,1,1};
	double trops_probs[3] = {0,0,0};
	long trops_total = trops_counter[0]+trops_counter[1]+trops_counter[2];
	long trops_id = 0;
	

    if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {

	printf("\n--------------------------------------------------------");
	fprintf(lenfp, "\n  Temperature:   Rearrangement: TreeStack size: Length:\n");
	printf("--------------------------------------------------------\n");
	}
		
	/*Writing output to table.tsv*/
    FILE * pFile;
    char change[10]="";
    if ((log_progress == LVB_TRUE) && (*current_iter == 0)) {
		if (rcstruct.verbose == LVB_TRUE)
		{
	   pFile = fopen ("changeAccepted.tsv","w");
	   fprintf (pFile, "Iteration\tAlgorithm\tAccepted\tLength\tTemperature\tTreestack Size\n");
	}
	}
	tree_minimum_length = (double) MSA->min_len_tree;
#ifndef parallel
	//Control parameter for parallel_selection==1(cluster SA)
	long p=2000;//for cluster SA
	double r=0.1;
	long iter_cluster=-1;
	if(rcstruct.parallel_selection!=1)
		p=-1;

	int rank,nprocs;	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	

#endif	
	   while (1) {
#ifndef parallel
		   iter_cluster++;

		   if(p >= 0)
		   {
			   if(iter_cluster==p)
			   {
			//	printf("\nbefore one step, p=%ld, rank=%d, aurrent_treelength: %ld\n",p,rank, current_tree_length);

				   Bcast_best_partial_tree_to_root(MSA, current_tree_length , rank, nprocs, p_current_tree, &root);//Bcast best tree with root
				   //printf("\niter_cluster=%ld,p=%ld\n",iter_cluster,p);
				   iter_cluster=-1;
				   p-=r*(double)iter;
				   //if(rank==0)
					   //printf("\niter=%d, p=%ld,r=%d\n",iter, p,r);
				current_tree_length=getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);		
			//	printf("\nAfter one step, p=%ld, rank=%d, aurrent_treelength: %ld\n",p,rank, current_tree_length);
			   }

			   if(p<0)//reset annealing by best initial tree
			   {
				   accepted = 0;		/* changes accepted */
				   failedcnt = 0; 	/* "failed count" for temperatures */
				   iter = 0;		/* iteration of mutate/evaluate loop */
				   *current_iter=0;
				   current_tree_length=getplen(MSA, p_current_tree, rcstruct, root, p_todo_arr, p_todo_arr_sum_changes, p_runs);			/* length of current tree */
				   t_n = 0;		/* ordinal number of current temperature */
				   proposed = 0;		/* trees proposed */
				   t = StartingTemperature(MSA, p_current_tree, rcstruct, root, log_progress);		/* current temperature */

				   w_changes_prop = 0;
				   w_changes_acc = 0;  
				   for(int i=0;i<3;i++)
				   {
					   trops_counter[i]=1;
					   trops_probs[i]=0;
				   }
				   trops_total = trops_counter[0]+trops_counter[1]+trops_counter[2];
				   trops_id = 0;
					
				   if(rank==0)
				   printf("\n\n Generate the best initial tree, current_tree_length=%ld, root=%ld, new start temperature: %f\n\n",current_tree_length,root,t);
			   }
			   //ClearTreestack(treestack_ptr);
			   //CompareTreeToTreestack(MSA, treestack_ptr, p_current_tree, root, LVB_FALSE);
		   }

#endif
        int changeAcc = 0;
    	*current_iter += 1;
		/* occasionally re-root, to prevent influence from root position */
		if ((*current_iter % REROOT_INTERVAL) == 0){
			root = arbreroot(MSA, p_current_tree, root);
			if ((log_progress == LVB_TRUE) && ((*current_iter % STAT_LOG_INTERVAL) == 0)) {
        		lenlog(lenfp, treestack_ptr, *current_iter, current_tree_length, t);
        	}
		}



#ifdef old

		/* send temperature to the master process*/
		if(*current_iter%STAT_LOG_INTERVAL==0)//由iter控制同步
		{
					if (request_handle_send != 0)// wait temp send in last interval 
					{ 
						MPI_Wait(&request_handle_send, MPI_STATUS_IGNORE); //等待成功发送同步信息，是上一个interval发送的
					}

					p_data_info_to_master->n_iterations = *current_iter;//确定上一个intervel的message已发送成功，即可准备这一轮的信息
					p_data_info_to_master->n_seed = p_rcstruct->seed;
					p_data_info_to_master->l_length = lenbest;
					p_data_info_to_master->temperature = t;
					/* printf("Process:%d   send temperature:%.3g   iterations:%ld\n", myMPIid, t, *current_iter); */
					MPI_Isend(p_data_info_to_master, 1, mpi_recv_data, MPI_MAIN_PROCESS, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, &request_handle_send);//发送本轮的温度

					/* now get the message to continue or not, but need in second iteration... */
					if (request_message_from_master != 0) {//wait manage recv in last interval
						MPI_Wait(&request_message_from_master, MPI_STATUS_IGNORE);
						MPI_Test(&request_message_from_master, &nFlag, &mpi_status);
						if (nFlag == 0) { printf("ERROR, mpi waiting is not working File:%s  Line:%d\n", __FILE__, __LINE__); }
						if (nFlag == 1)
						{
							//Anealing有两种可能方向，killed和finished，后面再根据finished后message判断是否restart
							//根据从master接受二点message判断是否继续
							if (p_data_info_from_master->n_is_to_continue == MPI_IS_TO_RESTART_ANNEAL){	/*it's there and need to restart*/
								MPI_Cancel(&request_handle_send);//刚才的send没用了
								request_message_from_master = 0;
								request_handle_send = 0;
								*p_n_state_progress = MESSAGE_ANNEAL_KILLED;
								break;
							}
							/* otherwise need to proceed... */
						}
					}
					/* printf("Process:%d   receive management\n", myMPIid); */
					/* need to get other message to proceed... */
					MPI_Irecv(p_data_info_from_master, 1, mpi_data_from_master, MPI_MAIN_PROCESS, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, &request_message_from_master);//接受这一轮的manage信息
		}
#endif

		lvb_assert(t > DBL_EPSILON);

		/* mutation: alternate between the two mutation functions */
		proposed_tree_root = root;
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
        		mutate_nni(MSA, p_proposed_tree, p_current_tree, root);	/* local change */
			strcpy(change,"NNI");
			if (rcstruct.algorithm_selection == 2)
			trops_id = 0;
		}
    	else if (random_val < trops_probs[0] + trops_probs[1]) {
        	mutate_spr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
            strcpy(change,"SPR");
			if (rcstruct.algorithm_selection == 2)
			trops_id = 1;
        }
    	else {
    	   	mutate_tbr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
            strcpy(change,"TBR");
			if (rcstruct.algorithm_selection == 2)
			trops_id = 2;
        }
		  }
		else
		{
		if (iter & 0x01) { mutate_spr(MSA, p_proposed_tree, p_current_tree, root);	/* global change */
		strcpy(change,"SPR"); }
		else { mutate_nni (MSA, p_proposed_tree, p_current_tree, root);	/* local change */
		strcpy(change,"NNI"); }
		}
		
		proposed_tree_length = getplen(MSA, p_proposed_tree, rcstruct, proposed_tree_root, p_todo_arr, p_todo_arr_sum_changes, p_runs);
		lvb_assert (proposed_tree_length >= 1L);
		tree_length_change = proposed_tree_length - current_tree_length;
		deltah = (tree_minimum_length / (double) current_tree_length) - (tree_minimum_length / (double) proposed_tree_length);
		if (deltah > 1.0) deltah = 1.0; /* MinimumTreeLength() problem with ambiguous sites */

		#ifdef LVB_MAPREDUCE  
			MPI_Bcast(&tree_length_change, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(&deltah,   1, MPI_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(&proposed_tree_length,  1, MPI_LONG, 0, MPI_COMM_WORLD);
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
				
				if(CompareMapReduceTrees(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, total_count,
					check_cmp, misc, mrTreeStack, mrBuffer) == 1) {
						accepted++; 
						MPI_Bcast(&accepted, 1, MPI_LONG, 0, MPI_COMM_WORLD);	
				}
				auto stop = high_resolution_clock::now();

				auto duration = duration_cast<microseconds>(stop - start_timer);

				if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

				FILE *timefunction = fopen("FunctionTimesMR","a+");
				fprintf (timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
				fclose(timefunction);
				}
				}
				#else
								if (proposed_tree_length <= best_tree_length)	/* store tree if new */
			{
				/*printf("%ld\n", *current_iter);*/
				if (proposed_tree_length < best_tree_length) {
					ClearTreestack(treestack_ptr);	/* discard old bests */
				}
					#ifdef LVB_HASH
					auto start_timer = high_resolution_clock::now();
						if(CompareHashTreeToHashstack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1)
					#else
					auto start_timer = high_resolution_clock::now();
						if(CompareTreeToTreestack(MSA, treestack_ptr, p_proposed_tree, proposed_tree_root, LVB_FALSE) == 1)
					#endif
				{
					accepted++;
				}
				auto stop = high_resolution_clock::now();

					auto duration = duration_cast<microseconds>(stop - start_timer);

					#ifdef LVB_HASH

					if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE *timefunction = fopen("FunctionTimesHASH","a+");
					fprintf (timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
					fclose(timefunction);
					}
					#else

					if ((log_progress == LVB_TRUE) && ((treestack_ptr->next % 10) == 0)) {

					FILE *timefunction = fopen("FunctionTimesNP","a+");
					fprintf (timefunction, "%ld\t%ld\t%ld\t%ld\n", *current_iter, duration.count(), treestack_ptr->next, best_tree_length);
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
			MPI_Bcast(&best_tree_length,  1, MPI_LONG, 0, MPI_COMM_WORLD);
			#endif
			}
			if (rcstruct.algorithm_selection == 1)
			changeAcc = 1;

		}
		else	/* poss. accept change for the worse */
		{
			if (rcstruct.algorithm_selection == 2)
			w_changes_prop ++; 
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
				(void) uni();
			}
			else	/* possibly accept the change */
			{
				pacc = exp_wrapper(-deltah/t);
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
		MPI_Bcast(&proposed,  1, MPI_LONG, 0, MPI_COMM_WORLD);
		#endif

		
		/* decide whether to reduce temperature */
		if (accepted >= maxaccept){	/* enough new trees */
      
			failedcnt = 0;  /* this temperature a 'success' */
			dect = LVB_TRUE;
		}
		else if (proposed >= maxpropose){	/* enough proposals */
			failedcnt++;
			#ifdef LVB_MAPREDUCE  
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
			else{	/* system not frozen, so further decrease temp. */
				dect = LVB_TRUE;
			}
		}

		if (dect == LVB_TRUE)
		{
			t_n++;	/* originally n is 0 */

			if (rcstruct.cooling_schedule == 0)  /* Geometric cooling */
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
		#ifdef LVB_MAPREDUCE  
		MPI_Bcast(&proposed,  1, MPI_LONG, 0, MPI_COMM_WORLD);
		accepted = 0;
		MPI_Bcast(&accepted,  1, MPI_LONG, 0, MPI_COMM_WORLD);
		dect = LVB_FALSE;
		#else
		accepted = 0;
		dect = LVB_FALSE;
		#endif
		if (rcstruct.algorithm_selection == 2)
			{ w_changes_prop = 0;
        w_changes_acc = 0; 
			}
		}

		iter++;

		if (rcstruct.n_number_max_trees > 0 && treestack_ptr->next >= rcstruct.n_number_max_trees){
			break;
		}

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
	if (rcstruct.verbose == LVB_TRUE)
	fprintf (pFile, "%ld\t%s\t%d\t%ld\t%lf\t%ld\n", iter, change, changeAcc, current_tree_length, t*10000, treestack_ptr->next);
		#ifdef LVB_MAPREDUCE  
			MPI_Barrier(MPI_COMM_WORLD);

	    }
	    print_sets(MSA, treestack_ptr, misc);
		#else 
    }
		#endif

#ifdef old
if (request_message_from_master != 0) //如果frozen，会有manage没receive到，所以request要取消
	MPI_Cancel(&request_message_from_master);
if (request_handle_send != 0) 
	MPI_Cancel(&request_handle_send);

	/* Send FINISH message to master */
if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED || *p_n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT)
{	/* it's necessary to send this message */
	/* send the MPI_ID then the root can translate for the number of tried_seed */
	int n_finish_message = MPI_FINISHED;
	MPI_Isend(&n_finish_message, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, &request_handle_send);
	MPI_Wait(&request_handle_send, MPI_STATUS_IGNORE); /* need to do this because the receiver is asynchronous */

	/* need to wait for information if is necessary to run another */
	/* this one need to print the result */
	MPI_Status mpi_status;
	MPI_Recv(p_data_info_from_master, 1, mpi_data_from_master, MPI_MAIN_PROCESS, MPI_TAG_SEND_RESTART, MPI_COMM_WORLD, &mpi_status); /* this one waits until the master receive all confirmations */

	//    	if (mpi_status.MPI_ERROR != MPI_SUCCESS) {
	//    	   char error_string[BUFSIZ];
	//    	   int length_of_error_string;
	//    	   MPI_Error_string(mpi_status.MPI_ERROR, error_string, &length_of_error_string);
	//    	   printf("Process:%d   %s\n", myMPIid, error_string);
	//    	}

	//从master接受信息判断是否继续
	if (p_data_info_from_master->n_is_to_continue == MPI_IS_NOT_TO_RESTART){
		if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED) *p_n_state_progress = MESSAGE_ANNEAL_KILLED_AND_NOT_REPEAT;
		else *p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT;
	}
	else{
		p_rcstruct->seed = p_data_info_from_master->n_seed;  /* new seed */
		if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED) *p_n_state_progress = MESSAGE_ANNEAL_KILLED_AND_REPEAT;
		else *p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_REPEAT;
	}
	*p_n_number_tried_seed = p_data_info_from_master->n_process_tried;  /* it's necessary to create a file with trees */
}


    free(p_data_info_to_master);
    free(p_data_info_from_master);
#endif




    /* free "local" dynamic heap memory */
	if (rcstruct.verbose == LVB_TRUE)
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
	
#ifndef parallel

	//TREESTACK_TREES best_tree;//the best tree from multiple runs
	TREESTACK best_treestack;//ts stand for treestack 


	long best_treelength=LONG_MAX;// the best final length from multiple runs
	int rank,nprocs;
	//MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	int nruns_local=rcstruct.nruns/nprocs;
	if(rank<rcstruct.nruns%nprocs)//distribute the remainder to first several ones
		nruns_local++;


	if(rcstruct.parallel_selection==0)
		log_progress=LVB_FALSE;//for default, no progressloged
	else if(rcstruct.parallel_selection==1)//only rank==0, printf progress
		if(rank!=0)
			log_progress=LVB_FALSE;



	//printf("\n\nThere are %d processes; rank %d , runs=%d sss\n\n\n", nprocs,rank, nruns_local);


#endif
	
    /* NOTE: These variables and their values are "dummies" and are no longer
     * used in the current version of LVB. However, in order to keep the
     * formatting of the output compatible with that of previous versions of
     * LVB these variables will continue to be used and written to the summary
     * files.  */
    long cyc = 0;	/* current cycle number */
    long start = 0;	/* current random (re)start number */

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

#ifdef test
 
    TREESTACK_TREE_NODES *test;
    test=treealloc(MSA, LVB_TRUE);
    
     PullRandomTree(MSA, tree);	/* initialise required variables */
    ss_init(MSA, tree, enc_mat);
    initroot = 0;

    treecopy(MSA,test,tree,LVB_FALSE);

	t0 = StartingTemperature(MSA, tree, rcstruct, initroot, log_progress);


    PullRandomTree(MSA, tree);	/* begin from scratch */
    ss_init(MSA, tree, enc_mat);
    initroot = 0;
    CompareTreeToTreestack(MSA, &treestack, tree, 0, LVB_FALSE);

    if(CompareTreeToTreestack(MSA, &treestack, test, 0, LVB_FALSE)==0)
	    printf("\n\n************they are same\n\n");
    else
	    printf("\n\n***********they are different\n\n\n");
    crash("just for test");

#endif


#ifndef parallel
    int seed[nruns_local];
seed[0]= rcstruct.seed;
/****Cursory random seed****/
seed[0]+=rank;

if(rcstruct.parallel_selection==1)
	nruns_local=1;//for select 2, only run once



//start  loop of annealing
for(int i=0;i<nruns_local;i++)
{

    rinit(seed[i]);//reset seed
    /* determine starting temperature */
    PullRandomTree(MSA, tree);	/* initialise required variables */
    ss_init(MSA, tree, enc_mat);
    initroot = 0;


	t0 = StartingTemperature(MSA, tree, rcstruct, initroot, log_progress);

	//ProgressInfo.tem=t0;
#endif

#ifdef old

	n_state_progress = 0;	/* there's no state in beginning */
	n_number_tried_seed = n_number_tried_seed_next;
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
    treelength = Anneal(MSA, &treestack, &stack_treevo, tree, rcstruct, initroot, t0, maxaccept,
    maxpropose, maxfail, stdout, iter_p, log_progress);
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

#ifdef old
    /* 		Several possible outputs */
    /*		ANNEAL_FINISHED_AND_NOT_REPEAT		0x01
		ANNEAL_FINISHED_AND_REPEAT			0x02
		ANNEAL_KILLED_AND_REPEAT			0x03
		ANNEAL_KILLED_AND_NOT_REPEAT		0x04 */

    /* is is killed is not necesary any data */ //正常frozen
    if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT || n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT){

	    /* work on the trees */
	    int l_pop = PullTreefromTreestack(MSA, tree, &initroot, bstack_overall, LVB_FALSE);
	    if (l_pop == 0){
		    printf("\nProcess:%d    Error: can't pop any tree from Treestack.   Rearrangements tried: %ld\n", myMPIid, l_iterations);
	    }
	    else{
		    CompareTreeToTreestack(MSA, bstack_overall, tree, initroot, LVB_FALSE);
		    treelength = deterministic_hillclimb(MSA, bstack_overall, tree, rcstruct, initroot, stdout, &l_iterations, myMPIid, log_progress);
		    /* save it */
		    sprintf(file_name_output, "%s_%d", rcstruct.file_name_out, n_number_tried_seed); /* name of output file for this process */
		    outtreefp = clnopen(file_name_output, "w");
		    trees_output = PrintTreestack(MSA, bstack_overall, outtreefp, LVB_FALSE);
		    clnclose(outtreefp, file_name_output);
		    printf("\nProcess:%d   Rearrangements tried: %ld\n", myMPIid, l_iterations);
		    if (trees_output == 1L) { printf("1 most parsimonious tree of length %ld written to file '%s'\n", treelength, file_name_output); }
		    else { printf("%ld equally parsimonious trees of length %ld written to file '%s'\n", trees_output, treelength, file_name_output); }
	    }
    }

    //要repeat，或不
    if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT || n_state_progress == MESSAGE_ANNEAL_KILLED_AND_REPEAT){
	    l_iterations = 0;		/* start iterations from zero */
	    free(tree);
	    FreeTreestackMemory(bstack_overall);
	    printf("Process:%d   try seed number process:%d   new seed:%d", myMPIid, n_number_tried_seed_next, rcstruct.seed);
	    rinit(rcstruct.seed); /* at this point the structure has a need see passed by master process */ //用从master得到的新seed重启， Main.c:479:生产新seed
    }
    else{
	    /* Save the finish state file */
	    if (rcstruct.n_flag_save_read_states == DO_SAVE_READ_STATES){
		    save_finish_state_file(&rcstruct, myMPIid);
	    }
	    break; /* it is not necessary to repeat again */
    }

    check_stdout();
    if (rcstruct.verbose == LVB_TRUE) printf("Ending start %ld cycle %ld\n", start, cyc);
    if (rcstruct.verbose == LVB_TRUE) clnclose(sumfp, SUMFNAM);
#endif


#ifndef parallel

/* compare the best tree with last best*/
    //printf("\n\n rank:%d  seed: %d;treelength=%d",rank, seed[i],treelength);

    if(treelength<best_treelength)
    {
	    best_treestack.size=treestack.size;
	    best_treestack.next=treestack.next;
	    best_treestack.stack=treestack.stack;

            treestack = CreateNewTreestack();

	    best_treelength=treelength;
	    
    }
    else
	    ClearTreestack(&treestack);

    if(i<nruns_local-1)
    {
	    seed[i+1]=get_other_seed_to_run_a_process();
    }


}
#endif



    /* "local" dynamic heap memory */
    free(tree);
	for (i = 0; i < MSA->n; i++) free(enc_mat[i]);
    free(enc_mat);



#ifndef parallel
/*Swap the best as return*/
	FreeTreestackMemory(MSA, &treestack);
	treestack.size=best_treestack.size;
	treestack.next=best_treestack.next;
	treestack.stack=best_treestack.stack;
	treelength=best_treelength;

	Root_get_best_treestack(MSA, best_treelength, 0, rank, nprocs, &treestack);


#endif


    //Terminate MPI
    MPI_Finalize();

    return treelength;

} /* end GetSoln() */
