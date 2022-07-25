#include "Para_lib.h"
#include "mpi.h"
#include "RandomNumberGenerator.h"



int get_other_seed_to_run_a_process() {
    //return (int)(rand() % (unsigned long)MAX_SEED);
    return randpint(900000000);
}

void Bcast_best_partial_tree_to_root(Dataptr MSA, long best_treelength, int rank, int nprocs, TREESTACK_TREE_NODES* BranchArray, long *tree_root)
{
    //TREESTACK* treestack_ptr, long iteration, long length, double temperature
    //先算出来最小的length，
    //再bcast最佳的，
    //最后每个proc用这个最佳initial，用不同的随机数生成器开始真正完整的独立的SA
    //可以选择把得到最佳initial整个封装到Para_lib的一个函数里，然后再主函数的anneal里面，真正的anneal前调用（现在就先传送brancharry即可，穿不穿栈之后再考虑）

    long* all_len;

    all_len = (long*)alloc(sizeof(long) * nprocs, "all_len");

    MPI_Allgather(&best_treelength, 1, MPI_LONG, all_len, 1, MPI_LONG, MPI_COMM_WORLD);

    int index_max = 0;

    for (int i = 0; i < nprocs; i++)
    {
        if (all_len[i] < all_len[index_max])
        {
            index_max = i;
        }
    }
    //传送用的结构
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

    Lvb_bit_length** old_sitestate = (Lvb_bit_length**)alloc(sizeof(Lvb_bit_length*) * MSA->numberofpossiblebranches, "array saving old sitestate");//save original sitestate
    for (int i = 0; i < MSA->numberofpossiblebranches; i++)
    {
        old_sitestate[i] = BranchArray[i].sitestate;
    }

    //MPI_Bcast(BranchArray, MSA->numberofpossiblebranches, MPI_BRANCH, index_max, MPI_COMM_WORLD);//广播最屌的树
    //use MPI_BYTE might have data layout problem(wait improvement)
    MPI_Bcast(BranchArray, MSA->tree_bytes, MPI_BYTE, index_max, MPI_COMM_WORLD);//广播最屌的树
    MPI_Bcast(tree_root, 1, MPI_LONG, index_max, MPI_COMM_WORLD);//广播最屌的树

//Code below for potential memory alignment(just in case)
    //unsigned char* TreeArray_uchar_star = (unsigned char*)BranchArray;
    //unsigned char* ss0_start = TreeArray_uchar_star + MSA->numberofpossiblebranches * sizeof(TREESTACK_TREE_NODES);

    //MPI_Bcast(BranchArray, MSA->numberofpossiblebranches, MPI_BRANCH, index_max, MPI_COMM_WORLD);


#ifdef test
printf("\n\n\n**************UNITL NOW, SUCCESS***********\n");
printf("\n\n\n**************best_rank=%d, my_rank=%ld root %d\n",index_max,rank,*tree_root);
for(int i=0;i<MSA->numberofpossiblebranches;i++)
{
	printf("parent:%ld, left:%ld,  right: %ld, changes:%d \n",BranchArray[i].parent,BranchArray[i].left,BranchArray[i].right,BranchArray[i].changes);
}


printf("\n\n\n**************UNITL NOW, SUCCESS***********\n\n\n\n");
#endif


    for (int i = 0; i < MSA->numberofpossiblebranches; i++)
    {
        BranchArray[i].sitestate = old_sitestate[i];//restore sitestate
    }
    for (int i = MSA->n; i < MSA->numberofpossiblebranches; i++) 
    {
        BranchArray[i].sitestate[0] == 0U;// make dirty
    }

    


    //这里应该pulltree，把最优的树付给现在的tree，再pushtree
    //PullTreefromTreestack(MSA, brancharray);//没搞完



}

void Root_get_best_treestack(Dataptr MSA, long best_treelength_local, int root, int rank, int nprocs, TREESTACK* treestack)
{
    long * all_len;

    //printf("\n\n\nlocal best length %ld,  rank %d ----\n\n",best_treelength, rank);
    //if(rank==1)
	    //best_treelength=1;

    all_len = (long*)alloc(sizeof(long) * nprocs,"all_len");

    MPI_Allgather(&best_treelength_local, 1, MPI_LONG, all_len, 1, MPI_LONG, MPI_COMM_WORLD);
#ifdef test
if(rank==0)
{
	printf("\n");
   for(int i=0;i<nprocs;i++)
	  printf("%ld,",all_len[i]);
}
#endif

    int index_max = 0;
    int nsame = 0;//equivalently best
    /*
    int is_best=1;//there might be several best
    for (int i = 0; i < nprocs; i++)
    {
        if (best_treelength[i] < best_treelength[is_best])
        {
            is_best = 0;
            break;
        }
    }

    return is_best;
    */

    for (int i = 0; i < nprocs; i++)
    {
        if (all_len[i] < all_len[index_max])//返回第一个最短的
        {
            index_max = i;
            nsame = 1;

        }

        if (all_len[i] == all_len[index_max])
            nsame++;
    }

    //if (rank == index_max)
        //printf("\n\n\nbest rank:%d, length %ld-------------\n\n", rank, all_len[index_max]);

    Send_best_treestack_to_root(MSA, rank, 0, index_max, nprocs, treestack);
    

}



void  Send_best_treestack_to_root(Dataptr MSA, int rank,int root, int best_rank, int nprocs, TREESTACK * treestack)
{
    //每个process保留使用的seeds，得到的best_length, iteration, 最终_temperature, 
    //（记得每一次都要rinit(seed)），并把最好的一个seed对应的treestack，temperature，
    //给保留，然后交流得到最好的那个，再把最好的那个treestack返回回去
    //注意stack里面brancharray的sitestate是不拷贝进去的
    


    //传送用的结构
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

    




    if (best_rank == root)//root already is best
        return;
    else if (best_rank == rank)
    {
        printf("\n\nbest_rank %d, has stack size %ld and next %ld\n\n", rank, treestack->size, treestack->next);
        MPI_Ssend(&treestack->next, 1, MPI_LONG, root,0, MPI_COMM_WORLD);

        MPI_Request req_1[treestack->next], req_2[treestack->next];
        MPI_Request req_3[MSA->nsets], req_4[MSA->nsets];
	

        //先送cnt，等对面分配好空间，再送objset
        for (int i = 0; i < treestack->next; i++)
        {
            //挨个送stack里的树,但是这样整个送的话，原本的sitestate也会被覆盖，要修改，但是原本anneal过程中就不需要把sitestate传送过去
                //setstate里每个set装的是long类型的叶子的编号（在MSA->row里的序号），
                //sitestate_2直接建立的大数组，把存sitestate的空间与指针连续，就像brancharray一样
                //但真正stack里的p_sitestate不是大数组，详情见copy_sitestate(Dataptr restrict MSA, Objset *p_sitestate_1)


            MPI_Issend(&treestack->stack[i].root, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &req_1[i]);//root
            MPI_Issend(treestack->stack[i].tree, MSA->numberofpossiblebranches, MPI_BRANCH, 0, 0, MPI_COMM_WORLD, &req_2[i]);//送树,不需要送bitetate，
        }
        
        MPI_Waitall(treestack->next, req_1, MPI_STATUSES_IGNORE);
        MPI_Waitall(treestack->next, req_2, MPI_STATUSES_IGNORE);

#ifdef test


printf("\n\n\n**************UNITL NOW, SUCCESS***********\n");
printf("\n\n\n**************bestrank=%d, last root %d, last tree:\n",rank,treestack->stack[treestack->next-1].root);
for(int i=0;i<MSA->numberofpossiblebranches;i++)
{
	printf("parent:%ld, left:%ld,  right: %ld, changes:%d \n",treestack->stack[treestack->next-1].tree[i].parent,treestack->stack[treestack->next-1].tree[i].left,treestack->stack[treestack->next-1].tree[i].right,treestack->stack[treestack->next-1].tree[i].changes);
}


printf("\n\n\n**************UNITL NOW, SUCCESS***********\n\n\n\n");
MPI_Ssend(&treestack->next, 1, MPI_LONG, root,0, MPI_COMM_WORLD);

#endif
        for (int i = 0; i < treestack->next; i++)
        {
            //开始送cnt和objset
            for (int j = 0; j < MSA->nsets; j++)
            {
                MPI_Issend(&treestack->stack[i].p_sitestate[j].cnt, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &req_3[j]);//这边可以结构体化一口气，只送cnt不送*set
                //等待对面创好set用的空间
            }
            MPI_Waitall(MSA->nsets, req_3, MPI_STATUSES_IGNORE);

            for (int j = 0; j < MSA->nsets; j++)
            {
                MPI_Issend(treestack->stack[i].p_sitestate[j].set, treestack->stack[i].p_sitestate[j].cnt, MPI_LONG, 0, 0, MPI_COMM_WORLD,&req_3[j]);//这边可以结构体化一口气，只送cnt不送*set
            }
            MPI_Waitall(MSA->nsets, req_3, MPI_STATUSES_IGNORE);



#ifdef test


	    printf("\n\nSUCCESS FOR CNT, on send side\n");
 for (int j = 0; j < MSA->nsets; j++)
            {
                printf("cnt:%d;  ",treestack->stack[i].p_sitestate[j].cnt);
		printf("this set comprise: ");
		for(int x=0;x<treestack->stack[i].p_sitestate[j].cnt;x++)
			printf("%ld,",treestack->stack[i].p_sitestate[j].set[x]);
		printf("\n");

            }

MPI_Ssend(&treestack->next, 1, MPI_LONG, root,0, MPI_COMM_WORLD);
#endif            


        }
    }

    else if (rank == 0)
    {
        MPI_Status status;
        long cnt;
        MPI_Recv(&treestack->next, 1, MPI_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        //空间不够就创空间

        //直接分配一堆空间,next指向的空间不一定可用，所以size==next需要increase，size则表示practically有几颗树的空间已被分配
        if (treestack->next > treestack->size)
        {
            treestack->stack = (TREESTACK_TREES*)realloc(treestack->stack, (treestack->next) * sizeof(TREESTACK_TREES));
            if (treestack->stack == NULL) {
                crash("out of memory: cannot increase allocation for best tree stack to %ld elements", treestack->size);
            }
            
            /* allocate space within stack */
            for (int i= treestack->size; i < treestack->next; i++)
            {
 
                treestack->stack[i].tree = treealloc(MSA, LVB_FALSE);//LVB_FALSE即不送sitestate过去
                treestack->stack[i].root = -1;
                int j;
                /* set memory for sitestate */
                treestack->stack[i].p_sitestate = (Objset*)alloc(MSA->nsets * sizeof(Objset), "object set object arrays");
                for (j = 0; j < MSA->nsets; j++) {
                    treestack->stack[i].p_sitestate[j].set = NULL; 
                    treestack->stack[i].p_sitestate[j].cnt = UNSET;
                }
                treestack->size++;
            }
        }


        //printf("\n\nrank %d, after creating treestack space, has stack size %ld and next %ld\n\n", rank, treestack->size, treestack->next);
        

        MPI_Request req_1[treestack->next], req_2[treestack->next];
        MPI_Request req_3[MSA->nsets];

        //收树
        for (int i = 0; i < treestack->next; i++)
        {


            MPI_Irecv(&treestack->stack[i].root, 1, MPI_LONG, best_rank, 0, MPI_COMM_WORLD, &req_1[i]);
            MPI_Irecv(treestack->stack[i].tree, MSA->numberofpossiblebranches, MPI_BRANCH, best_rank, 0, MPI_COMM_WORLD, &req_2[i]);//送树,其实不需要送bitetate，但还是送，防止节点还要重新推导中间节点

        }

        MPI_Waitall(treestack->next, req_1, MPI_STATUSES_IGNORE);
        MPI_Waitall(treestack->next, req_2, MPI_STATUSES_IGNORE);
#ifdef test


MPI_Recv(&treestack->next, 1, MPI_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
printf("\n\n\n**************UNITL NOW, SUCCESSi  2***********\n");
printf("\n\n\n**************rank=%d, last root %d, last tree:\n",rank,treestack->stack[treestack->next-1].root);

for(int i=0;i<MSA->numberofpossiblebranches;i++)
{
	printf("parent:%ld, left:%ld,  right: %ld, changes:%d \n",treestack->stack[treestack->next-1].tree[i].parent,treestack->stack[treestack->next-1].tree[i].left,treestack->stack[treestack->next-1].tree[i].right,treestack->stack[treestack->next-1].tree[i].changes);
}

printf("\n\n\n**************UNITL NOW, SUCCESS  2***********\n\n\n\n");
#endif

        for (int i = 0; i < treestack->next; i++)
        {
            for (int j = 0; j < MSA->nsets; j++)
            {
                MPI_Irecv(&treestack->stack[i].p_sitestate[j].cnt, 1, MPI_LONG, best_rank, 0, MPI_COMM_WORLD, &req_3[j]);//这边可以结构体化一口气，只送cnt不送*set

            }
            MPI_Waitall(MSA->nsets, req_3, MPI_STATUSES_IGNORE);

            for (int j = 0; j < MSA->nsets; j++)
            {

               //等待对面创好set用的空间
               long to_copy = treestack->stack[i].p_sitestate[j].cnt * sizeof(long);
               if (treestack->stack[i].p_sitestate[j].set == NULL) {	// need to alloc memory
                   treestack->stack[i].p_sitestate[j].set = (long*)alloc(to_copy, "object set object arrays");

               }
               else if (cnt != treestack->stack[i].p_sitestate[j].cnt) {
                   treestack->stack[i].p_sitestate[j].set = (long*)realloc(treestack->stack[i].p_sitestate[j].set, to_copy);
                   if (treestack->stack[i].p_sitestate[j].set == NULL) {
                       crash("out of memory: cannot increase allocation for best sitestate %ld elements", to_copy);
                   }
               }


               MPI_Irecv(treestack->stack[i].p_sitestate[j].set, treestack->stack[i].p_sitestate[j].cnt, MPI_LONG, best_rank, 0, MPI_COMM_WORLD,&req_3[j]); 
            
            }
            MPI_Waitall(MSA->nsets, req_3, MPI_STATUSES_IGNORE);
            
#ifdef test
MPI_Recv(&treestack->next, 1, MPI_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

	    printf("\n\nSUCCESS FOR CNT\n");
 for (int j = 0; j < MSA->nsets; j++)
            {
                printf("cnt:%d;  ",treestack->stack[i].p_sitestate[j].cnt);
		printf("this set comprise: ");
		for(int x=0;x<treestack->stack[i].p_sitestate[j].cnt;x++)
			printf("%d,",treestack->stack[i].p_sitestate[j].set[x]);
		printf("\n");

            }

#endif            


        }
    }




}

#ifndef old
void Slave_interval_reached(MPI_Request *request_handle_send,SendInfoToMaster *p_data_info_to_master, MPI_Datatype mpi_recv_data, MPI_Request *request_message_from_master, RecvInfoFromMaster * p_data_info_from_master, MPI_Datatype mpi_data_from_master, int *p_n_state_progress)//mpi_recv_data is slave-> master, mpi_data_from_master is master->slave
{
	int nFlag;
	MPI_Status status;
	/* send temperature to the master process*/
	if (request_handle_send != 0) 
	{ 
		MPI_Wait(request_handle_send, MPI_STATUS_IGNORE); 
	}
	p_data_info_to_master->n_iterations = *current_iter;
	p_data_info_to_master->n_seed = p_rcstruct->seed;
	p_data_info_to_master->l_length = lenbest;
	p_data_info_to_master->temperature = t;
	/* printf("Process:%d   send temperature:%.3g   iterations:%ld\n", myMPIid, t, *current_iter); */
	MPI_Isend(p_data_info_to_master, 1, mpi_recv_data, MPI_MAIN_PROCESS, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, request_handle_send);
	/* now get the message to continue or not, but need in second iteration... */
	if (request_message_from_master != 0) 
	{
		MPI_Wait(request_message_from_master, MPI_STATUS_IGNORE);
		MPI_Test(request_message_from_master, &nFlag, &mpi_status);
		if (nFlag == 0) { printf("ERROR, mpi waiting is not working File:%s  Line:%d\n", __FILE__, __LINE__); }
		if (nFlag == 1)
		{
			if (p_data_info_from_master->n_is_to_continue == MPI_IS_TO_RESTART_ANNEAL)
			{	/*it's there and need to restart*/
				MPI_Cancel(request_handle_send);
				*request_message_from_master = 0;
				*request_handle_send = 0;
				*p_n_state_progress = MESSAGE_ANNEAL_KILLED;
				break;
			}
			/* otherwise need to proceed... */
		}
	}
	/* printf("Process:%d   receive management\n", myMPIid); */
	/* need to get other message to proceed... */
	MPI_Irecv(p_data_info_from_master, 1, mpi_data_from_master, MPI_MAIN_PROCESS, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, request_message_from_master);
}
}



void Slave_wait_final_message(MPI_Request *request_message_from_master, MPI_Request *request_handle_send,int *p_n_state_progress, RecvInfoFromMaster *p_data_info_from_master, MPI_Datatype mpi_data_from_master, Parameters *p_rcstruct, int *p_n_number_tried_seed, long best_tree_length)
{
	if (request_message_from_master != 0) 
		MPI_Cancel(request_message_from_master);
	if (request_handle_send != 0) 
		MPI_Cancel(request_handle_send);

	/* Send FINISH message to master */
	if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED || *p_n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT)
	{	/* it's necessary to send this message */
		/* send the MPI_ID then the root can translate for the number of tried_seed */
		int n_finish_message = MPI_FINISHED;
		MPI_Isend(n_finish_message, 1, MPI_INT, MPI_MAIN_PROCESS, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, request_handle_send);
		MPI_Wait(&request_handle_send, MPI_STATUS_IGNORE); /* need to do this because the receiver is asynchronous */

		/*For collecting result sned send length if frozen, sned -1 of length if killed*/
		MPI_Isend(&best_tree_length, 1, MPI_LONG, MPI_MAIN_PROCESS, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, request_handle_send);
		MPI_Wait(&request_handle_send,MPI_STATUS_IGNORE);


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

		if (p_data_info_from_master->n_is_to_continue == MPI_IS_NOT_TO_RESTART)
		{
			if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED) 
				*p_n_state_progress = MESSAGE_ANNEAL_KILLED_AND_NOT_REPEAT;
			else 
				*p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT;
		}
		else{
			p_rcstruct->seed = p_data_info_from_master->n_seed;  /* new seed */
			if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED) 
				*p_n_state_progress = MESSAGE_ANNEAL_KILLED_AND_REPEAT;
			else 
				*p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_REPEAT;
		}
		*p_n_number_tried_seed = p_data_info_from_master->n_process_tried;  /* it's necessary to create a file with trees */
	}
}


void Slave_after_anneal_once(Dataptr MSA, TREESTACK_TREE_NODES *tree, int *n_state_progress, long * initroot, TREESTACK *treestack, 
		TREESTACK *best_treestack, int myMPIid, long *l_iterations, long *treelength, long *best_treelength)
{
    /* 		Several possible outputs */
    /*		ANNEAL_FINISHED_AND_NOT_REPEAT		0x01
		ANNEAL_FINISHED_AND_REPEAT			0x02
		ANNEAL_KILLED_AND_REPEAT			0x03
		ANNEAL_KILLED_AND_NOT_REPEAT		0x04 */

    /* is is killed is not necesary any data */ //frozen, need to out put the tree found
    if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT || n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT)
    {

	    /* work on the trees */
	    int l_pop = PullTreefromTreestack(MSA, tree, &initroot, treestack, LVB_FALSE);
	    if (l_pop == 0)
	    {
		    printf("\nProcess:%d    Error: can't pop any tree from Treestack.   Rearrangements tried: %ld\n", myMPIid, l_iterations);
	    }
	    else
	    {
		    CompareTreeToTreestack(MSA, treestack, tree, initroot, LVB_FALSE);
		    //no climbing in the newest version
		    //treelength = deterministic_hillclimb(MSA, bstack_overall, tree, rcstruct, initroot, stdout, &l_iterations, myMPIid, log_progress);
		    //Send and keep the best tree stack this process ever have
		    
		    
		    if(*treelength<*best_treelength)
		    {
			    //swap best treestack and the working treestacki, reduce complexity of creating new tree
			    long temp_size=best_treestack->size;			/* number of trees currently allocated for */
			    long temp_next=best_treestack->next;			/* next unused element of stack */
			    TREESTACK_TREES *temp_stack=best_treestack->stack;

			    best_treestack->size=treestack->size;
			    best_treestack->next=treestack->next;
			    best_treestack->stack=treestack->stack;

			    treestack->size=temp_size;
			    treestack->next=temp_next;
			    treestack->stack=temp_stack;

			    *best_treelength=*treelength; 
		    }
		    //clear working treestack for next iteration
		    ClearTreestack(treestack);

	    }
    }

    //repeat, or not
    if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT || n_state_progress == MESSAGE_ANNEAL_KILLED_AND_REPEAT)
    {
	    *l_iterations = 0;		/* start iterations from zero */
	    free(tree);
	    FreeTreestackMemory(treestack);
	    printf("Process:%d   try seed number process:%d   new seed:%d", myMPIid, n_number_tried_seed_next, rcstruct.seed);
	    rinit(rcstruct.seed); /* at this point the structure has a need see passed by master process */ //repeat using seeds from master (Main.c:479 produce new seed)
	    return 1;
    }
    else
    {
	    /* Save the finish state file */
	    //break; /* it is not necessary to repeat again */
	    return 0;//not repeat
    }


}

void get_temperature_and_control_process_from_other_process(int num_procs, int n_seeds_to_try){

		int nProcessFinished = 0;
		MPI_Request *pHandleTemperatureRecv;
		MPI_Request *pHandleFinishProcess;
		MPI_Request *pHandleManagementProcess;
		MPI_Status mpi_status;
		IterationTemperature *p_calc_iterations;
		int *pIntFinishedProcessChecked, *pIntFinishedProcess, *pIntProcessControlNumberRunning;
		int nFlag, i, n_seeds_tried = 0;


		/*collect final result length of -1 means killed*/
		SendInfoToMaster final_results[n_seeds_need_to_try]; 
		int idx=0;



		pHandleTemperatureRecv = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for temperature...");
		pHandleFinishProcess = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for finish process...");
		pHandleManagementProcess = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for management process...");
		pIntProcessControlNumberRunning = (int *) alloc(num_procs * sizeof(int), "Buffer with index between process and tried process, but the ones that are running...");
		pIntFinishedProcess = (int *) alloc(num_procs * sizeof(int), "Buffer for finished processes...");
		pIntFinishedProcessChecked = (int *) alloc(num_procs * sizeof(int), "Buffer for processes states...");
		memset(pIntFinishedProcessChecked, 0, num_procs * sizeof(int)); 	/* control the state of a process, but for all seeds tried */
		memset(pIntProcessControlNumberRunning, 0, num_procs * sizeof(int));		/* has the index number of the process that is running */
		memset(pIntFinishedProcess, 0, num_procs * sizeof(int));		/* control if a specific process is finished */
		for (i = 0; i < num_procs; i++) { *(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN; } /* all the process start with this state */

		/* structure to use sending temperature and number of interactions to master process */
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
		MPI_Datatype 	mpi_send_data;
		MPI_Aint     	displacements_2[1];
		displacements_2[0] = offsetof(RecvInfoFromMaster, n_seed);
		MPI_Type_create_struct(nItems, blocklengths_2, displacements_2, types_2, &mpi_send_data);
		MPI_Type_commit(&mpi_send_data);

		/* structure with temperature and interaction data */
		SendInfoToMaster **p_info_temperature;
		p_info_temperature = (SendInfoToMaster **) malloc(sizeof(SendInfoToMaster *) * (num_procs - 1));
		for (i = 0; i < num_procs - 1; i++){
			*(p_info_temperature + i) = (SendInfoToMaster *) malloc(sizeof(SendInfoToMaster));
		}

		/* structure with management data */
		RecvInfoFromMaster **p_info_manage;
		p_info_manage = (RecvInfoFromMaster **) malloc(sizeof(RecvInfoFromMaster *) * (num_procs - 1));
		for (i = 0; i < num_procs - 1; i++){
			*(p_info_manage + i) = (RecvInfoFromMaster *) malloc(sizeof(RecvInfoFromMaster));
		}

		/* first get handles for all receivers */
		for (i = 1; i < num_procs; i++) {
			MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
			MPI_Irecv(pIntFinishedProcess + i, 1, MPI_INT, i, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, pHandleFinishProcess + i);
			*(pHandleManagementProcess + i) = 0;
			*(pIntProcessControlNumberRunning + i) = i - 1;   /* the index number of ones that are running */
		}

		/* alloc memory main calc iterations*/
		p_calc_iterations = get_alloc_main_calc_iterations();
		n_seeds_tried = num_procs - 1;	/* we try always these number of seeds in beginning */
		while (1){

			for (i = 1; i < num_procs; i++) {

				/* if is equal to ANNEAL_STOP_PROCESS it is not necessary to do anything else to this process */
				/* because the limit of tried process is reached */
				if (*(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE && *(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS){

					MPI_Test(pHandleTemperatureRecv + i, &nFlag, &mpi_status);
					if (nFlag == 1){	/* message received */
						printf("Process:%d    main process getting temperature and length from process:%d   temperature:%-15.8g   length:%ld   iteration:%d\n", mpi_status.MPI_SOURCE,
							i, (*(p_info_temperature + (i - 1)))->temperature, (*(p_info_temperature + (i - 1)))->l_length,
							    (*(p_info_temperature + (i - 1)))->n_iterations);

						/* all the processes are synchronized by the number of iterations... */
						/* upload the data to the structure */
						/* make calculations to how many standard deviations has */
						/* if it is less than X standard deviations can continue */
						/* add data from this process to structure calc temperature iterations */
						add_temperature_cal_iterations(p_calc_iterations, *(p_info_temperature + (i - 1)), *(pIntProcessControlNumberRunning + i));
						if (is_possible_to_continue(p_calc_iterations, (*(p_info_temperature + (i - 1)))->temperature,
								(*(p_info_temperature + (i - 1)))->n_iterations, (*(p_info_temperature + (i - 1)))->l_length, num_procs, 0)){

							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_CONTINUE_ANNEAL;
							/* need to check this one, to see if the previous message was delivered... */
							if (*(pHandleManagementProcess + i) != 0){
								/* printf("Process:%d   wait management from process:%d\n", 0, i); */
								MPI_Wait(pHandleManagementProcess + i, MPI_STATUS_IGNORE);
							}
						}
						else{
							/* if the previous one was not delivered it can be canceled... */
							if (*(pHandleManagementProcess + i) != 0{
								MPI_Test(pHandleManagementProcess + i, &nFlag, &mpi_status);
								if (nFlag == 0) MPI_Cancel(pHandleManagementProcess + i);
							}

							*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE;
							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;
						}

						/* send management message */
						/*printf("Process:%d   send management to process:%d\n", 0, i); */
						MPI_Isend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);

						/* launch other receive temperature */
						/* printf("Process:%d   receive temperature from process:%d\n", 0, i); */
						if ((*(p_info_manage + (i - 1)))->n_is_to_continue != MPI_IS_TO_RESTART_ANNEAL)
							MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
					}
				}

				/* test the ones that finish, and all of them need to send this message...  */
				if (*(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS){
					MPI_Test(pHandleFinishProcess + i, &nFlag, &mpi_status);
					if (nFlag == 1 && *(pIntFinishedProcess + i) == MPI_FINISHED){
						MPI_Test(pHandleTemperatureRecv + i, &nFlag, &mpi_status);
						if (nFlag == 0) MPI_Cancel(pHandleTemperatureRecv + i);
						if (*(pHandleManagementProcess + i) != 0) MPI_Cancel(pHandleManagementProcess + i);
						printf("Process:%d    finish\n", i);
						
						/*receive final tree length*/
						MPI_Recv(&final_results[idx], 1, mpi_recv_data, i, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, &mpi_status); /* this one waits until the master receive all confirmations */
						idx++;
						printf("seed %d generates result: final temp iter\n",)






						if (n_seeds_tried < n_seeds_to_try){		/* try another seed */
							n_seeds_tried += 1;									/* try another seed */
							*(pIntProcessControlNumberRunning + i) = n_seeds_tried;
							*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN;

							(*(p_info_manage + (i - 1)))->n_seed = get_other_seed_to_run_a_process();
							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;
							(*(p_info_manage + (i - 1)))->n_process_tried = n_seeds_tried;
							printf("Process:%d    send to restart\n", i);
							MPI_Send(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_SEND_RESTART, MPI_COMM_WORLD);

							/* launch other wait messages */
							MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
							MPI_Irecv(pIntFinishedProcess + i, 1, MPI_INT, i, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, pHandleFinishProcess + i);
						}
						else{
							*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS;
							(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_NOT_TO_RESTART;
							printf("Process:%d    isn't to restart\n", i);
							MPI_Send(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_SEND_RESTART, MPI_COMM_WORLD);

							/* there's no need other MPI messages */
							*(pHandleTemperatureRecv + i) = 0;
							*(pHandleFinishProcess + i) = 0;
						}
						nProcessFinished += 1;
					}
				}
			}

			/* is it finished, all of then?  */
			if (nProcessFinished == n_seeds_to_try){
				/* All processes are finished */
				/* cancel the requests of temperature */
				for (i = 1; i < num_procs; i++){
					if (*(pHandleTemperatureRecv + i) != 0) MPI_Cancel(pHandleTemperatureRecv + i);
					if (*(pHandleManagementProcess + i) != 0) MPI_Cancel(pHandleManagementProcess + i);
				}
				break;
			}
		}

		/* release memory main calc iterations */
		release_main_calc_iterations(p_calc_iterations);

		/* release other memory */
		for (i = 0; i < num_procs - 1; i++){ free(*(p_info_temperature + i)); }
		free(p_info_temperature);
		for (i = 0; i < num_procs - 1; i++){ free(*(p_info_manage + i)); }
		free(p_info_manage);

		free(pHandleTemperatureRecv);		/* all asynchronous messages */
		free(pHandleManagementProcess);		/* all asynchronous messages */
		free(pHandleFinishProcess);			/* all asynchronous messages */
		free(pIntFinishedProcess);
		free(pIntProcessControlNumberRunning);
		free(pIntFinishedProcessChecked);
	}

#endif
