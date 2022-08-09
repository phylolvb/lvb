#include "RandomNumberGenerator.h"
#include "InteractionTemperature.h"
#include "Para_lib.h"



int get_other_seed_to_run_a_process() {
    //return (int)(rand() % (unsigned long)MAX_SEED);
    return randpint(900000000);
}


void Create_MPI_Datatype(MPI_Datatype *MPI_BRANCH, MPI_Datatype *MPI_SLAVEtoMASTER, MPI_Datatype* MPI_MASTERtoSLAVE)
{
	//Bcast_best_partial_tree 以及Send_best_treestack_to_root 需要MPI_BRANCH
	//get_temperature_and_control_process_from_other_process 与 SLAVE_anneal需要 send message 与recvmessage
	//这个函数可以在 getsol里调用
	//ps: 原本的slave里的branch还是单独创立的，因为如果想改成getsol传过去，还要改Solve.h太麻烦
		//传送用的结构
	int pad_size = sizeof(Lvb_bit_length*);
	int nItems = 2;
	int          	blocklengths_1[2] = { 4, pad_size };//content of sitestate doesn't matter
	MPI_Datatype 	types_1[2] = { MPI_LONG, MPI_BYTE };
	//MPI_Datatype 	MPI_BRANCH;
	MPI_Aint     	displacements_1[2];
	displacements_1[0] = offsetof(TREESTACK_TREE_NODES, parent);
	displacements_1[1] = offsetof(TREESTACK_TREE_NODES, sitestate);
	MPI_Type_create_struct(nItems, blocklengths_1, displacements_1, types_1, MPI_BRANCH);
	MPI_Type_commit(MPI_BRANCH);


	/* structure to use sending temperature and number of interactions to master process */
	nItems = 3;
	int          	blocklengths_2[3] = { 3, 1, 4 };
	MPI_Datatype 	types_2[3] = { MPI_INT, MPI_LONG, MPI_DOUBLE };
	MPI_Aint     	displacements_2[3];
	displacements_2[0] = offsetof(SendInfoToMaster, n_iterations);
	displacements_2[1] = offsetof(SendInfoToMaster, l_length);
	displacements_2[2] = offsetof(SendInfoToMaster, temperature);
	MPI_Type_create_struct(nItems, blocklengths_2, displacements_2, types_2, MPI_SLAVEtoMASTER);
	MPI_Type_commit(MPI_SLAVEtoMASTER);

	nItems = 2;
	int          	blocklengths_3[2] = { 3,1 };
	MPI_Datatype 	types_3[2] = { MPI_INT, MPI_DOUBLE};
	MPI_Aint     	displacements_3[2];
	displacements_3[0] = offsetof(RecvInfoFromMaster, n_seed);
	displacements_3[1] = offsetof(RecvInfoFromMaster, Critical_temp);
	MPI_Type_create_struct(nItems, blocklengths_3, displacements_3, types_3, MPI_MASTERtoSLAVE);
	MPI_Type_commit(MPI_MASTERtoSLAVE);

}


double Calculate_specific_heat(double* HI, int num_of_HI, double temp)
{

	double l_HI = 0.0;
	double l_avg, l_std, l_ss;

	/* calc average */
	for (int i=0; i<num_of_HI; i++) 
	{
		l_HI += HI[i];
	}

	l_avg = l_HI / (double)num_of_HI;

	/* calc std */
	l_ss = 0.0;
	for (int i = 0; i < num_of_HI; i++) 
	{
		l_ss += pow_wrapper(fabs(l_avg - HI[i]), 2);
	}

	return l_ss/ pow_wrapper(fabs(temp), 2);
}


void Bcast_best_partial_tree(Dataptr MSA, long best_treelength, int rank, int nprocs, TREESTACK_TREE_NODES* BranchArray, long *tree_root, int* from_rank, MPI_Datatype 	MPI_BRANCH)
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

	*from_rank = index_max;

	/*
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
	*/


    Lvb_bit_length** old_sitestate = (Lvb_bit_length**)alloc(sizeof(Lvb_bit_length*) * MSA->numberofpossiblebranches, "array saving old sitestate");//save original sitestate
    for (int i = 0; i < MSA->numberofpossiblebranches; i++)
    {
        old_sitestate[i] = BranchArray[i].sitestate;
    }


    MPI_Bcast(BranchArray, MSA->numberofpossiblebranches, MPI_BRANCH, index_max, MPI_COMM_WORLD);
    //use MPI_BYTE might have data layout problem(wait improvement)
    //MPI_Bcast(BranchArray, MSA->tree_bytes, MPI_BYTE, index_max, MPI_COMM_WORLD);//广播最屌的树,//广播bit会导致存叶子的顺序不同节点间是不同的
    MPI_Bcast(tree_root, 1, MPI_LONG, index_max, MPI_COMM_WORLD);


    for (int i = 0; i < MSA->numberofpossiblebranches; i++)
    {
        BranchArray[i].sitestate = old_sitestate[i];//restore sitestate
    }
    for (int i = MSA->n; i < MSA->numberofpossiblebranches; i++) 
    {
        BranchArray[i].sitestate[0] = 0U;// make dirty
    }


}

void Root_get_best_treestack(Dataptr MSA, long *best_treelength_local, int root, int rank, int nprocs, TREESTACK* treestack, MPI_Datatype MPI_BRANCH)
{
	long * all_len;

	//printf("\n\n\nlocal best length %ld,  rank %d ----\n\n",*best_treelength_local, rank);
	//if(rank==1)
	//best_treelength=1;

	all_len = (long*)alloc(sizeof(long) * nprocs,"all_len");

	MPI_Allgather(best_treelength_local, 1, MPI_LONG, all_len, 1, MPI_LONG, MPI_COMM_WORLD);
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

	Send_best_treestack_to_root(MSA, rank, 0, index_max, nprocs, treestack,MPI_BRANCH);

	if(rank==0)
		*best_treelength_local=all_len[index_max];

	free(all_len);

}



void  Send_best_treestack_to_root(Dataptr MSA, int rank,int root, int best_rank, int nprocs, TREESTACK * treestack, MPI_Datatype MPI_BRANCH)
{
    //每个process保留使用的seeds，得到的best_length, iteration, 最终_temperature, 
    //（记得每一次都要rinit(seed)），并把最好的一个seed对应的treestack，temperature，
    //给保留，然后交流得到最好的那个，再把最好的那个treestack返回回去
    //注意stack里面brancharray的sitestate是不拷贝进去的
    

	/*
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
	*/
    




    if (best_rank == root)//root already is best
        return;
    else if (best_rank == rank)
    {
        //printf("\n\nbest_rank %d, has stack size %ld and next %ld\n\n", rank, treestack->size, treestack->next);
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


#ifdef test
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

#endif

#ifndef old
void Slave_wait_final_message(MPI_Request *request_message_from_master, MPI_Request *request_handle_send,int *p_n_state_progress, 
	RecvInfoFromMaster *p_data_info_from_master, Parameters *p_rcstruct, int *p_n_number_tried_seed, SendInfoToMaster * p_data_info_to_master, double * critical_t,
	MPI_Datatype mpi_recv_data, MPI_Datatype mpi_data_from_master, int first_commu_to_master)
{
	

	/* Send FINISH message to master */


		//frozen,需要等temp与manage,需要现在发送temp与等待manage；killed则不需要
		if ( *p_n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT)
		{
			if (first_commu_to_master != 1)
			{
				MPI_Wait(request_handle_send, MPI_STATUS_IGNORE);
			}
			if (first_commu_to_master != 1)
			{
				//printf("\nwait,  Slave_anneal for recv_manage\n");
				MPI_Wait(request_message_from_master, MPI_STATUS_IGNORE);
			}

			//这里还需要判断一下上一轮的判定结果是让继续还是杀死，若杀死，则重置为kill的state，不kill则blocking 发送temp，和接受manage（新种子）
			if (p_data_info_from_master->n_is_to_continue == MPI_IS_TO_RESTART_ANNEAL || p_data_info_from_master->n_is_to_continue == MPI_IS_NOT_TO_RESTART)
			{
				*p_n_state_progress = MESSAGE_ANNEAL_KILLED;//还是被kill了，无需再发temp与收manage，因为master已知该slave、被kill了
				
			}
			else
			{
				p_data_info_to_master->n_finish_message = MPI_FINISHED;//frozen,not killed

				MPI_Issend(p_data_info_to_master, 1, mpi_recv_data, MPI_MAIN_PROCESS, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, request_handle_send);
				MPI_Wait(request_handle_send, MPI_STATUS_IGNORE);

				/* need to wait for information if is necessary to run another */
				MPI_Status mpi_status;
				MPI_Recv(p_data_info_from_master, 1, mpi_data_from_master, MPI_MAIN_PROCESS, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, &mpi_status); /* this one waits until the master receive all confirmations */
			}
			
		}

		*request_handle_send = 0;
		*request_message_from_master = 0;


		if (p_data_info_from_master->n_is_to_continue == MPI_IS_NOT_TO_RESTART)
		{
			if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED)
				*p_n_state_progress = MESSAGE_ANNEAL_KILLED_AND_NOT_REPEAT;
			else
				*p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT;
		}
		else 
		{
			p_rcstruct->seed = p_data_info_from_master->n_seed;  /* new seed */
			*critical_t = p_data_info_from_master->Critical_temp; /* critical temperature for new run*/
			if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED)
				*p_n_state_progress = MESSAGE_ANNEAL_KILLED_AND_REPEAT;
			else
				*p_n_state_progress = MESSAGE_ANNEAL_FINISHED_AND_REPEAT;

			*p_n_number_tried_seed = p_data_info_from_master->n_process_tried;  
		}


}


int Slave_after_anneal_once(Dataptr MSA, TREESTACK_TREE_NODES *tree, int n_state_progress, long initroot, TREESTACK *treestack, 
		TREESTACK *best_treestack, int myMPIid, long *l_iterations, long *treelength, long *best_treelength, 
		int n_number_tried_seed_next, Parameters rcstruct)
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
		    printf("\nProcess:%d    Error: can't pop any tree from Treestack.   Rearrangements tried: %ld\n", myMPIid, *l_iterations);
	    }
	    else
	    {
		    CompareTreeToTreestack(MSA, treestack, tree, initroot, LVB_FALSE);
		    //no climbing in the newest version
		    //treelength = deterministic_hillclimb(MSA, bstack_overall, tree, rcstruct, initroot, stdout, &l_iterations, myMPIid, log_progress);
		    //Send and keep the best tree stack this process ever have

			/* compare the best tree with last best*/
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
	    }
    }

    
	//clear working treestack for next iteration, whether it is frozen or killed
	ClearTreestack(treestack);
	
	//repeat, or not
    if (n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT || n_state_progress == MESSAGE_ANNEAL_KILLED_AND_REPEAT)
    {
	    *l_iterations = 0;		/* start iterations from zero */
	    //printf("Process:%d   try seed number process:%d   new seed:%d\n", myMPIid, n_number_tried_seed_next, rcstruct.seed);
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




void write_final_results(Info_record* record, Parameters rcstruct, double overall_time_taken)
{
	char statistics_output[200];
	sprintf(statistics_output, "%s_%d.sta", rcstruct.file_name_in, rcstruct.nruns); /* name of output file for this process */
	char* p = statistics_output;
	while (*p != '\0')
	{
		if (*p == '/'|| *p == '.')
			*p = '-';
		p++;
	}
	//printf("\n%s\n", statistics_output);
	FILE* out = (FILE*)alloc(sizeof(FILE), "alloc statistics output"); 
	out = fopen(statistics_output, "w");


	//calculate average and std and percentage of killing
	int frozen_count = 0;
	long l_len = 0;
	double l_temp = 0.0;
	double l_time = 0.0;
	double l_len_avg, l_temp_avg, l_time_avg, l_len_std, l_temp_std,l_time_std, l_len_ss, l_temp_ss, l_time_ss;

	/* calc average */
	for (int i = 0; i < rcstruct.nruns; i++)
	{
		if (record[i].frozen_or_kill == 1)
		{
			frozen_count++;
			l_len += record[i].result.l_length;
			l_temp += record[i].result.temperature;
			l_time+=record[i].time_consumed;
		}
	}

	l_len_avg = (double)l_len / (double)frozen_count;
	l_temp_avg = (double)l_temp / (double)frozen_count;
	l_time_avg = (double)l_time / (double)frozen_count;

	l_len_ss = 0.0;
	l_temp_ss = 0.0;
	l_time_ss = 0.0;

	l_len_std = 0.0;
	l_temp_std = 0.0;
	l_time_std = 0.0;

	for (int i = 0; i < rcstruct.nruns; i++) 
	{
		if (record[i].frozen_or_kill == 1)
		{
			l_len_ss += pow_wrapper(fabs(l_len_avg - (double)(record[i].result.l_length)), 2);
			l_temp_ss += pow_wrapper(fabs(l_temp_avg - (double)(record[i].result.temperature)), 2);
			l_time_ss += pow_wrapper(fabs(l_time_avg - (double)(record[i].time_consumed)), 2);
		}
		
	}

	l_len_std = sqrt(l_len_ss / (double)frozen_count);
	l_temp_std = sqrt(l_temp_ss / (double)frozen_count);
	l_time_std = sqrt(l_time_ss / (double)frozen_count);





	fprintf(out, "\nCALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS: %d\nCALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS:%d \n",
		CALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS, CALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS);

	fprintf(out, "\nAverage of temperature: %lf \n", l_temp_avg);
	fprintf(out, "\nDeviation of temperature: %lf \n", l_temp_std);
	fprintf(out, "\nAverage of length: %lf \n", l_len_avg);
	fprintf(out, "\nDeviation of length: %lf \n", l_len_std);
	fprintf(out, "\nAverage of time taken(only those frozen): %.2lf \n", l_time_avg);
	fprintf(out, "\nDeviation of time: %lf \n", l_time_std);
	fprintf(out, "\nPercentage of killing: %lf \n", (double)(rcstruct.nruns - frozen_count) / (double)rcstruct.nruns);

	fprintf(out, "\nTime: %.2lf \n",overall_time_taken);

	fprintf(out, "\n--------------------------------------------------------");
	fprintf(out, "\nNo.\tSeed used\tstart temp\tnumber of iterations\tK or F\tLength\ttemperature\tTime\tSlave time\tcommu cost\t \n");
	fprintf(out, "--------------------------------------------------------\n");

	for (int i = 0; i < rcstruct.nruns; i++)
	{
		fprintf(out, "%d\t%d\t%lf\t%d\t\t\t",i+1, record[i].result.n_seed, record[i].result.start_temperature, record[i].result.n_iterations);
		if(record[i].frozen_or_kill==0)
			fprintf(out, "Killed\t\t");
		else
			fprintf(out, "Frozen\t\t");
		fprintf(out, "%ld\t%lf\t%.2lf\t%.2lf\t\t%.2lf\n", 
			record[i].result.l_length, record[i].result.temperature,record[i].time_consumed,record[i].slave_time_consumed,record[i].slave_comm_cost);
	}
	fclose(out);
}

void Slave_send_record_to_Master(int depth, Slave_record* record_slave)
{
	//printf("\n%d, time:%.2lf, no_seed:%d\n", depth,record_slave->anneal_once, record_slave->no_seed);
	if (depth > 1)
		Slave_send_record_to_Master(depth - 1, record_slave->next);

	MPI_Ssend(&(record_slave->anneal_once), 1, MPI_DOUBLE, 0, record_slave->no_seed, MPI_COMM_WORLD);
	MPI_Ssend(&(record_slave->comm_cost), 1, MPI_DOUBLE, 0, record_slave->no_seed, MPI_COMM_WORLD);
	free(record_slave);
}

void Master_recv_record_from_Slave(Info_record* record, int nruns)
{

	MPI_Status mpi_status;
	for (int i = 0; i < nruns; i++)
	{
		
		MPI_Recv(&(record[i].slave_time_consumed), 1, MPI_DOUBLE, MPI_ANY_SOURCE, i+1, MPI_COMM_WORLD, &mpi_status);
		MPI_Recv(&(record[i].slave_comm_cost), 1, MPI_DOUBLE, MPI_ANY_SOURCE, i+1, MPI_COMM_WORLD, &mpi_status);
	}
	
}


void get_temperature_and_control_process_from_other_process(int num_procs, int n_seeds_to_try, Info_record * record, MPI_Datatype mpi_recv_data, MPI_Datatype mpi_send_data)
{

		int nProcessFinished = 0;
		MPI_Request *pHandleTemperatureRecv;
		//MPI_Request *pHandleFinishProcess;
		MPI_Request *pHandleManagementProcess;
		MPI_Status mpi_status;
		IterationTemperature *p_calc_iterations;

		/*删去*pIntFinishedProcess,*/
		int *pIntFinishedProcessChecked,  *pIntProcessControlNumberRunning;//check代表该proc应该是什么状态，
		//有三种可能，MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE，MESSAGE_ANNEAL_STOP_PROCESS，MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN
		int nFlag, i, n_seeds_tried = 0;


		/*collect final result length of -1 means killed*/
		//SendInfoToMaster Final_results[n_seeds_to_try]; 
		//int Final_frozen_or_kill[n_seeds_to_try];//0->killed,1->frozen
		//int idx=0;



		pHandleTemperatureRecv = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for temperature...");
		//pHandleFinishProcess = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for finish process...");
		pHandleManagementProcess = (MPI_Request *) alloc(num_procs * sizeof(MPI_Request), "MPI_request for management process...");
		pIntProcessControlNumberRunning = (int *) alloc(num_procs * sizeof(int), "Buffer with index between process and tried process, but the ones that are running...");
		//pIntFinishedProcess = (int *) alloc(num_procs * sizeof(int), "Buffer for finished processes...");
		pIntFinishedProcessChecked = (int *) alloc(num_procs * sizeof(int), "Buffer for processes states...");
		memset(pIntFinishedProcessChecked, 0, num_procs * sizeof(int)); 	/* control the state of a process, but for all seeds tried */
		memset(pIntProcessControlNumberRunning, 0, num_procs * sizeof(int));		/* has the index number of the process that is running */
		//memset(pIntFinishedProcess, 0, num_procs * sizeof(int));		/* control if a specific process is finished */
		for (i = 0; i < num_procs; i++) { *(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN; } /* all the process start with this state */

		
		/* structure to use sending temperature and number of interactions to master process */
		/*
		int				nItems = 3;
		int          	blocklengths[3] = {3, 1, 1};
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
		*/

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


		int killCall_counts[num_procs - 1];//counts of call to is_possible_to_continue, first CALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS calls won't kill

		/* first get handles for all receivers */
		for (i = 1; i < num_procs; i++) {//temp从0开始用，finish从1开始用，
			MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
			*(pHandleManagementProcess + i) = 0;
			*(pIntProcessControlNumberRunning + i) = i;   /* the index number of ones that are running *///从第0个seed开始算，应该是bug，应该从i开始
			killCall_counts[i - 1] = 0;
		}

		/* alloc memory main calc iterations*/
		p_calc_iterations = get_alloc_main_calc_iterations();
		n_seeds_tried = num_procs - 1;	/* we try always these number of seeds in beginning */


		/*Calculate critical temperature*/
		double sum_sh = 0.0, sum_critical_temp=0.0;//sum of specific heat
		int num_sh = 0;


		printf("\nCALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS: %d\nCALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS:%d \n", 
			CALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS, CALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS);

		printf("\n--------------------------------------------------------");
		printf("\n  No.\tSeed used\tStart temp\tnumber of iterations\tFinal esult\tLength\ttemperature\tRestart\n");
		printf("--------------------------------------------------------\n");

		
		while (1){

			for (i = 1; i < num_procs; i++) {

				/* if is equal to ANNEAL_STOP_PROCESS it is not necessary to do anything else to this process */
				/* because the limit of tried process is reached */
				if (*(pIntFinishedProcessChecked + i) != MESSAGE_ANNEAL_STOP_PROCESS)
				{
					MPI_Test(pHandleTemperatureRecv + i, &nFlag, &mpi_status);
					if (nFlag == 1) {	/* message received */
						//printf("\nProcess:%d    main process getting temperature and length from process:%d   temperature:%-15.8g   length:%ld   iteration:%d\n", mpi_status.MPI_SOURCE,
							//i, (*(p_info_temperature + (i - 1)))->temperature, (*(p_info_temperature + (i - 1)))->l_length,
								//(*(p_info_temperature + (i - 1)))->n_iterations);

						/* all the processes are synchronized by the number of iterations... */
						/* upload the data to the structure */
						/* make calculations to how many standard deviations has */
						/* if it is less than X standard deviations can continue */
						/* add data from this process to structure calc temperature iterations */

						//等待之前的manage到达
						if (*(pHandleManagementProcess + i) != 0)
						{
							MPI_Wait(pHandleManagementProcess + i, MPI_STATUS_IGNORE);							
						}
						
							

						//while 里interval到了，发送的
						if ((*(p_info_temperature + (i - 1)))->n_finish_message == MPI_IS_TO_CONTINUE)
						{
							add_temperature_cal_iterations(p_calc_iterations, *(p_info_temperature + (i - 1)), *(pIntProcessControlNumberRunning + i));
							if (is_possible_to_continue(p_calc_iterations, (*(p_info_temperature + (i - 1)))->temperature,
								(*(p_info_temperature + (i - 1)))->n_iterations, (*(p_info_temperature + (i - 1)))->l_length, num_procs, killCall_counts[i-1])) {
								(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_CONTINUE_ANNEAL;

								killCall_counts[i - 1]++;

								/* send management message */
								/*printf("Process:%d   send management to process:%d\n", 0, i); */
								MPI_Issend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);

								/* launch other receive temperature */
								/* printf("Process:%d   receive temperature from process:%d\n", 0, i); */
								//if ((*(p_info_manage + (i - 1)))->n_is_to_continue != MPI_IS_TO_RESTART_ANNEAL)
								MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);

							}
							else {//kill

								killCall_counts[i - 1] = 0;
								/*receive final tree length*/
								//Final_results[idx] = *(*(p_info_temperature + (i - 1)));
								//Final_frozen_or_kill[idx]=0;
								record[pIntProcessControlNumberRunning[i] - 1].result = *(*(p_info_temperature + (i - 1)));
								record[pIntProcessControlNumberRunning[i] - 1].frozen_or_kill = 0;
								record[pIntProcessControlNumberRunning[i] - 1].end = clock();
								record[pIntProcessControlNumberRunning[i] - 1].time_consumed 
									= ((double)(record[pIntProcessControlNumberRunning[i] - 1].end - record[pIntProcessControlNumberRunning[i] - 1].start)) / CLOCKS_PER_SEC;
								

								printf("%d\t%d\t%lf\t%d\t\t\t", *(pIntProcessControlNumberRunning + i), record[pIntProcessControlNumberRunning[i] - 1].result.n_seed,
									record[pIntProcessControlNumberRunning[i] - 1].result.start_temperature, record[pIntProcessControlNumberRunning[i] - 1].result.n_iterations);
								printf("Killed\t\t");
								printf("%ld\t%lf\t", record[pIntProcessControlNumberRunning[i] - 1].result.l_length, record[pIntProcessControlNumberRunning[i] - 1].result.temperature);
								//idx++;


								//*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE;
								//(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;

								if (n_seeds_tried < n_seeds_to_try) {		/* try another seed */
									n_seeds_tried += 1;									/* try another seed */
									*(pIntProcessControlNumberRunning + i) = n_seeds_tried;
									*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN;
									printf("Yes\n");


									(*(p_info_manage + (i - 1)))->n_seed = get_other_seed_to_run_a_process();
									(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;
									(*(p_info_manage + (i - 1)))->n_process_tried = n_seeds_tried;
									if (num_sh == SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS)
										(*(p_info_manage + (i - 1)))->Critical_temp = sum_critical_temp;
									else
										(*(p_info_manage + (i - 1)))->Critical_temp = -1;

									MPI_Issend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);
									MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);

									//record start time
									record[pIntProcessControlNumberRunning[i] - 1].start = clock();
								}
								else {
									*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS;
									(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_NOT_TO_RESTART;
									//printf("the proc is Process:%d    isn't to restart\n", i);
									printf("No\n");
									
									//Must wait for this slave receiving manage message
									//MPI_Ssend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD);
									MPI_Issend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);
									//printf("\nrank %d, stuck_2\n");

									*(pHandleTemperatureRecv + i) = 0;//无需在接收temp了，因为这个slave不会restart
								}
								nProcessFinished += 1;
							}
						}
						else if ((*(p_info_temperature + (i - 1)))->n_finish_message == MPI_FINISHED)//frozen送来的
						{
							killCall_counts[i - 1] = 0;
							/*receive final tree length*/
							record[pIntProcessControlNumberRunning[i] - 1].result = *(*(p_info_temperature + (i - 1)));
							record[pIntProcessControlNumberRunning[i] - 1].frozen_or_kill = 1;
							record[pIntProcessControlNumberRunning[i] - 1].end = clock();
							record[pIntProcessControlNumberRunning[i] - 1].time_consumed
								= ((double)(record[pIntProcessControlNumberRunning[i] - 1].end - record[pIntProcessControlNumberRunning[i] - 1].start)) / CLOCKS_PER_SEC;


							printf("%d\t%d\t%lf\t%d\t\t\t", *(pIntProcessControlNumberRunning + i), record[pIntProcessControlNumberRunning[i] - 1].result.n_seed,
								record[pIntProcessControlNumberRunning[i] - 1].result.start_temperature, record[pIntProcessControlNumberRunning[i] - 1].result.n_iterations);
							printf("Frozen\t\t");
							printf("%ld\t%lf\t", record[pIntProcessControlNumberRunning[i] - 1].result.l_length, record[pIntProcessControlNumberRunning[i] - 1].result.temperature);
							
							/*Only those frozen will be used to calculate critical temperature*/
							if (num_sh < SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS)
							{
								sum_sh += (*(p_info_temperature + (i - 1)))->peak_sh;
								sum_critical_temp += (*(p_info_temperature + (i - 1)))->Critical_temp;
								num_sh++;
								if (num_sh == SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS)
								{
									sum_sh /= SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS;
									sum_critical_temp /= SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS;
								}
							}


							if (n_seeds_tried < n_seeds_to_try) {		/* try another seed */
								n_seeds_tried += 1;									/* try another seed */
								*(pIntProcessControlNumberRunning + i) = n_seeds_tried;
								*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN;

								(*(p_info_manage + (i - 1)))->n_seed = get_other_seed_to_run_a_process();
								(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_TO_RESTART_ANNEAL;
								(*(p_info_manage + (i - 1)))->n_process_tried = n_seeds_tried;
								if (num_sh == SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS)
									(*(p_info_manage + (i - 1)))->Critical_temp = sum_critical_temp;
								else
									(*(p_info_manage + (i - 1)))->Critical_temp = -1;


								printf("Yes\n");

								//MPI_Ssend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD);
								MPI_Issend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);

								/* launch other wait messages */
								MPI_Irecv(*(p_info_temperature + (i - 1)), 1, mpi_recv_data, i, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, pHandleTemperatureRecv + i);
								//MPI_Irecv(pIntFinishedProcess + i, 1, MPI_INT, i, MPI_TAG_SEND_FINISHED, MPI_COMM_WORLD, pHandleFinishProcess + i);

								//record start time
								record[pIntProcessControlNumberRunning[i] - 1].start = clock();

							}
							else {
								*(pIntFinishedProcessChecked + i) = MESSAGE_ANNEAL_STOP_PROCESS;
								(*(p_info_manage + (i - 1)))->n_is_to_continue = MPI_IS_NOT_TO_RESTART;
								//printf("the proc is Process:%d    isn't to restart\n", i);
								printf("No\n");


								//MPI_Ssend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_SEND_RESTART, MPI_COMM_WORLD);
								MPI_Issend(*(p_info_manage + (i - 1)), 1, mpi_send_data, i, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, pHandleManagementProcess + i);
								//printf("\nrank %d, stuck_2\n");
								/* there's no need other MPI messages */
								*(pHandleTemperatureRecv + i) = 0;
								//*(pHandleFinishProcess + i) = 0;
							}
							nProcessFinished += 1;
						}
					}
				}
				//printf("\nwarn num of seeds used : %d,nProcessFinished: %d \n", n_seeds_tried, nProcessFinished);
				/* test the ones that finish, and all of them need to send this message...  */
				
				
				
			}

			/* is it finished, all of then?  */
			if (nProcessFinished == n_seeds_to_try)
			{
				//Wait all not-to-restart manage
				MPI_Waitall(num_procs - 1, pHandleManagementProcess + 1, MPI_STATUSES_IGNORE);
				/* All processes are finished */
				/* cancel the requests of temperature */
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
		//free(pHandleFinishProcess);			/* all asynchronous messages */
		//free(pIntFinishedProcess);
		free(pIntProcessControlNumberRunning);
		free(pIntFinishedProcessChecked);
	}

#endif
