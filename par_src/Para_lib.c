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
    //���������С��length��
    //��bcast��ѵģ�
    //���ÿ��proc��������initial���ò�ͬ���������������ʼ���������Ķ�����SA
    //����ѡ��ѵõ����initial������װ��Para_lib��һ�������Ȼ������������anneal���棬������annealǰ���ã����ھ��ȴ���brancharry���ɣ�������ջ֮���ٿ��ǣ�

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
    //�����õĽṹ
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

    //MPI_Bcast(BranchArray, MSA->numberofpossiblebranches, MPI_BRANCH, index_max, MPI_COMM_WORLD);//�㲥��ŵ���
    //use MPI_BYTE might have data layout problem(wait improvement)
    MPI_Bcast(BranchArray, MSA->tree_bytes, MPI_BYTE, index_max, MPI_COMM_WORLD);//�㲥��ŵ���
    MPI_Bcast(tree_root, 1, MPI_LONG, index_max, MPI_COMM_WORLD);//�㲥��ŵ���

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

    


    //����Ӧ��pulltree�������ŵ����������ڵ�tree����pushtree
    //PullTreefromTreestack(MSA, brancharray);//û����



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
        if (all_len[i] < all_len[index_max])//���ص�һ����̵�
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
    //ÿ��process����ʹ�õ�seeds���õ���best_length, iteration, ����_temperature, 
    //���ǵ�ÿһ�ζ�Ҫrinit(seed)����������õ�һ��seed��Ӧ��treestack��temperature��
    //��������Ȼ�����õ���õ��Ǹ����ٰ���õ��Ǹ�treestack���ػ�ȥ
    //ע��stack����brancharray��sitestate�ǲ�������ȥ��
    


    //�����õĽṹ
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
	

        //����cnt���ȶ������ÿռ䣬����objset
        for (int i = 0; i < treestack->next; i++)
        {
            //������stack�����,�������������͵Ļ���ԭ����sitestateҲ�ᱻ���ǣ�Ҫ�޸ģ�����ԭ��anneal�����оͲ���Ҫ��sitestate���͹�ȥ
                //setstate��ÿ��setװ����long���͵�Ҷ�ӵı�ţ���MSA->row�����ţ���
                //sitestate_2ֱ�ӽ����Ĵ����飬�Ѵ�sitestate�Ŀռ���ָ������������brancharrayһ��
                //������stack���p_sitestate���Ǵ����飬�����copy_sitestate(Dataptr restrict MSA, Objset *p_sitestate_1)


            MPI_Issend(&treestack->stack[i].root, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &req_1[i]);//root
            MPI_Issend(treestack->stack[i].tree, MSA->numberofpossiblebranches, MPI_BRANCH, 0, 0, MPI_COMM_WORLD, &req_2[i]);//����,����Ҫ��bitetate��
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
            //��ʼ��cnt��objset
            for (int j = 0; j < MSA->nsets; j++)
            {
                MPI_Issend(&treestack->stack[i].p_sitestate[j].cnt, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &req_3[j]);//��߿��Խṹ�廯һ������ֻ��cnt����*set
                //�ȴ����洴��set�õĿռ�
            }
            MPI_Waitall(MSA->nsets, req_3, MPI_STATUSES_IGNORE);

            for (int j = 0; j < MSA->nsets; j++)
            {
                MPI_Issend(treestack->stack[i].p_sitestate[j].set, treestack->stack[i].p_sitestate[j].cnt, MPI_LONG, 0, 0, MPI_COMM_WORLD,&req_3[j]);//��߿��Խṹ�廯һ������ֻ��cnt����*set
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
        //�ռ䲻���ʹ��ռ�

        //ֱ�ӷ���һ�ѿռ�,nextָ��Ŀռ䲻һ�����ã�����size==next��Ҫincrease��size���ʾpractically�м������Ŀռ��ѱ�����
        if (treestack->next > treestack->size)
        {
            treestack->stack = (TREESTACK_TREES*)realloc(treestack->stack, (treestack->next) * sizeof(TREESTACK_TREES));
            if (treestack->stack == NULL) {
                crash("out of memory: cannot increase allocation for best tree stack to %ld elements", treestack->size);
            }
            
            /* allocate space within stack */
            for (int i= treestack->size; i < treestack->next; i++)
            {
 
                treestack->stack[i].tree = treealloc(MSA, LVB_FALSE);//LVB_FALSE������sitestate��ȥ
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

        //����
        for (int i = 0; i < treestack->next; i++)
        {


            MPI_Irecv(&treestack->stack[i].root, 1, MPI_LONG, best_rank, 0, MPI_COMM_WORLD, &req_1[i]);
            MPI_Irecv(treestack->stack[i].tree, MSA->numberofpossiblebranches, MPI_BRANCH, best_rank, 0, MPI_COMM_WORLD, &req_2[i]);//����,��ʵ����Ҫ��bitetate���������ͣ���ֹ�ڵ㻹Ҫ�����Ƶ��м�ڵ�

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
                MPI_Irecv(&treestack->stack[i].p_sitestate[j].cnt, 1, MPI_LONG, best_rank, 0, MPI_COMM_WORLD, &req_3[j]);//��߿��Խṹ�廯һ������ֻ��cnt����*set

            }
            MPI_Waitall(MSA->nsets, req_3, MPI_STATUSES_IGNORE);

            for (int j = 0; j < MSA->nsets; j++)
            {

               //�ȴ����洴��set�õĿռ�
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


void interval_reached()
{
	/* send temperature to the master process*/
					if (request_handle_send != 0) { MPI_Wait(&request_handle_send, MPI_STATUS_IGNORE); }
					p_data_info_to_master->n_iterations = *current_iter;
					p_data_info_to_master->n_seed = p_rcstruct->seed;
					p_data_info_to_master->l_length = lenbest;
					p_data_info_to_master->temperature = t;
					/* printf("Process:%d   send temperature:%.3g   iterations:%ld\n", myMPIid, t, *current_iter); */
					MPI_Isend(p_data_info_to_master, 1, mpi_recv_data, MPI_MAIN_PROCESS, MPI_TAG_SEND_TEMP_MASTER, MPI_COMM_WORLD, &request_handle_send);

					/* now get the message to continue or not, but need in second iteration... */
					if (request_message_from_master != 0) {
						MPI_Wait(&request_message_from_master, MPI_STATUS_IGNORE);
						MPI_Test(&request_message_from_master, &nFlag, &mpi_status);
						if (nFlag == 0) { printf("ERROR, mpi waiting is not working File:%s  Line:%d\n", __FILE__, __LINE__); }
						if (nFlag == 1){
							if (p_data_info_from_master->n_is_to_continue == MPI_IS_TO_RESTART_ANNEAL){	/*it's there and need to restart*/
								MPI_Cancel(&request_handle_send);
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
					MPI_Irecv(p_data_info_from_master, 1, mpi_data_from_master, MPI_MAIN_PROCESS, MPI_TAG_MANAGEMENT_MASTER, MPI_COMM_WORLD, &request_message_from_master);
				}
}



void wait_final_message()
{
	if (request_message_from_master != 0) MPI_Cancel(&request_message_from_master);
	    if (request_handle_send != 0) MPI_Cancel(&request_handle_send);

	    /* Send FINISH message to master */
	    if (*p_n_state_progress == MESSAGE_ANNEAL_KILLED || *p_n_state_progress == MESSAGE_ANNEAL_FINISHED_AND_REPEAT){	/* it's necessary to send this message */
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
}
