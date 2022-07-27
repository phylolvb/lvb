#include "LVB.h"
#include "mpi.h"

/* MPI definitions... */
#define MPI_MAIN_PROCESS	0		/* main process */

#define	MPI_TAG_MATRIX				1
#define	MPI_TAG_NAME_AND_SEQ_DATA		2
#define	MPI_TAG_BINARY_DATA				3
#define MPI_TAG_PARAMS					4
#define MPI_TAG_SEND_TEMP_MASTER		5
#define MPI_TAG_SEND_FINISHED			6
#define MPI_TAG_SEND_RESTART			7
#define MPI_TAG_MANAGEMENT_MASTER		8


#define MPI_FINISHED							0x00
#define MPI_IS_TO_RESTART_ANNEAL				0x01
#define MPI_IS_TO_CONTINUE_ANNEAL				0x02
#define MPI_IS_NOT_TO_RESTART					0x03
#define MPI_IS_TO_CONTINUE						0x04

/* END MPI definitions... */

/* Anneal state *///在anneal里会有int *p_n_state_progress存放下面的属性
#define MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN			0x00
#define MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT				0x01
#define MESSAGE_ANNEAL_FINISHED_AND_REPEAT					0x02
#define MESSAGE_ANNEAL_KILLED								0x03
#define MESSAGE_ANNEAL_KILLED_AND_REPEAT					0x04
#define MESSAGE_ANNEAL_KILLED_AND_NOT_REPEAT				0x05
#define MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE		0x06
#define MESSAGE_ANNEAL_STOP_PROCESS							0x07
#define MESSAGE_BEGIN_CONTROL								0x08
/* Anneal state */

typedef struct
{
        int n_iterations;		/* number of iterations */
        int n_seed;				/* seed for this temperature and iteration */
        long l_length;			/* length of the tree */
        double temperature;		/* temperature */
} SendInfoToMaster;

/* structure to use sending if is to continue and new seed if it is */
typedef struct
{
        int n_seed;				/* new seed to start the process again */
        int n_is_to_continue;	/* if it is to start the process *//*可以等于上面MPI_开头的macro*/
        int n_process_tried;	/* id of the seed tried */
} RecvInfoFromMaster;



int get_other_seed_to_run_a_process();
//void Root_get_best_treestack(Dataptr MSA, long best_length, int root, int rank, int nprocs, TREESTACK* treestack);
//void  Send_best_treestack_to_root(Dataptr MSA, int rank, int root, int best_rank, int nprocs, long best_treelength, TREESTACK* treestack);
void  Send_best_treestack_to_root(Dataptr MSA, int rank, int root, int best_rank, int nprocs, TREESTACK* treestack);
void Root_get_best_treestack(Dataptr MSA, long *best_treelength_local, int root, int rank, int nprocs, TREESTACK* treestack);
void Bcast_best_partial_tree_to_root(Dataptr MSA, long best_treelength, int rank, int nprocs, TREESTACK_TREE_NODES* BranchArray, long* tree_root);
void Slave_interval_reached(MPI_Request *request_handle_send,SendInfoToMaster *p_data_info_to_master, MPI_Datatype mpi_recv_data, MPI_Request *request_message_from_master, RecvInfoFromMaster * p_data_info_from_master, MPI_Datatype mpi_data_from_master, int *p_n_state_progress);

void Slave_wait_final_message(MPI_Request *request_message_from_master, MPI_Request *request_handle_send,int *p_n_state_progress, RecvInfoFromMaster *p_data_info_from_master, Parameters *p_rcstruct, int *p_n_number_tried_seed, SendInfoToMaster * p_data_info_to_master, MPI_Datatype mpi_rec_data, MPI_Datatype mpi_data_from_master);

int Slave_after_anneal_once(Dataptr MSA, TREESTACK_TREE_NODES *tree, int n_state_progress, long * initroot, TREESTACK *treestack, 
		TREESTACK *best_treestack, int myMPIid, long *l_iterations, long *treelength, long *best_treelength, int n_number_tried_seed_next);
void get_temperature_and_control_process_from_other_process(int num_procs, int n_seeds_to_try);

