#include "LVB.h"
#include "mpi.h"


#ifndef STRUCT_PARA
#define STRUCT_PARA
#include "Clock.h"

/* test if is possible to continue */
#define CALC_ITERATION_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS 3
#define CALC_ITERATION_NUMBER_STD_TO_RESTART_PROCESS 1


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

//Master 里面 pIntFinishedProcess的值（初始版本只会是MPI_FINISHED）；也是(*(p_info_manage + (i - 1)))->n_is_to_continue的值有#define MPI_IS_TO_RESTART_ANNEAL				0x01
//#define MPI_IS_TO_CONTINUE_ANNEAL， #define MPI_IS_NOT_TO_RESTART	三种可能性
#define MPI_FINISHED							0x00
#define MPI_IS_TO_RESTART_ANNEAL				0x01
#define MPI_IS_TO_CONTINUE_ANNEAL				0x02
#define MPI_IS_NOT_TO_RESTART					0x03
#define MPI_IS_TO_CONTINUE						0x04

/* END MPI definitions... */

/* Anneal state *///在anneal里会有int *p_n_state_progress，以及pIntFinishedProcessChecked数组存放下面的属性
#define MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN			0x00//正在跑，或者master发了新seed，还不知道slave就没接受到
#define MESSAGE_ANNEAL_FINISHED_AND_NOT_REPEAT				0x01
#define MESSAGE_ANNEAL_FINISHED_AND_REPEAT					0x02
#define MESSAGE_ANNEAL_KILLED								0x03
#define MESSAGE_ANNEAL_KILLED_AND_REPEAT					0x04
#define MESSAGE_ANNEAL_KILLED_AND_NOT_REPEAT				0x05
#define MESSAGE_ANNEAL_STOP_PROCESS_WAIT_FINAL_MESSAGE		0x06//wait final message指在master kill这个proc了，然后等待slave发送FinishMessage了，
//finalmessage特指这个finishmessage，若还有种子，状态就会变成MESSAGE_ANNEAL_IS_RUNNING_OR_WAIT_TO_RUN，反之变为MESSAGE_ANNEAL_STOP_PROCESS
#define MESSAGE_ANNEAL_STOP_PROCESS							0x07//种子已经用完
#define MESSAGE_BEGIN_CONTROL								0x08
/* Anneal state */


/*Parallel selection*/
#define	PARALLEL_STATIC_AMI				0
#define	PARALLEL_CLUSTER	1
#define	PARALLEL_DYNAMIC_AMI_KILL			2
#define PARALLEL_DYNAMIC_AMI_PHASE_TRANSITION					3
#define PARALLEL_DYNAMIC_AMI_KILL_PHASE_TRANSITION 4

/*Control start from critical temperature*/
#define INTERVAL_FOR_SPECIFIC_HEAT 10
#define SPECIFIC_HEAT_ONLY_RELEASE_AFTER_NUMBER_CHUNCHS 10

typedef struct
{
        int n_iterations;		/* number of iterations */
        int n_seed;				/* seed for this temperature and iteration */
        int n_finish_message;  //MPI_FINISHED 指代frozen， MPI_IS_TO_CONTINUE指代单纯在每个interval到了，送的
        long l_length;			/* length of the tree */
        double temperature;		/* temperature */
        double start_temperature;
        double Critical_temp;
        double peak_sh;//peak specific heat
} SendInfoToMaster;

/* structure to use sending if is to continue and new seed if it is */
typedef struct
{
        int n_seed;				/* new seed to start the process again */
        int n_is_to_continue;	/* if it is to start the process *//*可以等于上面MPI_开头的macro*/
        int n_process_tried;	/* id of the seed tried */
        double Critical_temp;
} RecvInfoFromMaster;

/* structures for calculation of averages temperatures and std */
typedef struct IndividualTemperature
{
    double d_temperature;	/* temperature */
    int n_try_process;		/* number of process tried, is sequential...*/
    int n_seed;				/* seed for this temperature and iteration */
                            /* the seed is here because it's easier to perform the algorithm */
    long l_length;			/* length of the tree */
    struct IndividualTemperature* p_next_temperature; /* next memory structure */
}IndividualTemperature;

typedef struct IterationTemperature
{
    int n_iteration;
    IndividualTemperature* p_temperature;
    struct IterationTemperature* p_next_iteration;
}IterationTemperature;


typedef struct Info_record
{
    int frozen_or_kill; //0 ->killed, 1->frozen
    clock_t start;
    clock_t end;
    double time_consumed;//time consumed from Master's perspective
    double slave_time_consumed;//time consumed from Slave's perspective, which is basically same as time_consumed
    double slave_comm_cost;
    int proc; //which process/rank deal with
    SendInfoToMaster result;
}Info_record;

typedef struct Slave_record
{
    double anneal_once;//time  used for anneal once
    double comm_cost;//communication cost of anneal once
    int no_seed;//number of seed used in this anneal
    struct Slave_record* next;
}Slave_record;//Slaves use it to store info during running, and finally after all seeds are used, send these records to 

int get_other_seed_to_run_a_process();

void Root_get_best_treestack(Dataptr MSA, long *best_treelength_local, int root, int rank, int nprocs, TREESTACK* treestack, MPI_Datatype MPI_BRANCH);
void Send_best_treestack_to_root(Dataptr MSA, int rank, int root, int best_rank, int nprocs, TREESTACK* treestack, MPI_Datatype MPI_BRANCH);

void Bcast_best_partial_tree(Dataptr MSA, long best_treelength, int rank, int nprocs, TREESTACK_TREE_NODES* BranchArray, long* tree_root, int *from_rank, MPI_Datatype 	MPI_BRANCH);

void Create_MPI_Datatype(MPI_Datatype* MPI_BRANCH, MPI_Datatype* MPI_SLAVEtoMASTER, MPI_Datatype* MPI_MASTERtoSLAVE);

double Calculate_specific_heat(double* HI, int num_of_HI, double temp);

void Slave_wait_final_message(MPI_Request* request_message_from_master, MPI_Request* request_handle_send, int* p_n_state_progress,
    RecvInfoFromMaster* p_data_info_from_master, Parameters* p_rcstruct, int* p_n_number_tried_seed, SendInfoToMaster* p_data_info_to_master, double * critical_t,
    MPI_Datatype mpi_recv_data, MPI_Datatype mpi_data_from_master, int first_commu_to_master);
int Slave_after_anneal_once(Dataptr MSA, TREESTACK_TREE_NODES* tree, int n_state_progress, long initroot, TREESTACK* treestack,
    TREESTACK* best_treestack, int myMPIid, long* l_iterations, long* treelength, long* best_treelength,
    int n_number_tried_seed_next, Parameters rcstruct);

void get_temperature_and_control_process_from_other_process(int num_procs, int n_seeds_to_try, Info_record* record, MPI_Datatype mpi_recv_data, MPI_Datatype mpi_send_data);

void write_final_results(Info_record* record, Parameters rcstruct, double overall_time_taken);

void Slave_send_record_to_Master(int depth, Slave_record* record_slave);
void Master_recv_record_from_Slave(Info_record* record, int nruns);
#endif