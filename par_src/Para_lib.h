#include "LVB.h"

int get_other_seed_to_run_a_process();
//void Root_get_best_treestack(Dataptr MSA, long best_length, int root, int rank, int nprocs, TREESTACK* treestack);
//void  Send_best_treestack_to_root(Dataptr MSA, int rank, int root, int best_rank, int nprocs, long best_treelength, TREESTACK* treestack);
void  Send_best_treestack_to_root(Dataptr MSA, int rank, int root, int best_rank, int nprocs, TREESTACK* treestack);
void Root_get_best_treestack(Dataptr MSA, long best_treelength_local, int root, int rank, int nprocs, TREESTACK* treestack);
void Bcast_best_partial_tree_to_root(Dataptr MSA, long best_treelength, int rank, int nprocs, TREESTACK_TREE_NODES* BranchArray, long* tree_root);
