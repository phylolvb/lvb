#ifdef LVB_PARALLEL_SEARCH

/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Chang Sik Kim,
Maximilian Strobl and Martyn Winn
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


/* ********** store_states.h - interface for store_states.c ********** */

#include "lvb.h"

#define STATE_BLOCK_UNI						1
#define STATE_BLOCK_PARAMETERS				2
#define STATE_BLOCK_TREESTACK				3
#define STATE_BLOCK_ANNEAL					4
#define STATE_BLOCK_ONLY_TESTS_FAIL			50

#ifndef LVB_STORE_STATES_H
#define LVB_STORE_STATES_H


unsigned long CalculateBlockCRC32(unsigned long ulCount, unsigned char *ucBuffer, unsigned long previousUlCRC);
unsigned long checkpoint_params(FILE *fp, Params *p_rcstruct);
Lvb_bool compare_params(Params *p_rcstruct, Params *p_rcstruct_2, Lvb_bool b_test_seed);
Lvb_bool is_process_ended(char *file_name);
Lvb_bool is_process_ended_by_MPIid(int myMPIid);
Lvb_bool is_parameters_are_same_from_state(Params *p_rcstruct, int myMPIid, Lvb_bool b_test_different_seeds);
Lvb_bool is_state_file_exist_and_consistency(int myMPIid);
Lvb_bool point_file_pointer_to_block(FILE *fp, unsigned short type_block);
FILE *open_file_by_MPIid(int myMPIid, char *p_open_type, Lvb_bool b_temp_file_name);
void print_information_checkpoint(char *title, int n_bytes_to_write, unsigned long checksum);
void print_memory_hex(char *p_char, int n_size);
unsigned long restore_params(FILE *fp, Params *p_rcstruct);
void rename_file_name(int myMPIid);
void remove_checkpoint_file(int myMPIid);
void save_finish_state_file(Params *p_rcstruct, int myMPIid);
Lvb_bool test_block_data(FILE *fp);
Lvb_bool test_consistency_state_file(char *file_name, int myMPIid);


unsigned long checkpoint_anneal(FILE *fp, Dataptr matrix, long accepted, Lvb_bool dect, double deltah, long deltalen,
    long failedcnt, long iter, long current_iter, long len, long lenbest, long lendash, double ln_t,
    long t_n, double t0, double pacc, long proposed, double r_lenmin, long rootdash, double t, double grad_geom,
    double grad_linear, Branch *p_current_tree, Lvb_bool b_with_sset_current_tree,
	Branch *p_proposed_tree, Lvb_bool b_with_sset_proposed_tree);

unsigned long restore_anneal(FILE *fp, Dataptr matrix, long *accepted, Lvb_bool *dect, double *deltah, long *deltalen,
    long *failedcnt, long *iter, long *current_iter, long *len, long *lenbest, long *lendash, double *ln_t,
    long *t_n, double *t0, double *pacc, long *proposed, double *r_lenmin, long *rootdash, double *t, double *grad_geom,
    double *grad_linear, Branch *p_current_tree, Lvb_bool b_with_sset_current_tree,
	Branch *p_proposed_tree, Lvb_bool b_with_sset_proposed_tree);

#endif /* LVB_STORE_STATES_H */

#endif