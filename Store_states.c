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

/* ********** store_states.c - solving functions ********** */

#include "Lvb.h"
#include "Store_states.h"

#define CRC32_POLYNOMIAL 					0xEDB88320
#define MINIMUM_FILE_CHEKPOINT_SIZE			100	/* value in bytes, this is a  */


unsigned long CRC32Value(int i)
{
	int j;
	unsigned long ulCRC = i;
	for (j=8; j>0; j--) {
		if (ulCRC & 1) ulCRC = (ulCRC >> 1)^CRC32_POLYNOMIAL;
		else ulCRC >>= 1;
	}
	return ulCRC;
}

/* calculate the CDC32 checksum */
unsigned long CalculateBlockCRC32(unsigned long ulCount, unsigned char *ucBuffer, unsigned long previousUlCRC)
{
	unsigned long ulCRC = previousUlCRC;
	while (ulCount-- != 0) {
		ulCRC = ((ulCRC >> 8) & 0x00FFFFFFL)^CRC32Value(((int)ulCRC^*ucBuffer++)&0xff);
	}
	return ulCRC;
}


/* print number of bytes of a specific block and the check sum */
void print_information_checkpoint(char *title, int n_bytes_to_write, unsigned long checksum){
	printf("\n%s        Bytes:%d      checksum:%lu\n", title, n_bytes_to_write, checksum);
}


/* test one block of data */
Lvb_bool test_block_data(FILE *fp){

	unsigned long n_bytes_to_read, n_elements_read;
	unsigned char *read_data;
	unsigned long checksum = 0, checksum_read;

	n_elements_read = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp);
//printf("n_bytes_to_read:%lu\n", n_bytes_to_read);
	if (n_elements_read != 1) return LVB_FALSE; 	/* some error with the file, start again */
	checksum = CalculateBlockCRC32(sizeof(n_bytes_to_read), (unsigned char*) &n_bytes_to_read, checksum);
	read_data = (unsigned char *) alloc(n_bytes_to_read, "alloc data to read from checkpoint block");
	n_elements_read = fread(read_data, n_bytes_to_read, 1, fp);
	if (n_elements_read != 1) return LVB_FALSE; 	/* some error with the file, start again */
	checksum = CalculateBlockCRC32(n_bytes_to_read, read_data, checksum);
	n_elements_read = fread(&checksum_read, sizeof(unsigned long), 1, fp);
//print_information_checkpoint("test_block_data", n_bytes_to_read, checksum);
//printf("Check_read:%lu\n", checksum_read);
	free(read_data);
	if (n_elements_read != 1) return LVB_FALSE; 	/* some error with the file, start again */
	if (checksum != checksum_read) return LVB_FALSE; 	/* some error with the file, start again */
	return LVB_TRUE;
}

/* test the consistency of the file, checksums and number of blocks */
/* if OK the process can be continue from a specific state point, otherwise need to start form begining */
/* Several blocks in the file: */
/*		1: int is_process_finished;		*/
/*		2: int number of blocks;		*/
/*		3: uni struture					*/
/*		4: Params struture				*/
/*		5: tree stack struture			*/

/* Each structure/block starts with:			*/
/*		1: unsigned long n_length_bytes_block	*/
/*		2: unsigned short type_block, IDs are defined in store_states.h*/
/*		3: block/structure with data			*/
/*		4: unsigned long CRC32 checksum,		*/
/*		Important: the checksum is calculated with n_length_bytes_block, type_block
 * 			and block/structure data, always, so, each block has
 * 			your own checksum */


/* IMPORTANT, if you change the structure of the file you need to change the tests:
 * 	1) test_lib_checkpoint_file_finished
 * 	2) test_lib_checkpoint_file
 */

Lvb_bool test_consistency_state_file(char *file_name, int myMPIid){

	Lvb_bool is_file_OK = LVB_FALSE, is_file_OK_temp;
	unsigned long n_read_values;
	int n_number_blocks;
	FILE *fp;			/* checkpoint file */

	fp = fopen(file_name, "rb");
	if (fp == NULL) return is_file_OK; 	/* some error with the file, start again */
	fseek(fp, 0L, SEEK_END);
	long int l_file_size = ftell(fp);
	if (l_file_size < MINIMUM_FILE_CHEKPOINT_SIZE){
		printf("\tCheckpoint file to small for MPIid %d\n", myMPIid);
		fclose(fp);
		return is_file_OK; 	/* some error with the file, start again the process */
	}
	/* point to begin plus Flag process finish */
	fseek(fp, 0L + sizeof(int), SEEK_SET);
	n_read_values = fread(&n_number_blocks, sizeof(n_number_blocks), 1, fp);
	lvb_assert(n_read_values == 1);
	/* read block by blck and test the consistency */
	for (int i = 0; i < n_number_blocks; i++){
		is_file_OK_temp = test_block_data(fp);
		if (is_file_OK_temp == LVB_FALSE){
			printf("\tFail consistency of checkpoint file for MPIid %d\n", myMPIid);
			fclose(fp);
			return is_file_OK;		/* some error with the file, start again the process */
		}
	}
	fclose(fp);
	return LVB_TRUE;
}

Lvb_bool is_state_file_exist_and_consistency(int myMPIid){
	char filename[LVB_FNAMSIZE];
	sprintf(filename, "%s_%d", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	return test_consistency_state_file(filename, myMPIid);
}

/* set the FILE cursor to the position of a block, the function will open the FILE pointer */
Lvb_bool point_file_pointer_to_block(FILE *fp, unsigned short type_block){

	int n_number_blocks;
	unsigned long n_bytes_to_read, n_elements_read, n_pos_file;
	unsigned short type_block_read;

	/* point to begin plus Flag process finish */
	fseek(fp, 0L + sizeof(int), SEEK_SET);	// FLAGS: is_process_finished, number_of_bloks,
	n_elements_read = fread(&n_number_blocks, sizeof(n_number_blocks), 1, fp);
	if (n_elements_read != 1) return LVB_FALSE; 	/* some error with the file, start again */
	n_pos_file = sizeof(int) << 1;

	for (unsigned short i = 0; i < n_number_blocks; i++){
		n_elements_read = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp);
		n_elements_read = fread(&type_block_read, sizeof(type_block_read), 1, fp);
		if (type_block_read == type_block){
			fseek(fp, 0L + n_pos_file, SEEK_SET);	// previous position
			return LVB_TRUE;
		}
		n_pos_file += n_bytes_to_read + sizeof(unsigned long) /* this belong to n_bytes_to_read */
						+ sizeof(unsigned long); /* this belong to checksum */
		fseek(fp, 0L + n_pos_file, SEEK_SET);	// FLAGS: is_process_finished, number_of_bloks,
	}
	return LVB_FALSE;
}

FILE *open_file_by_MPIid(int myMPIid, char *p_open_type, Lvb_bool b_temp_file_name){
	char filename[LVB_FNAMSIZE];
	if (b_temp_file_name == LVB_TRUE) sprintf(filename, "%s_%d.temp", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	else sprintf(filename, "%s_%d", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	FILE *fp = fopen(filename, p_open_type);
	lvb_assert(fp != NULL);
	return fp;
}

Lvb_bool is_process_ended(char *file_name){
	unsigned long n_elements_read;
	int is_process_finished;
	FILE *fp = fopen(file_name, "rb");
	if (fp == NULL) return LVB_FALSE; 	/* some error with the file, start again */

	n_elements_read = fread(&is_process_finished, sizeof(is_process_finished), 1, fp);
	if (n_elements_read != 1) return LVB_FALSE; 	/* some error with the file, start again */

	if (is_process_finished == CHECK_POINT_PROCESS_FINISHED) return LVB_TRUE;
	/* other value is not assumed */
	return LVB_FALSE;
}

Lvb_bool is_process_ended_by_MPIid(int myMPIid){
	char filename[LVB_FNAMSIZE];
	sprintf(filename, "%s_%d", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	return is_process_ended(filename);
}
/* each structure has the number of bytes to read in the first position in the file */
/* the last one is a checksum unsigned long, it not summed in the structure */
unsigned long checkpoint_params(FILE *fp, Params *p_rcstruct)
{
	unsigned long n_bytes_to_write = sizeof(long) + 4 * sizeof(int) + LVB_FNAMSIZE + sizeof(unsigned short); /* only yhe name of the file */

#ifndef MAP_REDUCE_SINGLE
	n_bytes_to_write += sizeof(int);	// n_seeds_need_to_try variable
#endif
	unsigned long checksum = 0;
	unsigned short type_block = STATE_BLOCK_PARAMETERS;
	fwrite(&n_bytes_to_write, sizeof(n_bytes_to_write), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_write), (unsigned char *) &n_bytes_to_write, checksum);
	fwrite(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
	fwrite(&p_rcstruct->verbose, sizeof(p_rcstruct->verbose), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->verbose), (unsigned char *) &p_rcstruct->verbose, checksum);
	fwrite(&p_rcstruct->seed, sizeof(p_rcstruct->seed), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->seed), (unsigned char *) &p_rcstruct->seed, checksum);
    fwrite(&p_rcstruct->cooling_schedule, sizeof(p_rcstruct->cooling_schedule), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->cooling_schedule), (unsigned char *) &p_rcstruct->cooling_schedule, checksum);
    fwrite(&p_rcstruct->n_file_format, sizeof(p_rcstruct->n_file_format), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_file_format), (unsigned char *) &p_rcstruct->n_file_format, checksum);
    fwrite(&p_rcstruct->n_processors_available, sizeof(p_rcstruct->n_processors_available), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_processors_available), (unsigned char *) &p_rcstruct->n_processors_available, checksum);
#ifndef MAP_REDUCE_SINGLE
    fwrite(&p_rcstruct->n_seeds_need_to_try, sizeof(p_rcstruct->n_seeds_need_to_try), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_seeds_need_to_try), (unsigned char *) &p_rcstruct->n_seeds_need_to_try, checksum);
#endif
    // don't do these two because can change in time
    //fwrite(&p_rcstruct->n_flag_save_states, sizeof(p_rcstruct->n_flag_save_states), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_flag_save_states), (unsigned char *) &p_rcstruct->n_flag_save_states, checksum);
    //fwrite(&p_rcstruct->n_flag_is_finished_process, sizeof(p_rcstruct->n_flag_is_files_read_state_valid), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_flag_is_files_read_state_valid), (unsigned char *) &p_rcstruct->n_flag_is_files_read_state_valid, checksum);
    fwrite(p_rcstruct->file_name_in, sizeof(p_rcstruct->file_name_in), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->file_name_in), (unsigned char *) p_rcstruct->file_name_in, checksum);
    fwrite(&checksum, sizeof(unsigned long), 1, fp);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(fflush(fp) == 0);
//   print_information_checkpoint("Save data params", n_bytes_to_write, checksum);
    return checksum;
}

unsigned long restore_params(FILE *fp, Params *p_rcstruct)
{
	unsigned long n_bytes_to_write = sizeof(long) + 4 * sizeof(int) + LVB_FNAMSIZE + sizeof(unsigned short), n_bytes_to_read = 0;
#ifndef MAP_REDUCE_SINGLE
	n_bytes_to_write += sizeof(int);	// n_seeds_need_to_try variable
#endif
	unsigned long checksum = 0, checksum_read, n_read_values;
	unsigned short type_block;
	n_read_values = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_read), (unsigned char *) &n_bytes_to_read, checksum);
	n_read_values = fread(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
	n_read_values = fread(&p_rcstruct->verbose, sizeof(p_rcstruct->verbose), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->verbose), (unsigned char *) &p_rcstruct->verbose, checksum);
	n_read_values = fread(&p_rcstruct->seed, sizeof(p_rcstruct->seed), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->seed), (unsigned char *) &p_rcstruct->seed, checksum);
	n_read_values = fread(&p_rcstruct->cooling_schedule, sizeof(p_rcstruct->cooling_schedule), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->cooling_schedule), (unsigned char *) &p_rcstruct->cooling_schedule, checksum);
	n_read_values = fread(&p_rcstruct->n_file_format, sizeof(p_rcstruct->n_file_format), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_file_format), (unsigned char *) &p_rcstruct->n_file_format, checksum);
	n_read_values = fread(&p_rcstruct->n_processors_available, sizeof(p_rcstruct->n_processors_available), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_processors_available), (unsigned char *) &p_rcstruct->n_processors_available, checksum);
#ifndef MAP_REDUCE_SINGLE
	n_read_values = fread(&p_rcstruct->n_seeds_need_to_try, sizeof(p_rcstruct->n_seeds_need_to_try), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_seeds_need_to_try), (unsigned char *) &p_rcstruct->n_seeds_need_to_try, checksum);
#endif
    // don't do these two because can change in time
	//n_read_values = fread(&p_rcstruct->n_flag_save_states, sizeof(p_rcstruct->n_flag_save_states), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_flag_save_states), (unsigned char *) &p_rcstruct->n_flag_save_states, checksum);
	//n_read_values = fread(&p_rcstruct->n_flag_is_finished_process, sizeof(p_rcstruct->n_flag_is_files_read_state_valid), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->n_flag_is_files_read_state_valid), (unsigned char *) &p_rcstruct->n_flag_is_files_read_state_valid, checksum);
	n_read_values = fread(p_rcstruct->file_name_in, sizeof(p_rcstruct->file_name_in), 1, fp); checksum = CalculateBlockCRC32(sizeof(p_rcstruct->file_name_in), (unsigned char *) p_rcstruct->file_name_in, checksum);
	n_read_values = fread(&checksum_read, sizeof(unsigned long), 1, fp);
	lvb_assert(n_read_values == 1);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(n_bytes_to_read == n_bytes_to_write);
    lvb_assert(checksum_read == checksum);
//    print_information_checkpoint("Read data params", n_bytes_to_write, checksum);
    return checksum;
}


unsigned long checkpoint_anneal(FILE *fp, Dataptr restrict matrix, long accepted, Lvb_bool dect, double deltah, long deltalen,
    long failedcnt, long iter, long current_iter, long len, long lenbest, long lendash, double ln_t,
    long t_n, double t0, double pacc, long proposed, double r_lenmin, long rootdash, double t, double grad_geom,
    double grad_linear, Branch *p_current_tree, Lvb_bool b_with_sset_current_tree,
	Branch *p_proposed_tree, Lvb_bool b_with_sset_proposed_tree)
{
	unsigned long n_bytes_to_write = 11 * sizeof(long) + 8 * sizeof(double) + sizeof(unsigned short) + sizeof(Lvb_bool);
	unsigned long checksum = 0;
	unsigned short type_block = STATE_BLOCK_ANNEAL;

	if (b_with_sset_current_tree == LVB_TRUE) n_bytes_to_write += matrix->tree_bytes;
	else n_bytes_to_write += matrix->tree_bytes_without_sset;
	if (b_with_sset_proposed_tree == LVB_TRUE) n_bytes_to_write += matrix->tree_bytes;
	else n_bytes_to_write += matrix->tree_bytes_without_sset;
	fwrite(&n_bytes_to_write, sizeof(n_bytes_to_write), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_write), (unsigned char *) &n_bytes_to_write, checksum);
	fwrite(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
    fwrite(&accepted, sizeof(accepted), 1, fp); checksum = CalculateBlockCRC32(sizeof(accepted), (unsigned char *) &accepted, checksum);
    fwrite(&dect, sizeof(dect), 1, fp); checksum = CalculateBlockCRC32(sizeof(dect), (unsigned char *) &dect, checksum);
    fwrite(&deltah, sizeof(deltah), 1, fp); checksum = CalculateBlockCRC32(sizeof(deltah), (unsigned char *) &deltah, checksum);
    fwrite(&deltalen, sizeof(deltalen), 1, fp); checksum = CalculateBlockCRC32(sizeof(deltalen), (unsigned char *) &deltalen, checksum);
    fwrite(&failedcnt, sizeof(failedcnt), 1, fp); checksum = CalculateBlockCRC32(sizeof(failedcnt), (unsigned char *) &failedcnt, checksum);
    fwrite(&iter, sizeof(iter), 1, fp); checksum = CalculateBlockCRC32(sizeof(iter), (unsigned char *) &iter, checksum);
    fwrite(&current_iter, sizeof(iter), 1, fp); checksum = CalculateBlockCRC32(sizeof(current_iter), (unsigned char *) &current_iter, checksum);
    fwrite(&len, sizeof(len), 1, fp); checksum = CalculateBlockCRC32(sizeof(len), (unsigned char *) &len, checksum);
    fwrite(&lenbest, sizeof(lenbest), 1, fp); checksum = CalculateBlockCRC32(sizeof(lenbest), (unsigned char *) &lenbest, checksum);
    fwrite(&lendash, sizeof(lendash), 1, fp); checksum = CalculateBlockCRC32(sizeof(lendash), (unsigned char *) &lendash, checksum);
    fwrite(&ln_t, sizeof(ln_t), 1, fp); checksum = CalculateBlockCRC32(sizeof(ln_t), (unsigned char *) &ln_t, checksum);
    fwrite(&t_n, sizeof(t_n), 1, fp); checksum = CalculateBlockCRC32(sizeof(t_n), (unsigned char *) &t_n, checksum);
    fwrite(&t0, sizeof(t0), 1, fp); checksum = CalculateBlockCRC32(sizeof(t0), (unsigned char *) &t0, checksum);
    fwrite(&pacc, sizeof(pacc), 1, fp); checksum = CalculateBlockCRC32(sizeof(pacc), (unsigned char *) &pacc, checksum);
    fwrite(&proposed, sizeof(proposed), 1, fp); checksum = CalculateBlockCRC32(sizeof(proposed), (unsigned char *) &proposed, checksum);
    fwrite(&r_lenmin, sizeof(r_lenmin), 1, fp); checksum = CalculateBlockCRC32(sizeof(r_lenmin), (unsigned char *) &r_lenmin, checksum);
    fwrite(&rootdash, sizeof(rootdash), 1, fp); checksum = CalculateBlockCRC32(sizeof(rootdash), (unsigned char *) &rootdash, checksum);
    fwrite(&t, sizeof(t), 1, fp); checksum = CalculateBlockCRC32(sizeof(t), (unsigned char *) &t, checksum);
    fwrite(&grad_geom, sizeof(grad_geom), 1, fp); checksum = CalculateBlockCRC32(sizeof(grad_geom), (unsigned char *) &grad_geom, checksum);
    fwrite(&grad_linear, sizeof(grad_linear), 1, fp); checksum = CalculateBlockCRC32(sizeof(grad_linear), (unsigned char *) &grad_linear, checksum);

    /* copy all space in used in treealloc */
    if (b_with_sset_current_tree == LVB_TRUE){
    	fwrite(p_current_tree, matrix->tree_bytes, 1, fp);
    	checksum = CalculateBlockCRC32(matrix->tree_bytes, (unsigned char *) p_current_tree, checksum);
    }
    else{
    	fwrite(p_current_tree, matrix->tree_bytes_without_sset, 1, fp);
    	checksum = CalculateBlockCRC32(matrix->tree_bytes_without_sset, (unsigned char *) p_current_tree, checksum);
    }

    if (b_with_sset_proposed_tree == LVB_TRUE){
    	fwrite(p_proposed_tree, matrix->tree_bytes, 1, fp);
    	checksum = CalculateBlockCRC32(matrix->tree_bytes, (unsigned char *) p_proposed_tree, checksum);
    }
    else{
    	fwrite(p_proposed_tree, matrix->tree_bytes_without_sset, 1, fp);
    	checksum = CalculateBlockCRC32(matrix->tree_bytes_without_sset, (unsigned char *) p_proposed_tree, checksum);
    }
    fwrite(&checksum, sizeof(unsigned long), 1, fp);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(fflush(fp) == 0);
    //   print_information_checkpoint("Save data params", n_bytes_to_write, checksum);
    return checksum;
}

unsigned long restore_anneal(FILE *fp, Dataptr restrict matrix, long *accepted, Lvb_bool *dect, double *deltah, long *deltalen,
    long *failedcnt, long *iter, long *current_iter, long *len, long *lenbest, long *lendash, double *ln_t,
    long *t_n, double *t0, double *pacc, long *proposed, double *r_lenmin, long *rootdash, double *t, double *grad_geom,
    double *grad_linear, Branch *p_current_tree, Lvb_bool b_with_sset_current_tree,
	Branch *p_proposed_tree, Lvb_bool b_with_sset_proposed_tree)
{
	unsigned long n_bytes_to_write = 11 * sizeof(long) + 8 * sizeof(double) + sizeof(unsigned short) + sizeof(Lvb_bool), n_bytes_to_read = 0;
	unsigned long checksum = 0, checksum_read, n_read_values;
	unsigned short type_block;
	Lvb_bit_length **p_array;

	if (b_with_sset_current_tree == LVB_TRUE || b_with_sset_proposed_tree == LVB_TRUE){
		p_array = (Lvb_bit_length **) alloc(matrix->nbranches * sizeof(Lvb_bit_length *), "alloc array Lvb_bit_length");
	}

	if (b_with_sset_current_tree == LVB_TRUE) n_bytes_to_write += matrix->tree_bytes;
	else n_bytes_to_write += matrix->tree_bytes_without_sset;
	if (b_with_sset_proposed_tree == LVB_TRUE) n_bytes_to_write += matrix->tree_bytes;
	else n_bytes_to_write += matrix->tree_bytes_without_sset;

	n_read_values = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_read), (unsigned char *) &n_bytes_to_read, checksum);
	n_read_values = fread(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
	n_read_values = fread(accepted, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) accepted, checksum);
	n_read_values = fread(dect, sizeof(Lvb_bool), 1, fp); checksum = CalculateBlockCRC32(sizeof(Lvb_bool), (unsigned char *) dect, checksum);
	n_read_values = fread(deltah, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) deltah, checksum);
	n_read_values = fread(deltalen, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) deltalen, checksum);
	n_read_values = fread(failedcnt, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) failedcnt, checksum);
	n_read_values = fread(iter, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) iter, checksum);
	n_read_values = fread(current_iter, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) current_iter, checksum);
	n_read_values = fread(len, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) len, checksum);
	n_read_values = fread(lenbest, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) lenbest, checksum);
	n_read_values = fread(lendash, sizeof(lendash), 1, fp); checksum = CalculateBlockCRC32(sizeof(lendash), (unsigned char *) lendash, checksum);
	n_read_values = fread(ln_t, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) ln_t, checksum);
	n_read_values = fread(t_n, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) t_n, checksum);
	n_read_values = fread(t0, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) t0, checksum);
	n_read_values = fread(pacc, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) pacc, checksum);
	n_read_values = fread(proposed, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) proposed, checksum);
	n_read_values = fread(r_lenmin, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) r_lenmin, checksum);
	n_read_values = fread(rootdash, sizeof(long), 1, fp); checksum = CalculateBlockCRC32(sizeof(long), (unsigned char *) rootdash, checksum);
	n_read_values = fread(t, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) t, checksum);
	n_read_values = fread(grad_geom, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) grad_geom, checksum);
	n_read_values = fread(grad_linear, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) grad_linear, checksum);

	if (b_with_sset_current_tree == LVB_TRUE){
		for (int i = 0; i < matrix->nbranches; i++) *(p_array + i) = p_current_tree[i].sset;
		n_read_values = fread(p_current_tree, matrix->tree_bytes, 1, fp);
		checksum = CalculateBlockCRC32(matrix->tree_bytes, (unsigned char *) p_current_tree, checksum);
		for (int i = 0; i < matrix->nbranches; i++) p_current_tree[i].sset = *(p_array + i);
	}
	else{
		n_read_values = fread(p_current_tree, matrix->tree_bytes_without_sset, 1, fp);
		checksum = CalculateBlockCRC32(matrix->tree_bytes_without_sset, (unsigned char *) p_current_tree, checksum);
	}
	if (b_with_sset_proposed_tree == LVB_TRUE){
		for (int i = 0; i < matrix->nbranches; i++) *(p_array + i) = p_proposed_tree[i].sset;
		n_read_values = fread(p_proposed_tree, matrix->tree_bytes, 1, fp);
		checksum = CalculateBlockCRC32(matrix->tree_bytes, (unsigned char *) p_proposed_tree, checksum);
		for (int i = 0; i < matrix->nbranches; i++) p_proposed_tree[i].sset = *(p_array + i);
	}
	else{
		n_read_values = fread(p_proposed_tree, matrix->tree_bytes_without_sset, 1, fp);
		checksum = CalculateBlockCRC32(matrix->tree_bytes_without_sset, (unsigned char *) p_proposed_tree, checksum);
	}
	n_read_values = fread(&checksum_read, sizeof(unsigned long), 1, fp);

	if (b_with_sset_current_tree == LVB_TRUE || b_with_sset_proposed_tree == LVB_TRUE) free(p_array);
    lvb_assert(n_read_values == 1);
	lvb_assert(ferror(fp) == 0);
	lvb_assert(n_bytes_to_read == n_bytes_to_write);
	lvb_assert(checksum_read == checksum);
	return checksum;
}

Lvb_bool compare_params(Params *p_rcstruct, Params *p_rcstruct_2, Lvb_bool b_test_seed){

	if (p_rcstruct->verbose != p_rcstruct_2->verbose) return LVB_FALSE;
	if (b_test_seed == LVB_TRUE && p_rcstruct->seed != p_rcstruct_2->seed) return LVB_FALSE;

	if (p_rcstruct->cooling_schedule != p_rcstruct_2->cooling_schedule) return LVB_FALSE;
	if (p_rcstruct->n_file_format != p_rcstruct_2->n_file_format) return LVB_FALSE;
	if (p_rcstruct->n_processors_available != p_rcstruct_2->n_processors_available) return LVB_FALSE;
#ifndef MAP_REDUCE_SINGLE
	if (p_rcstruct->n_seeds_need_to_try != p_rcstruct_2->n_seeds_need_to_try) return LVB_FALSE;
#endif
	// don't do these two because can change in time
	//if (p_rcstruct->n_flag_save_read_states != p_rcstruct_2->n_flag_save_read_states) return LVB_FALSE;
	//if (p_rcstruct->n_flag_is_finished_process != p_rcstruct_2->n_flag_is_finished_process) return LVB_FALSE;
	if (memcmp(p_rcstruct->file_name_in, p_rcstruct_2->file_name_in, LVB_FNAMSIZE) != 0) return LVB_FALSE;
	return LVB_TRUE;
}

Lvb_bool is_parameters_are_same_from_state(Params *p_rcstruct, int myMPIid, Lvb_bool b_test_different_seeds){

	Params rcstruct;
	FILE *fp = open_file_by_MPIid(myMPIid, "rb", LVB_FALSE);
	lvb_assert(point_file_pointer_to_block(fp, STATE_BLOCK_PARAMETERS) == LVB_TRUE);

	/* read parameters block and compare with the one that was saved before */
	restore_params(fp, &rcstruct);
	lvb_assert(fclose(fp) == 0);

	/* test different seeds */
    Lvb_bool b_is_same = compare_params(p_rcstruct, &rcstruct, LVB_FALSE);
    if (b_is_same == LVB_TRUE && b_test_different_seeds == LVB_TRUE && rcstruct.seed != p_rcstruct->seed){
    	printf("Warning: The seed for this process %d is different in the checkpoint file read from previous state.\n", myMPIid);
    	printf("         The seed will be replaced by this one :%d\n", rcstruct.seed);
    	p_rcstruct->seed = rcstruct.seed;
    }
    return b_is_same;
}

void save_finish_state_file(Params *p_rcstruct, int myMPIid){
	char filename[LVB_FNAMSIZE], filename_temp[LVB_FNAMSIZE];
	int n_number_blocks, is_process_finished;
	sprintf(filename_temp, "%s_%d.temp", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	sprintf(filename, "%s_%d", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */

	FILE *fp;		/* checkpoint file */

	/* start write */
	fp = fopen(filename_temp, "wb");
	lvb_assert(fp != NULL);
	is_process_finished = CHECK_POINT_PROCESS_FINISHED;
	fwrite(&is_process_finished, sizeof(is_process_finished), 1, fp);
	n_number_blocks = 1;
	fwrite(&n_number_blocks, sizeof(n_number_blocks), 1, fp);
	checkpoint_params(fp, p_rcstruct);
   	lvb_assert(fclose(fp) == 0);

	/* atomic operation
	 *	Change temp name for real name in atomic operation.
	 */
	rename(filename_temp, filename);
}


/* rename temp file name to normal name */
void rename_file_name(int myMPIid){
	char filename[LVB_FNAMSIZE], filename_temp[LVB_FNAMSIZE];
	sprintf(filename_temp, "%s_%d.temp", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	sprintf(filename, "%s_%d", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	/* atomic operation
	 *	Change temp name for real name in atomic operation.
	 */
	rename(filename_temp, filename);
}


/* rename temp file name to normal name */
void remove_checkpoint_file(int myMPIid){
	char filename[LVB_FNAMSIZE];
	sprintf(filename, "%s_%d", CHECKPOINT_FNAM_BASE, myMPIid); /* name of output file for this process */
	/* atomic operation
	 *	Change temp name for real name in atomic operation.
	 */
	remove(filename);
}


void print_memory_hex(char *p_char, int n_size){

	printf("#################################################\n###############################\n");
	int n_change_line = 40;
	for (int i = 0; i < n_size; i ++){
		printf("%x", *(p_char + i) & 0xff);
		if ((i % n_change_line) == 0) printf("\n");
	}
	printf("\n#################################################\n###############################\n");
}


