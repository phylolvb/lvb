/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
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

/**********

=head1 NAME

parsim.c - tree evaluation

=cut

**********/

#include "lvb.h"

long getplen(Dataptr restrict matrix, Branch *barray,
		const long root, const long *restrict weights, long *p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs)
{
	unsigned current_ss;		/* current state set */

    long branch;			/* current branch number */
    long changes = 0;			/* tree length (number of changes) */
    long n_changes_temp;	/* temp variable to count the changes */
    long done = 0;			/* count of branches "done" */
	long i;					/* loop counter */
	long k;					/* current character number */
    long left;				/* current left child number */
    long right;				/* current right child number */
    long todo_cnt = 0;			/* count of branches "to do" */
    long l_end = 0;
//#define PRINT_PRINTF	// print extra information

    /* calculate state sets and changes where not already known */
#ifdef COMPILE_OPEN_MP
#ifdef	PRINT_PRINTF
	printf("start openMP\n");
#endif

	memset(p_runs, 0, matrix->n_threads_getplen * (matrix->nbranches - matrix->n) * sizeof(int));
    for (i = matrix->n; i < matrix->nbranches; i++) {
		if (barray[i].sset[0] == 0U){
#ifdef	PRINT_PRINTF
			printf("touch: %d   left:%d  right:%d\n", i, barray[i].left, barray[i].right);
#endif
			*(p_todo_arr + todo_cnt++) = i;
		}
		else{
			changes += barray[i].changes;
			for (k = 0; k < matrix->n_threads_getplen; k++){
				*(p_runs + ((i - matrix->n) * matrix->n_threads_getplen) + k) = 1;
			}
		}
	}

    omp_set_dynamic(0);	  /* disable dinamic threathing */
	#pragma omp parallel num_threads(matrix->n_threads_getplen) private(done, l_end, k, i, left, right, branch, n_changes_temp, current_ss)
	{
		done = 0;
		l_end = matrix->n_slice_size_getplen * (omp_get_thread_num() + 1);
		if (matrix->n_threads_getplen == (omp_get_thread_num() + 1)) l_end += matrix->m - (matrix->n_slice_size_getplen * matrix->n_threads_getplen);
#ifdef	PRINT_PRINTF
		printf("1 : Thread# %d: begin = %d    l_end: %d\n", omp_get_thread_num(), matrix->n_slice_size_getplen * omp_get_thread_num(), l_end);
#endif
		while (done < todo_cnt) {
			for (i = 0; i < todo_cnt; i++) {
				branch = *(p_todo_arr + i);
#ifdef	PRINT_PRINTF
				printf("1 : Thread# %d: try to make branch: %d     is_possible_to_run: %d\n", omp_get_thread_num(), branch,
						*(p_runs + ((branch - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()));
#endif
				if (*(p_runs + ((branch - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()) == 0) /* "dirty" */
				{
					left = barray[branch].left;
					right = barray[branch].right;
#ifdef	PRINT_PRINTF
					printf("1 : Thread# %d: left = %d  is_possible_to_run: %d\n", omp_get_thread_num(), left, *(p_runs + ((left - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()));
					printf("1 : Thread# %d: right = %d  is_possible_to_run: %d\n", omp_get_thread_num(), right, *(p_runs + ((right - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()));
#endif
					if ((left < matrix->n || *(p_runs + ((left - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()) == 1) &&
							(right < matrix->n || *(p_runs + ((right - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()) == 1))
					{
#ifdef	PRINT_PRINTF
						printf("1 : Thread# %d: make branch: %d\n", omp_get_thread_num(), branch);
#endif
						n_changes_temp = 0;
						for (k = matrix->n_slice_size_getplen * omp_get_thread_num();k < l_end; k++){
							barray[branch].sset[k] = barray[left].sset[k] & barray[right].sset[k];
							if (barray[branch].sset[k] == 0U){
								barray[branch].sset[k] = barray[left].sset[k] | barray[right].sset[k];
								n_changes_temp += weights[k];
							}
						}
						*(p_todo_arr_sum_changes + (i * matrix->n_threads_getplen) + omp_get_thread_num()) = n_changes_temp;

						done++;
						*(p_runs + ((branch - matrix->n) * matrix->n_threads_getplen) + omp_get_thread_num()) = 1;
					}
				}
			}
		}

		/* count the changes to the root one */
		n_changes_temp = 0;
		left = barray[root].left;
		right = barray[root].right;
		for (k = matrix->n_slice_size_getplen * omp_get_thread_num();k < l_end; k++) {
			current_ss = barray[left].sset[k] & barray[right].sset[k];
			if (current_ss == 0U){
				current_ss = barray[left].sset[k] | barray[right].sset[k];
				n_changes_temp += weights[k];
			}
			if ((current_ss & barray[root].sset[k]) == 0U) n_changes_temp += weights[k];
		}
		*(p_todo_arr_sum_changes + (todo_cnt * matrix->n_threads_getplen) + omp_get_thread_num()) = n_changes_temp;
	}


	/* sum the changes */
    for (i = 0; i < todo_cnt; i++) {
    	barray[*(p_todo_arr + i)].changes = *(p_todo_arr_sum_changes + (i * matrix->n_threads_getplen));
    	for(k = 1; k < matrix->n_threads_getplen; k++){
    		barray[*(p_todo_arr + i)].changes += *(p_todo_arr_sum_changes + (i * matrix->n_threads_getplen) + k);
    	}
    	changes += barray[*(p_todo_arr + i)].changes;
    }
    /* sum the changes to the root */
	for(k = 0; k < matrix->n_threads_getplen; k++){
		changes += *(p_todo_arr_sum_changes + (todo_cnt * matrix->n_threads_getplen) + k);
	}

#else

	long ch;					/* partial changes */
	long j;						/* loop counter */
	uint32_t not_u;				/* complement of u */
	uint32_t shifted;			/* ~u, shifted in partial len calcs */
	uint32_t u;					/* for s. set and length calcs */
	uint32_t x;					/* batch of 8 left state sets */
	uint32_t y;					/* batch of 8 right state sets */
	uint32_t z;

    for (i = matrix->n; i < matrix->nbranches; i++) {
		if (barray[i].sset[0] == 0U){
			*(p_todo_arr + todo_cnt++) = i;
//			printf("touch:%d   left:%d  right:%d\n", i, barray[i].left, barray[i].right);
			barray[i].changes = 0;
		}
		else{
			changes += barray[i].changes;
		}
	}

//    printf("start\n");
    while (done < todo_cnt) {
		for (i = 0; i < todo_cnt; i++) {
			branch = *(p_todo_arr + i);
//			printf("touch:%d\n", branch);
			if (barray[branch].sset[0] == 0U)	/* "dirty" */
			{
				left = barray[branch].left;
				right = barray[branch].right;
//				printf("touch:%d   left:%d  right:%d\n", branch, left, right);
				if (barray[left].sset[0] && barray[right].sset[0])
				{
					uint32_t *restrict l_ssets = barray[left].sset;
					uint32_t *restrict r_ssets = barray[right].sset;
					for (j = 0; j < matrix->nwords; j++){
						x = l_ssets[j];
						y = r_ssets[j];
						u = ((((x & y & 0x77777777U) + 0x77777777U) | (x & y)) & 0x88888888U) >> 3;
						z = (x & y) | ((x | y) & ((u + 0x77777777U) ^ 0x88888888U));
						not_u = ~u;
						ch = 0;
						for (k = 0; k < 8; k++){
							/*shifted = not_u >> k * NIBBLE_WIDTH;
							ch += (shifted & 1U) ? (weights[j * 8 + k]) : 0;*/
							shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
							ch += (shifted & 1U) ? (weights[(j >> 3) + k]) : 0;
						}
//						printf("before branch:8  pos:0  memory:%1x   pos:%d    branch:7  memory:0x%1x\n", barray[8].sset[0], j, barray[branch].sset[j]);
						barray[branch].sset[j] = z;
//						printf("before branch:8  pos:0  memory:%1x   pos:%d   branch:7  memory:0x%1x\n", barray[8].sset[0], j, barray[branch].sset[j]);
						barray[branch].changes += ch;
						changes += ch;
					}
					done++;
				}
			}
		}
	}

    /* root: add length for root branch structure, and also for true root which
        * lies outside the LVB tree data structure; all without altering the
        * "root" struct statesets (since these represent actual data for the
        * leaf) */
	left = barray[root].left;
	right = barray[root].right;
	for (j = 0; j < matrix->nwords; j++){
		x = barray[left].sset[j];
		y = barray[right].sset[j];
		u = ((((x & y & 0x77777777U) + 0x77777777U) | (x & y)) & 0x88888888U) >> 3;
		z = (x & y) | ((x | y) & ((u + 0x77777777U) ^ 0x88888888U));
		not_u = ~u;

		for (k = 0; k < 8; k++){
			shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
			changes += (shifted & 1U) ? (weights[(j >> 3) + k]) : 0;
			/*shifted = not_u >> k * NIBBLE_WIDTH;
			changes += (shifted & 1U) * (weights[j * 8 + k]);*/
		}

		x = z;
		y = barray[root].sset[j];
		u = ((((x & y & 0x77777777U) + 0x77777777U) | (x & y)) & 0x88888888U) >> 3;
		z = (x & y) | ((x | y) & ((u + 0x77777777U) ^ 0x88888888U));
		not_u = ~u;

		for (k = 0; k < 8; k++) {
			shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
			changes += (shifted & 1U) ? (weights[(j >> 3) + k]) : 0;
			/*shifted = not_u >> k * NIBBLE_WIDTH;
			changes += (shifted & 1U) ? (weights[j * 8 + k]) : 0;*/
		}
	}

#endif

    lvb_assert(changes > 0);
#ifdef	PRINT_PRINTF
    printf("changes: %d    finish\n", changes);
#endif
 //   exit(1);
    return changes;

} /* end getplen() */
