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

/* ========== parsim.c - tree evaluation ========== */

#include "Lvb.h"

#ifndef NP_Implementation
#ifdef __PPC64__
#include <Popcntll_macro.h>
#endif

long getplen(Dataptr restrict matrix, Branch *barray, Params rcstruct, const long root,
	     long *restrict p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs)
#else
long getplen(Dataptr restrict matrix, Branch *barray, Params rcstruct,
		const long root, long *restrict p_todo_arr,
		long *p_todo_arr_sum_changes, int *p_runs, const long *restrict weights)
#endif
{
    long branch;				/* current branch number */
    long changes = 0;			/* tree length (number of changes) */
    long n_changes_temp;		/* temp variable to count the changes */
    long done = 0;				/* count of branches "done" */
    long i;						/* loop counter */
    long k;						/* current character number */
    long left;					/* current left child number */
    long right;					/* current right child number */
    long todo_cnt = 0;			/* count of branches "to do" */
    long l_end = 0;

//#define	PRINT_PRINTF

    /* calculate state sets and changes where not already known */
    if (matrix->n_threads_getplen > 1){	/* only if is greather than 1 that use the thread version */

#ifdef	PRINT_PRINTF
	printf("start openMP\n");
#endif

		/* get the branches to touch */
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
		#pragma omp parallel num_threads(matrix->n_threads_getplen) private(done, l_end, k, i, left, right, branch, n_changes_temp)
		{
			long j;								/* loop counter */
			long ch;							/* partial changes */
			#ifdef NP_Implementation
			Lvb_bit_length not_u;				// complement of u
			Lvb_bit_length shifted;				// ~u, shifted in partial len calcs
			#endif
			Lvb_bit_length u;					/* for s. set and length calcs */
			Lvb_bit_length x;					/* batch of 8 left state sets */
			Lvb_bit_length y;					/* batch of 8 right state sets */

			done = 0;
			l_end = matrix->n_slice_size_getplen * (omp_get_thread_num() + 1);
			if (matrix->n_threads_getplen == (omp_get_thread_num() + 1)) l_end += matrix->nwords - (matrix->n_slice_size_getplen * matrix->n_threads_getplen);
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
							Lvb_bit_length *restrict l_ssets = barray[left].sset;
							Lvb_bit_length *restrict r_ssets = barray[right].sset;
							for (j = matrix->n_slice_size_getplen * omp_get_thread_num(); j < l_end; j++){
								x = l_ssets[j];
								y = r_ssets[j];

								/* because bootstrap change the weights of the positions it is necessary to look one by one */
								/* MDW: weights gone, still necessary? */
								u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
#ifdef COMPILE_64_BITS
							#if (defined(__x86_64__) && defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
								__asm__ ("popcnt %1, %0" : "=r" (ch) : "0" (u));
							#else
								#ifndef NP_Implementation
								#ifdef __PPC64__
									ch = u;
									LVB_POPCNT_LL(ch);
								#else
								#endif
									ch = __builtin_popcountll(u);
								#endif
								#ifndef NP_Implementation
							#endif
							#endif
#else
								ch = __builtin_popcount(u);
#endif
								ch = LENGTH_WORD - ch;
								u >>= 3;
								#ifndef NP_Implementation
								n_changes_temp += ch;
								#else
								if (rcstruct.bootstraps > 0 && ch > 0){
									not_u = ~u;
									for (k = 0; k < LENGTH_WORD; k++){
										shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
										n_changes_temp += (shifted & 1U) ? (weights[(j >> LENGTH_WORD_BITS_MULTIPLY) + k]) : 0;
									}
								}
								else{
									n_changes_temp += ch;
								}
								#endif
								barray[branch].sset[j] = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));;
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
			for (j = matrix->n_slice_size_getplen * omp_get_thread_num(); j < l_end; j++){
				x = barray[left].sset[j];
				y = barray[right].sset[j];

				u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);

#ifdef COMPILE_64_BITS
			#if (defined(__x86_64__) && defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
				__asm__ ("popcnt %1, %0" : "=r" (ch) : "0" (u));
			#else
				#ifndef NP_Implementation
				#ifdef __PPC64__
					ch = u;
					LVB_POPCNT_LL(ch);
				#else
				#endif
					ch = __builtin_popcountll(u);
				#endif
			#ifndef NP_Implementation
			#endif
			#endif
#else
				ch = __builtin_popcount(u);
#endif

				ch = LENGTH_WORD - ch;
				u >>= 3;
				#ifndef NP_Implementation
				n_changes_temp += ch;
				#else
				if (rcstruct.bootstraps > 0 && ch > 0){
					not_u = ~u;
					for (k = 0; k < LENGTH_WORD; k++){
						shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
						n_changes_temp += (shifted & 1U) ? (weights[(j >> LENGTH_WORD_BITS_MULTIPLY) + k]) : 0;
					}
				}
				else{
					n_changes_temp += ch;
				}
				#endif

				x = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
				y = barray[root].sset[j];

				u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);

#ifdef COMPILE_64_BITS
			#if (defined(__x86_64__) && defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
				__asm__ ("popcnt %1, %0" : "=r" (ch) : "0" (u));
			#else
			#ifndef NP_Implementation
				#ifdef __PPC64__
					ch = u;
					LVB_POPCNT_LL(ch);
				#else
				#endif
					ch = __builtin_popcountll(u);
				#endif
			#ifndef NP_Implementation
			#endif
			#endif
#else
				ch = __builtin_popcount(u);
#endif
				ch = LENGTH_WORD - ch;
				#ifndef NP_Implementation
				n_changes_temp += ch;
			}
			#else
			if (rcstruct.bootstraps > 0 && ch > 0){
					u >>= 3;
					not_u = ~u;

					for (k = 0; k < LENGTH_WORD; k++) {
						shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
						n_changes_temp += (shifted & 1U) ? (weights[(j >> LENGTH_WORD_BITS_MULTIPLY) + k]) : 0;
					}
				}
				else{
					n_changes_temp += ch;
				}
			}
			#endif
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
		/* END of threading code */
    }
	else{	/* code to the orginal version, without threading */

		long ch;					/* partial changes */
		long j;						/* loop counter */

		#ifndef NP_Implementation
		Lvb_bit_length u;				/* for s. set and length calcs */
		Lvb_bit_length x;				/* batch of 8 left state sets */
		Lvb_bit_length y;				/* batch of 8 right state sets */
		#else
		Lvb_bit_length not_u;				/* complement of u */
		Lvb_bit_length shifted;			/* ~u, shifted in partial len calcs */
		Lvb_bit_length u;					/* for s. set and length calcs */
		Lvb_bit_length x;					/* batch of 8 left state sets */
		Lvb_bit_length y;					/* batch of 8 right state sets */
		#endif

		for (i = matrix->n; i < matrix->nbranches; i++) {
			if (barray[i].sset[0] == 0U){
#ifdef	PRINT_PRINTF
			printf("touch: %d   left:%d  right:%d\n", i, barray[i].left, barray[i].right);
#endif
				*(p_todo_arr + todo_cnt++) = i;
				barray[i].changes = 0;
			}
			else{
				changes += barray[i].changes;
			}
		}

		while (done < todo_cnt) {
			for (i = 0; i < todo_cnt; i++) {
				branch = *(p_todo_arr + i);
				if (barray[branch].sset[0] == 0U)	/* "dirty" */
				{
					left = barray[branch].left;
					right = barray[branch].right;
					if (barray[left].sset[0] && barray[right].sset[0])
					{
						Lvb_bit_length *restrict l_ssets = barray[left].sset;
						Lvb_bit_length *restrict r_ssets = barray[right].sset;
						for (j = 0; j < matrix->nwords; j++){
							x = l_ssets[j];
							y = r_ssets[j];
#ifdef	PRINT_PRINTF
	printf("       branch: %d   left = %d\n", branch, left);
	printf("       branch: %d   right = %d\n", branch, right);
#endif
							u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
#ifdef	PRINT_PRINTF
	#ifdef COMPILE_64_BITS
		printf("		u: 0x%016llX    count_bits: %d\n", u, __builtin_popcountll(u));
	#else
		printf("		u: 0x%X    count_bits: %d   x&y: 0x%X\n", u, __builtin_popcount(u), x & y);
	#endif
#endif

#ifdef COMPILE_64_BITS
						#if (defined(__x86_64__) && defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
							__asm__ ("popcnt %1, %0" : "=r" (ch) : "0" (u));
						#else
						#ifndef NP_Implementation
							#ifdef __PPC64__
								ch = u;
								LVB_POPCNT_LL(ch);
							#else
						#endif
								ch = __builtin_popcountll(u);
							#endif
						#ifndef NP_Implementation
						#endif
						#endif
#else
							ch = __builtin_popcount(u);
#endif
							ch = LENGTH_WORD - ch;

							u >>= 3;
							#ifdef NP_Implementation
							if (rcstruct.bootstraps != 0 && ch > 0){
								not_u = ~u;
								ch = 0;
								for (k = 0; k < LENGTH_WORD; k++){
									shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
									ch += (shifted & 1U) ? (weights[(j >> LENGTH_WORD_BITS_MULTIPLY) + k]) : 0;
								}
							}
							#endif

							barray[branch].sset[j] = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
#ifdef	PRINT_PRINTF
	#ifdef COMPILE_64_BITS
		printf("      branch:%d   0x%016llX\n", branch, barray[branch].sset[j]);
	#else
		printf("      branch:%d   0x%X\n", branch, barray[branch].sset[j]);
	#endif
#endif
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

			u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);

#ifdef COMPILE_64_BITS
		#if (defined(__x86_64__) && defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
			__asm__ ("popcnt %1, %0" : "=r" (ch) : "0" (u));
		#else
		#ifndef NP_Implementation
			#ifdef __PPC64__
				ch = u;
				LVB_POPCNT_LL(ch);
			#else
			#endif
				ch = __builtin_popcountll(u);
			#endif
		#ifndef NP_Implementation
		#endif
		#endif
#else
			ch = __builtin_popcount(u);
#endif
			ch = LENGTH_WORD - ch;
			u >>= 3;
			#ifndef NP_Implementation
			changes += ch;
			#else
			if (rcstruct.bootstraps != 0 && ch > 0){
				not_u = ~u;
				for (k = 0; k < LENGTH_WORD; k++){
					shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
					changes += (shifted & 1U) ? (weights[(j >> LENGTH_WORD_BITS_MULTIPLY) + k]) : 0;
				}
			}
			else{
				changes += ch;
			}
			#endif

			x = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
			y = barray[root].sset[j];

			u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);

#ifdef COMPILE_64_BITS
		#if (defined(__x86_64__) && defined(__GNUC__) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
			__asm__ ("popcnt %1, %0" : "=r" (ch) : "0" (u));
		#else
		#ifndef NP_Implementation
			#ifdef __PPC64__
				ch = u;
				LVB_POPCNT_LL(ch);
			#else
			#endif
				ch = __builtin_popcountll(u);
			#endif
		#ifndef NP_Implementation
		#endif
		#endif
#else
			ch = __builtin_popcount(u);
#endif
			ch = LENGTH_WORD - ch;
			#ifndef NP_Implementation
			changes += ch;
			#else
			if (rcstruct.bootstraps != 0 && ch > 0){
				u >>= 3;
				not_u = ~u;
				for (k = 0; k < LENGTH_WORD; k++) {
					shifted = not_u >> (k << NIBBLE_WIDTH_BITS);
					changes += (shifted & 1U) ? (weights[(j >> LENGTH_WORD_BITS_MULTIPLY) + k]) : 0;
				}
			}
			else{
				changes += ch;
			}
			#endif
		}
    }

    lvb_assert(changes > 0);
#ifdef	PRINT_PRINTF
    printf("changes: %d\n", changes);
#endif
 //   exit(1);
    return changes;

} /* end getplen() */
