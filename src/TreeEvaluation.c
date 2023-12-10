/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== TreeEvaluation.c - tree evaluation ========== */

#include "TreeEvaluation.h"

long getplen(Dataptr restrict MSA, TREESTACK_TREE_NODES *BranchArray, Arguments args, const long root,
			 long *restrict p_todo_arr, long *p_todo_arr_sum_changes, int *p_runs)
{
	long branch;						  /* current branch number */
	long changes = 0;					  /* tree length (number of changes) */
	long n_changes_temp;				  /* temp variable to count the changes */
	long done = 0;						  /* count of branches "done" */
	long i;								  /* loop counter */
	long k;								  /* current character number */
	long left;							  /* current left child number */
	long right;							  /* current right child number */
	long number_of_internal_branches = 0; /* count of branches "to do" */
	long l_end = 0;

	/* calculate state sets and changes where not already known */
	if (MSA->n_threads_getplen > 1)
	{ /* only if is greather than 1 that use the thread version */

		/* get the branches to touch */
		memset(p_runs, 0, MSA->n_threads_getplen * (MSA->numberofpossiblebranches - MSA->n) * sizeof(int));
		for (i = MSA->n; i < MSA->numberofpossiblebranches; i++)
		{
			if (BranchArray[i].sitestate[0] == 0U)
			{
				*(p_todo_arr + number_of_internal_branches++) = i;
			}
			else
			{
				changes += BranchArray[i].changes;
				for (k = 0; k < MSA->n_threads_getplen; k++)
				{
					*(p_runs + ((i - MSA->n) * MSA->n_threads_getplen) + k) = 1;
				}
			}
		}

		omp_set_dynamic(0); /* disable dinamic threathing */
#pragma omp parallel num_threads(MSA->n_threads_getplen) private(done, l_end, k, i, left, right, branch, n_changes_temp)
		{
			long j;			  /* loop counter */
			long ch;		  /* partial changes */
			Lvb_bit_length u; /* for s. set and length calcs */
			Lvb_bit_length x; /* batch of 8 left state sets */
			Lvb_bit_length y; /* batch of 8 right state sets */

			done = 0;
			l_end = MSA->n_slice_size_getplen * (omp_get_thread_num() + 1);
			if (MSA->n_threads_getplen == (omp_get_thread_num() + 1))
				l_end += MSA->nwords - (MSA->n_slice_size_getplen * MSA->n_threads_getplen);

			while (done < number_of_internal_branches)
			{
				for (i = 0; i < number_of_internal_branches; i++)
				{
					branch = *(p_todo_arr + i);
					if (*(p_runs + ((branch - MSA->n) * MSA->n_threads_getplen) + omp_get_thread_num()) == 0) /* "dirty" */
					{
						left = BranchArray[branch].left;
						right = BranchArray[branch].right;
						if ((left < MSA->n || *(p_runs + ((left - MSA->n) * MSA->n_threads_getplen) + omp_get_thread_num()) == 1) &&
							(right < MSA->n || *(p_runs + ((right - MSA->n) * MSA->n_threads_getplen) + omp_get_thread_num()) == 1))
						{
							n_changes_temp = 0;
							Lvb_bit_length *restrict l_sitestates = BranchArray[left].sitestate;
							Lvb_bit_length *restrict r_sitestates = BranchArray[right].sitestate;
							for (j = MSA->n_slice_size_getplen * omp_get_thread_num(); j < l_end; j++)
							{
								x = l_sitestates[j];
								y = r_sitestates[j];
								u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
								__asm__("popcnt %1, %0"
										: "=r"(ch)
										: "0"(u));
								ch = LENGTH_WORD - ch;
								u >>= 3;
								n_changes_temp += ch;
								BranchArray[branch].sitestate[j] = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
								;
							}
							*(p_todo_arr_sum_changes + (i * MSA->n_threads_getplen) + omp_get_thread_num()) = n_changes_temp;
							done++;
							*(p_runs + ((branch - MSA->n) * MSA->n_threads_getplen) + omp_get_thread_num()) = 1;
						}
					}
				}
			}

			/* count the changes to the root one */
			n_changes_temp = 0;
			left = BranchArray[root].left;
			right = BranchArray[root].right;
			for (j = MSA->n_slice_size_getplen * omp_get_thread_num(); j < l_end; j++)
			{
				x = BranchArray[left].sitestate[j];
				y = BranchArray[right].sitestate[j];
				u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
				__asm__("popcnt %1, %0"
						: "=r"(ch)
						: "0"(u));

				ch = LENGTH_WORD - ch;
				u >>= 3;
				n_changes_temp += ch;

				x = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
				y = BranchArray[root].sitestate[j];
				u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
				__asm__("popcnt %1, %0"
						: "=r"(ch)
						: "0"(u));
				ch = LENGTH_WORD - ch;
				n_changes_temp += ch;
			}
			*(p_todo_arr_sum_changes + (number_of_internal_branches * MSA->n_threads_getplen) + omp_get_thread_num()) = n_changes_temp;
		}

		/* sum the changes */
		for (i = 0; i < number_of_internal_branches; i++)
		{
			BranchArray[*(p_todo_arr + i)].changes = *(p_todo_arr_sum_changes + (i * MSA->n_threads_getplen));
			for (k = 1; k < MSA->n_threads_getplen; k++)
			{
				BranchArray[*(p_todo_arr + i)].changes += *(p_todo_arr_sum_changes + (i * MSA->n_threads_getplen) + k);
			}
			changes += BranchArray[*(p_todo_arr + i)].changes;
		}
		/* sum the changes to the root */
		for (k = 0; k < MSA->n_threads_getplen; k++)
		{
			changes += *(p_todo_arr_sum_changes + (number_of_internal_branches * MSA->n_threads_getplen) + k);
		}
		/* END of threading code */
	}
	else
	{ /* code to the orginal version, without threading */

		long ch;		  /* partial changes */
		long j;			  /* loop counter */
		Lvb_bit_length u; /* for s. set and length calcs */
		Lvb_bit_length x; /* batch of 8 left state sets */
		Lvb_bit_length y; /* batch of 8 right state sets */

		for (i = MSA->n; i < MSA->numberofpossiblebranches; i++)
		{
			if (BranchArray[i].sitestate[0] == 0U)
			{
				*(p_todo_arr + number_of_internal_branches++) = i;
				BranchArray[i].changes = 0;
			}
			else
			{
				changes += BranchArray[i].changes;
			}
		}

		while (done < number_of_internal_branches)
		{
			for (i = 0; i < number_of_internal_branches; i++)
			{
				branch = *(p_todo_arr + i);
				if (BranchArray[branch].sitestate[0] == 0U) /* "dirty" */
				{
					left = BranchArray[branch].left;
					right = BranchArray[branch].right;
					if (BranchArray[left].sitestate[0] && BranchArray[right].sitestate[0])
					{
						Lvb_bit_length *restrict l_sitestates = BranchArray[left].sitestate;
						Lvb_bit_length *restrict r_sitestates = BranchArray[right].sitestate;
						for (j = 0; j < MSA->nwords; j++)
						{
							x = l_sitestates[j];
							y = r_sitestates[j];
							u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
							__asm__("popcnt %1, %0"
									: "=r"(ch)
									: "0"(u));
							ch = LENGTH_WORD - ch;

							u >>= 3;
							BranchArray[branch].sitestate[j] = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
							BranchArray[branch].changes += ch;
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
		left = BranchArray[root].left;
		right = BranchArray[root].right;
		for (j = 0; j < MSA->nwords; j++)
		{
			x = BranchArray[left].sitestate[j];
			y = BranchArray[right].sitestate[j];
			u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
			__asm__("popcnt %1, %0"
					: "=r"(ch)
					: "0"(u));
			ch = LENGTH_WORD - ch;
			u >>= 3;
			changes += ch;

			x = (x & y) | ((x | y) & ((u + MASK_SEVEN) ^ MASK_EIGHT));
			y = BranchArray[root].sitestate[j];
			u = ((((x & y & MASK_SEVEN) + MASK_SEVEN) | (x & y)) & MASK_EIGHT);
			__asm__("popcnt %1, %0"
					: "=r"(ch)
					: "0"(u));
			ch = LENGTH_WORD - ch;
			changes += ch;
		}
	}

	lvb_assert(changes > 0);
	return changes;

} /* end getplen() */
