/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== wrapper.c - LVB to PHYLIP interface ========== */

#include "Lvb.h"

#ifndef NP_Implementation
#ifdef MAP_REDUCE_SINGLE
	#include "InputOptions.h"
#else
	int read_file(char *file_name, int n_file_format, Dataptr p_lvbmat, DataSeqPtr lvbmat_seq);
	void phylip_mat_dims_in_external(char *file_name, int n_file_format, long *species_ptr, long *sites_ptr, int *max_length_name);
#endif
#else
	int read_file(char *file_name, int n_file_format, Dataptr p_lvbmat);
	void phylip_mat_dims_in_external(char *file_name, int n_file_format, long *species_ptr, long *sites_ptr, int *max_length_name);
#endif

	#ifndef NP_Implementation
	int phylip_dna_matrin(char *p_file_name, int n_file_format, Dataptr lvbmat, DataSeqPtr lvbmat_seq)
	#else
	void phylip_dna_matrin(char *p_file_name, int n_file_format, Dataptr lvbmat)
	#endif
	{
		#ifndef NP_Implementation
		int n_error_code = read_file(p_file_name, n_file_format, lvbmat, lvbmat_seq);
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
		#else
		read_file(p_file_name, n_file_format, lvbmat);
		#endif

		/* check number of sequences is in range for LVB */
	    if (lvbmat->n < MIN_N) crash("The data matrix must have at least %ld sequences.", MIN_N);
	    else if (lvbmat->n > MAX_N) crash("The data matrix must have no more than %ld sequences.", MAX_N);
	    /* check number of sites is in range for LVB */
	    else if (lvbmat->m < MIN_M) crash("The data matrix must have at least %ld sites.", MIN_M);
	    else if (lvbmat->m > MAX_M) crash("The data matrix must have no more than %ld sites.", MAX_M);
		#ifndef NP_Implementation
	    return EXIT_SUCCESS;
		#else
		lvb_assert (lvbmat->nsets <= (MAX_N - 3));
		#endif
	} /* end phylip_dna_matrin() */

void phylip_mat_dims_in(char *p_file_name, int n_file_format, long *species_ptr, long *sites_ptr, int *max_length_name){

	phylip_mat_dims_in_external(p_file_name, n_file_format, species_ptr, sites_ptr, max_length_name);
}