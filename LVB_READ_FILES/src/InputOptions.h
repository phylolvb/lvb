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

/* ========== ReadFile.h - interface for ReadFile.cpp ========== */

#ifndef READFILE_H_
#define READFILE_H_

#include "MSAInput.h"
#include "../../DataStructure.h"
#include <stdio.h>
#include <string.h>
#include <getopt.h>
using namespace std;

//#define NP_Implementation
// #define MPI_Implementation

#ifdef NP_Implementation

extern "C" void read_file(char *file_name, int n_file_type, Dataptr p_lvbmat);
extern "C" void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);
extern "C" void read_parameters(Params *prms, int argc, char **argv);


void read_file(char *file_name, int n_file_type, DataStructure *p_lvbmat);
void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);
void free_lvbmat_structure(DataStructure *p_lvbmat);
long brcnt(long n) { return (n << 1) - 3; }; /* return number of branches in unrooted binary tree structure containing n tips */

#endif // #ifdef NP_Implementation

#ifdef MPI_Implementation

#ifdef MAP_REDUCE_SINGLE
    #include "../../Lvb.h"
	#ifdef __cplusplus
		extern "C" int read_file(char *file_name, int n_file_type, Dataptr p_lvbmat, DataSeqPtr lvbmat_seq);
		extern "C" void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);
		extern "C" int read_parameters(Params *prms, int argc, char **argv);
	#endif

	int read_file(char *file_name, int n_file_type, Dataptr p_lvbmat, DataSeqPtr lvbmat_seq);
	void free_lvbmat_structure(DataStructure *p_lvbmat);
	int read_parameters(Params *prms, int argc, char **argv);
#else
	extern "C" int read_file(char *file_name, int n_file_type, Dataptr p_lvbmat, DataSeqPtr lvbmat_seq);
	extern "C" void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);
	extern "C" int read_parameters(Params *prms, int argc, char **argv);
	int read_file(char *file_name, int n_file_type, DataStructure *p_lvbmat, DataSeqStructure *p_lvbmat_seq);
	void free_lvbmat_structure(DataSeqStructure *p_lvbmat_seq, int n_size);
	int read_parameters(Params *prms, int argc, char **argv);
#endif

void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);
long brcnt(long n); /* return number of branches in unrooted binary tree structure containing n tips */

#endif // #ifdef MPI_Implementation

#endif /* READFILE_H_ */