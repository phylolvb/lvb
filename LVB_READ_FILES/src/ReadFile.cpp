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

#include "ReadFile.h"

void read_file(char *file_name, int n_file_type, Dataptr p_lvbmat){

	CReadFiles readFiles = CReadFiles();
	/// read file
	std::string sz_file_name = std::string(file_name);
	readFiles.read_file(file_name, n_file_type);

    /* scalars */
	//cout << (long) readFiles.get_length_sequences() << endl;
    p_lvbmat->m = (long) readFiles.get_length_sequences();
    p_lvbmat->original_m = p_lvbmat->m;
    p_lvbmat->n = (long) readFiles.get_number_seqs();
    p_lvbmat->nbranches = brcnt(p_lvbmat->n); 		/* possible number of braches */

    /* array for row title strings */
    p_lvbmat->rowtitle = (char **) malloc((p_lvbmat->n) * sizeof(char *));
    if (p_lvbmat->rowtitle == NULL) readFiles.exit_error(1 , "Fail to allocate memory...");

    /* array for row strings */
    p_lvbmat->row = (char **) malloc((p_lvbmat->n) * sizeof(char *));
    if (p_lvbmat->row == NULL) readFiles.exit_error(1 , "Fail to allocate memory...");

    /* we want null-terminated strings, so we cannot simply point to
     * the same, non-null-terminated arrays as are found in PHYLIP's
     * data structures */
    for (int i = 0; i < p_lvbmat->n; i++){
    	p_lvbmat->rowtitle[i] = (char*) malloc(sizeof(char) * (readFiles.get_max_length_seq_name() + 1));
    	p_lvbmat->row[i] = (char*) malloc(sizeof(char) * (p_lvbmat->m + 1));
    }
    for (int i = 0; i < p_lvbmat->n; i++) {
        for (int j = 0; j < p_lvbmat->m; j++) p_lvbmat->row[i][j] = readFiles.get_char_sequences(i, j);
        p_lvbmat->row[i][p_lvbmat->m] = '\0';
    }
    for (int i = 0; i < p_lvbmat->n; i++) {
        for (int j = 0; j < readFiles.get_length_seq_name(i); j++) p_lvbmat->rowtitle[i][j] = readFiles.get_char_seq_name(i, j);
        p_lvbmat->rowtitle[i][readFiles.get_length_seq_name(i)] = '\0';
    }
/*	std::string file_name_out = "/home/mmp/Downloads/file_nexus_nex_out.fas";
	readFiles.save_file(file_name_out);*/
}


void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name){

	CReadFiles readFiles = CReadFiles();
	/// read file
	std::string sz_file_name = std::string(file_name);
	readFiles.read_file(file_name, n_file_type);

	*species_ptr = (long) readFiles.get_number_seqs();
	*sites_ptr = (long) readFiles.get_length_sequences();
	*max_length_name = readFiles.get_max_length_seq_name();
//	std::string file_name_out = "/home/mmp/workspace/LVB_READ_FILES/aln_files/file_nexus_nex_out.fas";
//	readFiles.save_file(file_name_out);

}


void free_lvbmat_structure(DataStructure *p_lvbmat){

	if (p_lvbmat->row != NULL){
		for(int i = 0; i < p_lvbmat->n; ++i) free(p_lvbmat->row[i]);
		free(p_lvbmat->row);
		p_lvbmat->row = NULL;
	}
	if (p_lvbmat->rowtitle != NULL){
		for(int i = 0; i < p_lvbmat->n; ++i) free(p_lvbmat->rowtitle[i]);
		free(p_lvbmat->rowtitle);
		p_lvbmat->rowtitle = NULL;
	}
//	free(p_lvbmat);
//	p_lvbmat = NULL;
}


