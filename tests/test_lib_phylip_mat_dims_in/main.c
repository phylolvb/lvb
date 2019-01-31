/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Maximilian Strobl, Chris Wood, and Fernando Guntoro.
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

#include "lvb.h"

/* basic test of phylip_mat_dims_in() */

/* these constants are given in the infile */
#define EXPECTED_N 10
#define EXPECTED_M 63

int main(void)
{
    long m;				/* sites */
    long n;				/* sequences */
    int max_length_name;		/* mas name length */
    Lvb_bool success = LVB_FALSE;	/* test passed */
    Params rcstruct;		/* configurable parameters */

    lvb_initialize();

//    rcstruct.file_name_in = "infile";
    strcpy(rcstruct.file_name_in, "infile");
    rcstruct.n_file_format = FORMAT_PHYLIP;

    phylip_mat_dims_in(rcstruct.file_name_in, rcstruct.n_file_format, &n, &m, &max_length_name);
    if ((n == EXPECTED_N) && (m == EXPECTED_M))
    {
    	/* try it again and check it still works */
	phylip_mat_dims_in(rcstruct.file_name_in, rcstruct.n_file_format, &n, &m, &max_length_name);
	if ((n == EXPECTED_N) && (m == EXPECTED_M))
	    success = LVB_TRUE;
    }

    if (success == LVB_TRUE)
        printf("test passed\n");
    else
    	printf("test failed\n");

    return 0;
}
