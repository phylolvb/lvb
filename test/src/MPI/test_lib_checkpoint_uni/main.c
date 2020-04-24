/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maxi
milian Strobl and Chris Wood.
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


#include "lvb.h"


/* test of checkpoint_uni() and restore_uni() */

#define RAND_NOS 27201		/* some large number */
#define SEED 40291989		/* arbitrary */

#ifdef CHECKPOINT_INTERVAL
	#undef CHECKPOINT_INTERVAL
#endif
#define CHECKPOINT_INTERVAL 3	/* arbitrary - low for a thorough test */

static double no_checkpoint[RAND_NOS];
static double with_checkpoint[RAND_NOS];

int main(int argc, char **argv)
{
    FILE *fp;			/* checkpoint file */
    long i;			/* loop counter */
    int my_id;			/* MPI process ID */

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (my_id == 0) {
    	rinit(SEED);
    	for (i = 0; i < RAND_NOS; i++){
    		no_checkpoint[i] = uni();
    	}
    	rinit(SEED);
    	for (i = 0; i < RAND_NOS; i++)
    	{
    		if ((i % CHECKPOINT_INTERVAL) == 0) {
    			fp = fopen("uni_checkpoint", "wb");
    			checkpoint_uni(fp);

    			lvb_assert(fclose(fp) == 0);
    			fp = fopen("uni_checkpoint", "rb");
    			restore_uni(fp);

    			/* point to begin again */
    			fseek(fp, 0L, SEEK_SET);
    			lvb_assert(test_block_data(fp) == LVB_TRUE);		/* test if block is correct */
    			lvb_assert(fclose(fp) == 0);
    		}
    		with_checkpoint[i] = uni();
    	}
    	remove("uni_checkpoint");

    	if (memcmp(no_checkpoint, with_checkpoint, RAND_NOS * sizeof(double)) == 0) {
    		printf("test passed\n");
    	}
    	else {
    		printf("FATAL ERROR\n");
    	}
    }
    MPI_Finalize();
}
