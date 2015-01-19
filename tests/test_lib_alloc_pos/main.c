/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
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

/* An array of TEST_ALLOC_LONGS long ints is allocated and the
 * test expects this allocation to succeed. */
#define TEST_ALLOC_LONGS 1000000

/* Subsequently, that array is freed and an array of TEST_ALLOC_CHARS
 * characters is allocated, and the test expects this to succeed as
 * well. */
#define TEST_ALLOC_CHARS 70000

long *lp = NULL;	/* pointer to start of an allocated array */
char *cp = NULL;	/* pointer to start of an allocated array */

int main(void)
{
    long i;	/* loop counter */

    lvb_initialize();

    /* Check that a zero byte allocation returns a NULL pointer. */
    lp = alloc(0, "zero-byte array");
    lvb_assert(lp == NULL);

    /* Check that a fairly large allocation succeeds and that the
     * memory allocated may be written to. */
    lp = alloc(TEST_ALLOC_LONGS * sizeof(long), "test array 1");
    lvb_assert(lp != NULL);
    for (i = 0; i < TEST_ALLOC_LONGS; i++)
	lp[i] = 1;

    free(lp);

    /* Basic check that heap is OK after free: allocate something else
     * and check we may write to it */
    cp = alloc(TEST_ALLOC_CHARS, "test array 2");
    for (i = 0; i < TEST_ALLOC_CHARS; i++)
    	cp[i] = 'X';

    free(cp);
    
    printf("test passed\n");
    return EXIT_SUCCESS;
}
