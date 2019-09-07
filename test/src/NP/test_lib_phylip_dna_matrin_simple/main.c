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

#include "Lvb.h"

/* Positive test that a simple (one line per sequence) matrix may be
 * read, both as a sequential and as an interleaved matrix. Example is
 * taken from the PHYLIP 3.6a documentation. */

static const char *name_expected[6] =
{
    "Archaeopt",
    "Hesperorni",
    "Baluchithe",
    "B. virgini",
    "Brontosaur",
    "B.subtilis"
};

static const char *sequence_expected[6] =
{
    "CGATGCTTACCGC",
    "CGTTACTCGTTGT",
    "TAATGTTAATTGT",
    "TAATGTTCGTTGT",
    "CAAAACCCATCAT",
    "GGCAGCCAATCAC"
};

static void check(Dataptr matrix)
/* check matrix has been input correctly, or crash if not */
{
    long i;		/* loop counter */
    lvb_assert(matrix->m == 13);
    lvb_assert(matrix->n == 6);

    for (i = 0; i < 6; i++)
    {
        lvb_assert(strlen(matrix->row[i]) == 13);
 	lvb_assert(strlen(matrix->rowtitle[i]) == strlen(name_expected[i]));
	lvb_assert(strcmp(matrix->row[i], sequence_expected[i]) == 0);
	lvb_assert(strcmp(matrix->rowtitle[i], name_expected[i]) == 0);
    }
}

int main(void)
{
    Dataptr matrix1;	/* data matrix as input (interleaved) */
    Dataptr matrix2;	/* data matrix as input (sequential) */

    Params rcstruct;		/* configurable parameters */
    strcpy(rcstruct.file_name_in, "infile");
    rcstruct.n_file_format = FORMAT_PHYLIP;

    lvb_initialize();

    matrix1 = malloc(sizeof(DataStructure));
    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, matrix1);
    check(matrix1);

    matrix2 = malloc(sizeof(DataStructure));
    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, matrix2);
    check(matrix2);

    rowfree(matrix1);
    rowfree(matrix2);
    printf("test passed\n");
    return 0;
}
