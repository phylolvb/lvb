/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and 
Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2020 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
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

#include "src/LVB.h"

/* Positive test that an interleaved MSA may be read. Example is taken
 * from the PHYLIP 3.6a documentation. */

static const char *name_expected[5] =
{
    "Turkey",
    "Salmo gair",
    "H. Sapiens",
    "Chimp",
    "Gorilla"
};

static const char *sequence_expected[5] =
{
    "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT",
    "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT",
    "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA",
    "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT",
    "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA"
};

int main(void)
{
    Dataptr MSA;	/* data MSA as input */
    long i;		/* loop counter */

    lvb_initialize();
    Parameters rcstruct;		/* configurable parameters */
    strcpy(rcstruct.file_name_in, "infile");
    rcstruct.n_file_format = FORMAT_PHYLIP;

    lvb_initialize();
    MSA = (data *) malloc(sizeof(DataStructure));
    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, MSA);
    lvb_assert(MSA->m == 42);
    lvb_assert(MSA->n == 5);

    for (i = 0; i < 5; i++)
    {
        lvb_assert(strlen(MSA->row[i]) == 42);
        lvb_assert(strlen(MSA->rowtitle[i]) == strlen(name_expected[i]));
	lvb_assert(strcmp(MSA->row[i], sequence_expected[i]) == 0);
	lvb_assert(strcmp(MSA->rowtitle[i], name_expected[i]) == 0);
    }

    rowfree(MSA);
    printf("test passed\n");
    return 0;
}
