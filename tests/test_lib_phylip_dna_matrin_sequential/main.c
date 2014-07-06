/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include "lvb.h"

/* Positive test that a sequential matrix may be read. Example is taken
 * from the PHYLIP 3.6a documentation. */

static const char *name_expected[5] =
{
    "Turkey    ",
    "Salmo gair",
    "H. Sapiens",
    "Chimp     ",
    "Gorilla   "
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
    Dataptr matrix;	/* data matrix as input */
    long i;		/* loop counter */

    lvb_initialize();
    matrix = phylip_dna_matrin(LVB_FALSE);
    lvb_assert(matrix->m == 42);
    lvb_assert(matrix->n == 5);

    for (i = 0; i < 5; i++)
    {
        lvb_assert(strlen(matrix->row[i]) == 42);
        lvb_assert(strlen(matrix->rowtitle[i]) == 10);
	lvb_assert(strcmp(matrix->row[i], sequence_expected[i]) == 0);
	lvb_assert(strcmp(matrix->rowtitle[i], name_expected[i]) == 0);
    }

    printf("test passed\n");
    return 0;
}
