/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include "lvb.h"

/* Positive test that a simple (one line per sequence) matrix may be
 * read, both as a sequential and as an interleaved matrix. Example is
 * taken from the PHYLIP 3.6a documentation. */

static const char *name_expected[6] =
{
    "Archaeopt ",
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
        lvb_assert(strlen(matrix->rowtitle[i]) == 10);
	lvb_assert(strcmp(matrix->row[i], sequence_expected[i]) == 0);
	lvb_assert(strcmp(matrix->rowtitle[i], name_expected[i]) == 0);
    }
}

int main(void)
{
    Dataptr matrix1;	/* data matrix as input (interleaved) */
    Dataptr matrix2;	/* data matrix as input (sequential) */

    lvb_initialize();
    matrix1 = phylip_dna_matrin(LVB_TRUE);
    check(matrix1);
    matrix2 = phylip_dna_matrin(LVB_FALSE);
    check(matrix2);

    printf("test passed\n");
    return 0;
}
