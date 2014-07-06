/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

wrapper.c - LVB to PHYLIP interface

=cut

**********/

#include "lvb.h"
#include "../PHYLIP_FOR_LVB/src/phylip.h"
#include "../PHYLIP_FOR_LVB/src/seq.h"

/**********

=head1 phylip_dna_matrin - READ PHYLIP-FORMAT DNA DATA MATRIX

=head2 SYNOPSIS

    Dataptr phylip_dna_matrin(Lvb_bool ileaved);

=head2 DESCRIPTION

Read a DNA data matrix in PHYLIP 3.6 format from file. The file name is
given by the macro MATFNAM in F<lvb.h>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item ileaved

C<LVB_TRUE> indicates that the matrix is in PHYLIP interleaved format.
Otherwise, it is assumed to be in PHYLIP sequential format.

=back

=head2 RETURN

Returns a pointer to a new, dynamically allocated LVB data matrix
structure containing the data matrix.

=cut

**********/

Dataptr phylip_dna_matrin(Lvb_bool ileaved)
{
    long i;		/* loop counter */
    long j;		/* loop counter */
    Dataptr lvbmat;	/* return value */

    /* indicate to PHYLIP whether the matrix is in interleaved or
     * sequential format */
    interleaved = ileaved;

    /* read the matrix and associated information in to PHYLIP's
     * globals spp, chars, y and nayme */
    dnapars_wrapper();

    /* check number of sequences is in range for LVB */
    if (spp < MIN_N)
	crash("The data matrix must have at least %ld sequences.",
	 MIN_N);
    if (spp > MAX_N)
	crash("The data matrix must have no more than %ld sequences.",
	 MAX_N);

    /* check number of sites is in range for LVB */
    else if (chars < MIN_M)
	crash("The data matrix must have at least %ld sites.", MIN_N);
    else if (chars > MAX_M)
	crash("The data matrix must have no more than %ld sites.",
	 MAX_M);

    /* transfer copies of PHYLIP data to LVB data structure */

    lvbmat = matalloc(spp);

    /* we want null-terminated strings, so we cannot simply point to
     * the same, non-null-terminated arrays as are found in PHYLIP's
     * data structures */
    for (i = 0; i < spp; i++)
    {
        lvbmat->rowtitle[i] = salloc(nmlngth, "sequence names");
        lvbmat->row[i] = salloc(chars, "sequences");
    }
    for (i = 0; i < spp; i++)
    {
        for (j = 0; j < nmlngth; j++)
	    lvbmat->rowtitle[i][j] = nayme[i][j];
	lvbmat->rowtitle[i][nmlngth] = '\0';
    }
    for (i = 0; i < spp; i++)
    {
        for (j = 0; j < chars; j++)
	    lvbmat->row[i][j] = y[i][j];
	lvbmat->row[i][chars] = '\0';
    }

    /* scalars */
    lvbmat->m = chars;
    lvbmat->n = spp;

    return lvbmat;
} /* end phylip_dna_matrin() */

/**********

=head1 phylip_mat_dims_in - READ PHYLIP MATRIX DIMENSIONS

=head2 SYNOPSIS

void phylip_mat_dims_in(long *species_ptr, long *sites_ptr);

=head2 DESCRIPTION

Reads the number of sequences and number of sites per sequence in a
PHYLIP-format matrix file.

=head2 PARAMETERS

=head3 OUTPUT

=over 4

=item species_ptr

On return, *C<species_ptr> is set to the number of sequences.

=item sites_ptr

On return, *C<sites_ptr> is set to the number of sites per
sequence.

=back

=head2 BUGS

C<phylip_mat_dims_in()> should not be called if the matrix file is
open.

=cut

**********/

void phylip_mat_dims_in(long *species_ptr, long *sites_ptr)
{
    FILE *old_infile;		/* temp. holder for global's value */
    long nonodes_dummy; 	/* required by inputnumbers */
    const long n_dummy = 1;	/* required by inputnumbers */

    old_infile = infile;	/* preserve original value */
    infile = clnopen(INFILE, "r");	/* assign to PHYLIP global */
    inputnumbers(species_ptr, sites_ptr, &nonodes_dummy, n_dummy);
    clnclose(infile, INFILE);
    infile = old_infile;	/* restore value */

} /* end phylip_mat_dimsin() */
