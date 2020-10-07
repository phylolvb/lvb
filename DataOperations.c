#ifdef LVB_NP

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

/* ========== DataOperations.c - data MSA operations ========== */

#include "LVB.h"

static long constchar(Dataptr restrict MSA, Lvb_bool *const togo, const Lvb_bool verbose);
static void cutcols(Dataptr MSA, const Lvb_bool *const tocut, long n_columns_to_change);
static void logcut(const Lvb_bool *const cut, const long m);

static char *getstatev(const Dataptr MSA, const long k)
/* return pointer to string containing 1 instance of each character state in
 * column k of MSA, or NULL if more than MAXSTATES states are
 * found; does not include '-', '?', 'N' or 'X'; ignores the special meaning
 * of other (partial) ambiguity codes, so can return NULL if these are present;
 * N.B. string is static and will be overwritten by later calls */
{
    static char statev[MAXSTATES + 1];	/* array of states */
    long statec;			/* number of states */
    long i;				/* loop counter */

    /* clear record of states */
    statev[0] = '\0';
    statec = 0;

    /* update record of states for column k */
    for (i = 0; i < MSA->n; ++i){
	if (strchr(statev, (int) MSA->row[i][k]) == NULL){	/* new state */
	    if ((MSA->row[i][k] != '-') && (MSA->row[i][k] != '?')
		&& (MSA->row[i][k] != 'N') && (MSA->row[i][k] != 'X')) {
		statev[statec++] = MSA->row[i][k];
	    }
	    if (statec > MAXSTATES) return NULL;
	    statev[statec] = '\0';	/* for strchr() */
	}
    }
	
    return statev;

} /* end getstatev() */

long MinimumTreeLength(const Dataptr MSA)
/* return minimum length of any tree based on MSA; FIXME not quite right
 * with ambiguity codes */
{
    long minlen = 0;	/* return value */
    char *statev;	/* list of states in current character */
    long k;		/* loop counter */

    for (k = 0; k < MSA->m; ++k) {
    	statev = getstatev(MSA, k);
    	if (statev == NULL) minlen += MAXSTATES;
    	else minlen += strlen(statev) - 1;
    }
    return minlen;

} /* end MinimumTreeLength() */

/**********

=head1 DNAToBinary - CONVERT DNA TEXT MATRIX TO BINARY STATESET MATRIX

=head2 SYNOPSIS

    void DNAToBinary(const Dataptr mat, Lvb_bool fifthstate,
     unsigned char **enc_mat);

=head2 DESCRIPTION

Converts a MSA of sequence strings to a MSA of binary-encoded
statesets, where each of A, C, T, G and O (deletion) is represented by
a different bit. Ambiguous bases are converted to the union of all the
bases they may represent. C<?> is treated as totally ambiguous and
C<-> is either treated as <?> or as <O>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item mat

C<mat>C<->E<gt>C<m> and C<mat>C<->E<gt>C<n> give the number of bases in
each sequence and the number of sequences, respectively. Member
C<mat>C<->E<gt>C<row> points to the first element in an array of
pointers, each of which points to a sequence stored as a text string.

=item fifthstate

If C<LVB_TRUE>, treat gaps indicated by C<-> as identical to C<O>. Otherwise,
treat gaps indicated by C<-> as identical to <?>, i.e., totally ambiguous.

=back

=head3 OUTPUT

=over 4

=item enc_mat

C<enc_mat> must point to the first element in an array of
C<mat>C<->E<gt>C<n> pointers, each of which points to an allocated
array of C<mat>C<->E<gt>C<n> elements. On return,
C<enc_mat>[I<i>][I<j>] will give the binary-encoded stateset for
C<mat>C<->E<gt>C<row>[I<i>][I<j>], where I<i> is in the interval
[0..C<mat>C<->E<gt>C<n>-1] and I<j> is in the interval
[0..C<mat>C<->E<gt>C<m>-1].

=back

=cut

**********/

void DNAToBinary(Dataptr restrict mat, Lvb_bit_length **enc_mat)
/* convert MSA from string form to binary-encoded form, in which each
 * biological character occupies half a byte; the MSA is padded with
 * ambiguous data as required to ensure all bytes are initialised, but padding
 * will not contribute to tree length - allowing the optimization of White and
 * Holland (2011) */
{
    long i;		/* loop counter */
    long j;		/* loop counter */
    long k;		/* loop counter */
    long mat_offset;	/* current position within MSA row */
    char base;		/* current base as text character */
    Lvb_bit_length sitestate = 0U;	/* binary-encoded single state set */
    Lvb_bit_length enc_sitestates;	/* set of four binary-encoded state sets */

    for (i = 0; i < mat->n; i++)
    {
        for (j = 0; j < mat->nwords; j++)
		{
			enc_sitestates = 0U;
			for (k = 0; k < LENGTH_WORD; k++)
			{
				mat_offset = (j << LENGTH_WORD_BITS_MULTIPLY) + k;
				if (mat_offset >= mat->m)	/* padding required */
					base = 'N';
				else
					base = mat->row[i][(j << LENGTH_WORD_BITS_MULTIPLY) + k];	/* observed base required */

				/* unambiguous bases */
				if (base == 'A')
					sitestate = A_BIT;
				else if (base == 'C')
					sitestate = C_BIT;
				else if (base == 'G')
					sitestate = G_BIT;
				else if (base == 'T')
					sitestate = T_BIT;
				else if (base == 'U')	/* treat the same as 'U' */
					sitestate = T_BIT;

				/* ambiguous bases */
				else if (base == 'Y')
					sitestate = C_BIT | T_BIT;
				else if (base == 'R')
					sitestate = A_BIT | G_BIT;
				else if (base == 'W')
					sitestate = A_BIT | T_BIT;
				else if (base == 'S')
					sitestate = C_BIT | G_BIT;
				else if (base == 'K')
					sitestate = T_BIT | G_BIT;
				else if (base == 'M')
					sitestate = C_BIT | A_BIT;
				else if (base == 'B')
					sitestate = C_BIT | G_BIT | T_BIT;
				else if (base == 'D')
					sitestate = A_BIT | G_BIT | T_BIT;
				else if (base == 'H')
					sitestate = A_BIT | C_BIT | T_BIT;
				else if (base == 'V')
					sitestate = A_BIT | C_BIT | G_BIT;
				else if (base == 'N')
					sitestate = A_BIT | C_BIT | G_BIT | T_BIT;
				else if (base == 'X')
					sitestate = A_BIT | C_BIT | G_BIT | T_BIT;

				/* total ambiguity or deletion - now always the same as 'N' */
				else if ((base == '?') || (base == '-'))
					sitestate = A_BIT | C_BIT | G_BIT | T_BIT;

				/* all other bases are not allowed */
				else
				{
					fprintf(stderr, "bad base symbol in data MSA: ");
					fputc(base, stderr);
					fputc('\n', stderr);
					crash("cannot read data MSA");
				}

				lvb_assert(sitestate != 0U);
				enc_sitestates |= sitestate << (k << NIBBLE_WIDTH_BITS);
			}
			enc_mat[i][j] = enc_sitestates;
		}
    }
} /* end DNAToBinary() */


void rowfree(Dataptr MSA)
/* free memory used for row strings and array of row strings in MSA,
 * and make the array of row title strings NULL;
 * or, if the array of row title strings is already NULL, do nothing */
{
    long i;	/* loop counter */

    if (MSA->row != NULL) {
    	for(i = 0; i < MSA->n; ++i){
    		free(MSA->row[i]);
    		free(MSA->rowtitle[i]);
    	}
    	free(MSA->row);
    	free(MSA->rowtitle);
    	MSA->row = NULL;
    }

} /* end rowfree() */


static long constchar(Dataptr restrict MSA, Lvb_bool *const togo, const Lvb_bool verbose)
/* Make sure MSA->m-element array togo is LVB_TRUE where MSA column
 * contains only one character state;
 * log details of new columns to ignore if verbose is LVB_TRUE.
 * scratch must point to the first element of an array of at least MSA->m
 * elements or arbitrary (even uninitialised) contents. It will be left with
 * arbitrary contents on return. */
{
    long k;		/* loop counter */
    long i;		/* loop counter */
    long n_columns = 0;

    /* discover variable columns */
    for (k = 0; k < MSA->m; ++k){
    	for (i = 1; i < MSA->n; ++i){
			if (MSA->row[i][k] != MSA->row[0][k]){
				togo[k] = LVB_TRUE;
				n_columns += 1;
				break;
			}
		}
    }

    if (verbose == LVB_TRUE){
    	printf("Constant columns excluded from analysis: ");
    	if (n_columns == 0) printf(" none found.\n");
    	else logcut(togo, MSA->m);
    }
    return n_columns;
} /* end constchar() */

void matchange(Dataptr MSA, const Parameters rcstruct)
/* change and remove columns in MSA, partly in response to rcstruct,
 * verbosely or not according to value of verbose */
{
    Lvb_bool *togo;	/* LVB_TRUE where column must go */
    long n_columns_to_change = 0;
    /* Allocate memory: this will be free'd just before we return.
     * Dynamic allocation is used because otherwise, each of the
     * arrays would have to have MAX_M elements. That would either
     * limit the program too much, or causes massive waste of
     * address space. */

    togo = (Lvb_bool *) alloc(MSA->m * sizeof(Lvb_bool), "'togo' array");

    /* initialize all elements to LVB_FALSE ('don't ignore') */
    for(n_columns_to_change = 0; n_columns_to_change < MSA->m; n_columns_to_change ++) *(togo + n_columns_to_change) = LVB_FALSE;

    n_columns_to_change = constchar(MSA, togo, (Lvb_bool) rcstruct.verbose);	/* compuslory cut */

    /* N.B. a function to mark autapomorphic characters for cutting
     * could be called at this point. The effect would be more noticeable
     * with unrealistically small test matrices than with real data */

    /* cut the cols as indicated, and crash verbosely if too few remain */
    if (n_columns_to_change != MSA->m){
    	cutcols(MSA, togo, n_columns_to_change);	/* make changes to MSA */
    }
    else{
        MSA->bytes = bytes_per_row(MSA->m);
        MSA->nwords = words_per_row(MSA->m);
        MSA->tree_bytes = tree_bytes(MSA);
        MSA->tree_bytes_without_sitestate = tree_bytes_without_sitestate(MSA);
		MSA->min_len_tree = MinimumTreeLength(MSA);
    }
    if (MSA->m < MIN_M)
    	crash("after constant columns are ignored, data MSA has\n"
    			"%ld columns, which is less than LVB's lower limit of\n"
    			"%ld columns.\n", MSA->m, MIN_M);
    else{
    	if (rcstruct.verbose == LVB_TRUE) printf("\nIn total, %ld columns are excluded from the analysis\n\n", MSA->original_m - MSA->m);
    }

    /* free "local" dynamic heap memory */
    free(togo);

} /* end matchange() */

static void cutcols(Dataptr MSA, const Lvb_bool *const tocut, long n_columns_to_change)
/* remove columns in MSA for which the corresponding element of
MSA->m-element array tocut is LVB_TRUE, and update MSA->m;
return the number of columns cut */
{
    char **newrow;			/* rows of reduced MSA */
    long i;				/* loop counter */
    long k;				/* loop counter */
    long newk;				/* current column of reduced MSA */

    long uun = MSA->n;
    long uum = MSA->m;
    /* memory for new MSA row array */
    newrow = (char **) alloc((size_t) uun * sizeof(char *), "pointers to new row strings");
    for (i = 0; i < uun; ++i) newrow[i] =  (char*) alloc(sizeof(char) * (n_columns_to_change + 1), "new row strings");

    newk = 0;
    for (k = 0; k < uum; ++k){	/* for every column */
		if (tocut[k] == LVB_TRUE){	/* keep this column */
			for (i = 0; i < uun; ++i)	/* fill for ea. row */
				newrow[i][newk] = MSA->row[i][k];
			++newk;	/* fill next row next time */
		}
    }

    /* trap impossible condition */
    lvb_assert(newk == n_columns_to_change);

    /* terminate new row strings */
    for (i = 0; i < uun; ++i) newrow[i][newk] = '\0';

    /* update MSA structure */
    MSA->row = newrow;
    MSA->m = n_columns_to_change;
    MSA->bytes = bytes_per_row(MSA->m);
    MSA->nwords = words_per_row(MSA->m);
    MSA->tree_bytes = tree_bytes(MSA);
    MSA->tree_bytes_without_sitestate = tree_bytes_without_sitestate(MSA);
	MSA->min_len_tree = MinimumTreeLength(MSA);
} /* end cutcols() */

static void logcut(const Lvb_bool *const cut, const long m)
/* log message saying columns for which m-element array cut is LVB_TRUE are
 * being cut */
{
    long k;				/* loop counter */
    long noperln = 0;			/* no. of numbers on current line */
    const long max_noperln = 8; 	/* max. numbers written per line */
    Lvb_bool newline = LVB_FALSE;	/* last number followed by '\n' */

    printf("\n");

    /* give formatted list of columns to go */
    for (k = 0; k < m; ++k){
		if (cut[k] == LVB_FALSE){
			printf("%ld", k + 1L);
			++noperln;
			if (noperln == max_noperln){	/* end line */
				noperln = 0;
				printf("\n");
				newline = LVB_TRUE;
			}
			else{	/* just put some space on this line */
				printf("\t");
				newline = LVB_FALSE;
			}
		}
    }
    if (newline == LVB_FALSE)
	printf("\n");

    if (fflush(stdout) != 0)
	crash("write error on standard output"); /* FIXME: helpful? */

} /* end logcut() */

long words_per_row(const long m)
/* return the number of 32-bit words required to hold an encoded row of the
 * data MSA for m state sets, assuming half a byte per state set rounded up
 * to the nearest 32-bit word, which allows for the optimization of White and
 * Holland (2011, Bioinformatics 27:1359-1367, specifically Section 2.10) */
{
    long words;		/* 32-bit words required */

    words = m >> LENGTH_WORD_BITS_MULTIPLY;
    if (m % LENGTH_WORD) words += 1;
    return words;
}

long bytes_per_row(const long m)
/* return the number of 8-bit bytes required to hold an encoded row of the
 * data MSA for m state sets, assuming half a byte per state set rounded up
 * to the nearest 32-bit word - which allows for the optimization of White and
 * Holland (2011, Bioinformatics 27:1359-1367, specifically Section 2.10) */
{
	return words_per_row(m) * sizeof(Lvb_bit_length);
}

#elif LVB_PARALLEL_SEARCH

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

/* ========== DataOperations.c - data MSA operations ========== */

#include "LVB.h"

#ifdef MPI_Implementation

static long constchar(Dataptr restrict MSA, Lvb_bool *const togo, const Lvb_bool verbose);
	static void cutcols(Dataptr MSA, const Lvb_bool *const tocut, long n_columns_to_change);
	static void logcut(const Lvb_bool *const cut, const long m);

static char *getstatev(const Dataptr MSA, const long k)
	/* return pointer to string containing 1 instance of each character state in
	 * column k of MSA, or NULL if more than MAXSTATES states are
	 * found; does not include '-', '?', 'N' or 'X'; ignores the special meaning
	 * of other (partial) ambiguity codes, so can return NULL if these are present;
	 * N.B. string is static and will be overwritten by later calls */
	{
	    static char statev[MAXSTATES + 1];	/* array of states */
	    long statec;			/* number of states */
	    long i;				/* loop counter */

	    /* clear record of states */
	    statev[0] = '\0';
	    statec = 0;

	    /* update record of states for column k */
	    for (i = 0; i < MSA->n; ++i){
			if (strchr(statev, (int) MSA->row[i][k]) == NULL){	/* new state */
				if ((MSA->row[i][k] != '-') && (MSA->row[i][k] != '?') && (MSA->row[i][k] != 'N') && (MSA->row[i][k] != 'X')) {
					statev[statec++] = MSA->row[i][k];
	    		}
	    		if (statec > MAXSTATES) return NULL;
	    		statev[statec] = '\0';	/* for strchr() */
	    	}
	    }
	    return statev;
	} /* end getstatev() */

	long MinimumTreeLength(const Dataptr MSA)
	/* return minimum length of any tree based on MSA; FIXME not quite right
	 * with ambiguity codes */
	{
		long minlen = 0;	/* return value */
		char *statev;	/* list of states in current character */
		long k;		/* loop counter */

		for (k = 0; k < MSA->m; ++k) {
			statev = getstatev(MSA, k);
			if (statev == NULL) minlen += MAXSTATES;
			else minlen += strlen(statev) - 1;
		}
		return minlen;

	} /* end MinimumTreeLength() */




/**********

=head1 DNAToBinary - CONVERT DNA TEXT MATRIX TO BINARY STATESET MATRIX

=head2 SYNOPSIS

    void DNAToBinary(const Dataptr mat, Lvb_bool fifthstate,
     unsigned char **enc_mat);

=head2 DESCRIPTION

Converts a MSA of sequence strings to a MSA of binary-encoded
statesets, where each of A, C, T, G and O (deletion) is represented by
a different bit. Ambiguous bases are converted to the union of all the
bases they may represent. C<?> is treated as totally ambiguous and
C<-> is either treated as <?> or as <O>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item mat

C<mat>C<->E<gt>C<m> and C<mat>C<->E<gt>C<n> give the number of bases in
each sequence and the number of sequences, respectively. Member
C<mat>C<->E<gt>C<row> points to the first element in an array of
pointers, each of which points to a sequence stored as a text string.

=item fifthstate

If C<LVB_TRUE>, treat gaps indicated by C<-> as identical to C<O>. Otherwise,
treat gaps indicated by C<-> as identical to <?>, i.e., totally ambiguous.

=back

=head3 OUTPUT

=over 4

=item enc_mat

C<enc_mat> must point to the first element in an array of
C<mat>C<->E<gt>C<n> pointers, each of which points to an allocated
array of C<mat>C<->E<gt>C<n> elements. On return,
C<enc_mat>[I<i>][I<j>] will give the binary-encoded stateset for
C<mat>C<->E<gt>C<row>[I<i>][I<j>], where I<i> is in the interval
[0..C<mat>C<->E<gt>C<n>-1] and I<j> is in the interval
[0..C<mat>C<->E<gt>C<m>-1].

=back

=cut

**********/

void DNAToBinary(Dataptr restrict mat, Lvb_bit_length **enc_mat)

/* convert MSA from string form to binary-encoded form, in which each
 * biological character occupies half a byte; the MSA is padded with
 * ambiguous data as required to ensure all bytes are initialised, but padding
 * will not contribute to tree length - allowing the optimization of White and
 * Holland (2011) */
{
    long i;		/* loop counter */
    long j;		/* loop counter */
    long k;		/* loop counter */
    long mat_offset;	/* current position within MSA row */
    char base;		/* current base as text character */
    Lvb_bit_length sitestate = 0U;	/* binary-encoded single state set */
    Lvb_bit_length enc_sitestates;	/* set of four binary-encoded state sets */

    for (i = 0; i < mat->n; i++)
    {
        for (j = 0; j < mat->nwords; j++)
		{
			enc_sitestates = 0U;
			for (k = 0; k < LENGTH_WORD; k++)
			{
				mat_offset = (j << LENGTH_WORD_BITS_MULTIPLY) + k;
				if (mat_offset >= mat->m)	/* padding required */
					base = 'N';
				else
					base = mat->row[i][(j << LENGTH_WORD_BITS_MULTIPLY) + k];	/* observed base required */

				/* unambiguous bases */
				if (base == 'A')
					sitestate = A_BIT;
				else if (base == 'C')
					sitestate = C_BIT;
				else if (base == 'G')
					sitestate = G_BIT;
				else if (base == 'T')
					sitestate = T_BIT;
				else if (base == 'U')	/* treat the same as 'U' */
					sitestate = T_BIT;

				/* ambiguous bases */
				else if (base == 'Y')
					sitestate = C_BIT | T_BIT;
				else if (base == 'R')
					sitestate = A_BIT | G_BIT;
				else if (base == 'W')
					sitestate = A_BIT | T_BIT;
				else if (base == 'S')
					sitestate = C_BIT | G_BIT;
				else if (base == 'K')
					sitestate = T_BIT | G_BIT;
				else if (base == 'M')
					sitestate = C_BIT | A_BIT;
				else if (base == 'B')
					sitestate = C_BIT | G_BIT | T_BIT;
				else if (base == 'D')
					sitestate = A_BIT | G_BIT | T_BIT;
				else if (base == 'H')
					sitestate = A_BIT | C_BIT | T_BIT;
				else if (base == 'V')
					sitestate = A_BIT | C_BIT | G_BIT;
				else if (base == 'N')
					sitestate = A_BIT | C_BIT | G_BIT | T_BIT;
				else if (base == 'X')
					sitestate = A_BIT | C_BIT | G_BIT | T_BIT;

				/* total ambiguity or deletion - now always the same as 'N' */
				else if ((base == '?') || (base == '-'))
					sitestate = A_BIT | C_BIT | G_BIT | T_BIT;

				/* all other bases are not allowed */
				else
				{
					fprintf(stderr, "bad base symbol in data MSA: ");
					fputc(base, stderr);
					fputc('\n', stderr);
					crash("cannot read data MSA");
				}

				lvb_assert(sitestate != 0U);
				enc_sitestates |= sitestate << (k << NIBBLE_WIDTH_BITS);
			}
			enc_mat[i][j] = enc_sitestates;
		}
    }
} /* end DNAToBinary() */


	void rowfree(Dataptr MSA, int n_lines)
	/* free memory used for row strings and array of row strings in MSA,
	 * and make the array of row title strings NULL;
	 * or, if the array of row title strings is already NULL, do nothing */
	{
		long i;	/* loop counter */

		if (MSA->row != NULL) {
			for(i = 0; i < n_lines; ++i){
				if (MSA->row != 0L) free(MSA->row[i]);
				free(MSA->rowtitle[i]);
			}
			if (MSA->row != 0L) free(MSA->row);
			free(MSA->rowtitle);
			MSA->row = NULL;
		}

	} /* end rowfree() */


	static long constchar(Dataptr restrict MSA, Lvb_bool *const togo, const Lvb_bool verbose)
/* Make sure MSA->m-element array togo is LVB_TRUE where MSA column
 * contains only one character state;
 * log details of new columns to ignore if verbose is LVB_TRUE.
 * scratch must point to the first element of an array of at least MSA->m
 * elements or arbitrary (even uninitialised) contents. It will be left with
 * arbitrary contents on return. */
{
    long k;		/* loop counter */
    long i;		/* loop counter */
    long n_columns = 0;

    /* discover variable columns */
    for (k = 0; k < MSA->m; ++k){
    	for (i = 1; i < MSA->n; ++i){
    		if (MSA->row[i][k] != MSA->row[0][k]){
				togo[k] = LVB_TRUE;
				n_columns += 1;
				break;
			}
		}
    }

    if (verbose == LVB_TRUE){
    	printf("Ignoring constant columns\n");
    	if (n_columns == 0) printf("... none found.\n");
    	else logcut(togo, MSA->m);
    }
    return n_columns;
} /* end constchar() */

    void matchange(Dataptr MSA, const Parameters rcstruct)
/* change and remove columns in MSA, partly in response to rcstruct,
 * verbosely or not according to value of verbose */
{
    Lvb_bool *togo;	/* LVB_TRUE where column must go */
    long n_columns_to_change = 0;
    /* Allocate memory: this will be free'd just before we return.
     * Dynamic allocation is used because otherwise, each of the
     * arrays would have to have MAX_M elements. That would either
     * limit the program too much, or causes massive waste of
     * address space. */

    togo = (Lvb_bool *) alloc(MSA->m * sizeof(Lvb_bool), "'togo' array");

    /* initialize all elements to LVB_FALSE ('don't ignore') */
    for(n_columns_to_change = 0; n_columns_to_change < MSA->m; n_columns_to_change ++) *(togo + n_columns_to_change) = LVB_FALSE;

    n_columns_to_change = constchar(MSA, togo, (Lvb_bool) rcstruct.verbose);	/* compuslory cut */
    /* N.B. a function to mark autapomorphic characters for cutting
	 * could be called at this point. The effect would be more noticable
	 * with unrealistically small test matrices than with real data */

	/* cut the cols as indicated, and crash verbosely if too few remain */
	if (n_columns_to_change != MSA->m){
		cutcols(MSA, togo, n_columns_to_change);	/* make changes to MSA */
	}
	else{
		MSA->bytes = bytes_per_row(MSA->m);
		MSA->nwords = words_per_row(MSA->m);
		MSA->tree_bytes = tree_bytes(MSA);
		MSA->tree_bytes_without_sitestate = tree_bytes_without_sitestate(MSA);
		MSA->min_len_tree = MinimumTreeLength(MSA);
	}
	if (MSA->m < MIN_M)
		crash("after constant columns are ignored, data MSA has\n"
				"%ld columns, which is less than LVB's lower limit of\n"
				"%ld columns.\n", MSA->m, MIN_M);
	else{
		if (rcstruct.verbose == LVB_TRUE) printf("A total of %ld columns will be ignored\n", MSA->original_m - MSA->m);
	}

    /* free "local" dynamic heap memory */
    free(togo);

} /* end matchange() */

    static void cutcols(Dataptr MSA, const Lvb_bool *const tocut, long n_columns_to_change)
/* remove columns in MSA for which the corresponding element of
MSA->m-element array tocut is LVB_TRUE, and update MSA->m;
return the number of columns cut */
{
    char **newrow;			/* rows of reduced MSA */
    long i;				/* loop counter */
    long k;				/* loop counter */
    long newk;				/* current column of reduced MSA */

    long uun = MSA->n;
    long uum = MSA->m;
    /* memory for new MSA row array */
    newrow = (char **) alloc((size_t) uun * sizeof(char *), "pointers to new row strings");
    for (i = 0; i < uun; ++i) newrow[i] =  (char*) alloc(sizeof(char) * (n_columns_to_change + 1), "new row strings");

    newk = 0;
    for (k = 0; k < uum; ++k){	/* for every column */
		if (tocut[k] == LVB_TRUE){	/* keep this column */
			for (i = 0; i < uun; ++i)	/* fill for ea. row */
				newrow[i][newk] = MSA->row[i][k];
			++newk;	/* fill next row next time */
		}
    }

    /* trap impossible condition */
    lvb_assert(newk == n_columns_to_change);

    /* terminate new row strings */
    for (i = 0; i < uun; ++i) newrow[i][newk] = '\0';

    /* update MSA structure */
    MSA->m = n_columns_to_change;
    MSA->bytes = bytes_per_row(MSA->m);
    MSA->nwords = words_per_row(MSA->m);
    MSA->tree_bytes = tree_bytes(MSA);
    MSA->tree_bytes_without_sitestate = tree_bytes_without_sitestate(MSA);
    MSA->row = newrow;
    MSA->min_len_tree = MinimumTreeLength(MSA);
} /* end cutcols() */


void get_bootstrap_weights(long *weight_arr, long m, long extras)
/* Fill first m elements of array whose first element is pointed to by
 * weight_arr with weights for a single bootstrap resample. This is
 * obtained on the assumption that extras constant characters were in
 * the original sequence, but are not represented in weight_arr. This
 * gives a bootstrap sample with these constant characters effectively
 * included. */
{
    long samples = 0;	/* size of the sample so far */
    long site;		/* number of current site to add to sample */

    memset(weight_arr, 0, m * sizeof(long));

    while (samples < (m + extras)){
    	site = randpint(m + extras - 1);
    	if (site < m) weight_arr[site] += 1;
    	samples++;
    }

} /* end get_bootstrap_weights() */

static void logcut(const Lvb_bool *const cut, const long m)
/* log message saying columns for which m-element array cut is LVB_TRUE are
 * being cut */
{
    long k;				/* loop counter */
    long noperln = 0;			/* no. of numbers on current line */
    const long max_noperln = 8; 	/* max. numbers written per line */
    Lvb_bool newline = LVB_FALSE;	/* last number followed by '\n' */

    printf("... will ignore column numbers:\n");

    /* give formatted list of columns to go */
    for (k = 0; k < m; ++k){
		if (cut[k] == LVB_FALSE){
			printf("%ld", k + 1L);
			++noperln;
			if (noperln == max_noperln){	/* end line */
				noperln = 0;
				printf("\n");
				newline = LVB_TRUE;
			}
			else{	/* just put some space on this line */
				printf("\t");
				newline = LVB_FALSE;
			}
		}
    }
    if (newline == LVB_FALSE)
	printf("\n");

    if (fflush(stdout) != 0)
	crash("write error on standard output"); /* FIXME: helpful? */

} /* end logcut() */


long words_per_row(const long m)
/* return the number of 32-bit words required to hold an encoded row of the
 * data MSA for m state sets, assuming half a byte per state set rounded up
 * to the nearest 32-bit word, which allows for the optimization of White and
 * Holland (2011, Bioinformatics 27:1359-1367, specifically Section 2.10) */
{
    long words;		/* 32-bit words required */

    words = m >> LENGTH_WORD_BITS_MULTIPLY;
    if (m % LENGTH_WORD) words += 1;
    return words;
}



long bytes_per_row(const long m)
/* return the number of 8-bit bytes required to hold an encoded row of the
 * data MSA for m state sets, assuming half a byte per state set rounded up
 * to the nearest 32-bit word - which allows for the optimization of White and
 * Holland (2011, Bioinformatics 27:1359-1367, specifically Section 2.10) */
{
	return words_per_row(m) * sizeof(Lvb_bit_length);
}

#endif


#endif