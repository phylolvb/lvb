/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

datops.c - data matrix operations

=cut

**********/

#include "lvb.h"

static long constchar(const Dataptr matrix, Lvb_bool *const togo, const Lvb_bool verbose);
static void cutcols(Dataptr matrix, const Lvb_bool *const tocut, long n_columns_to_change);
static void logcut(const Lvb_bool *const cut, const long m);

static char *getstatev(const Dataptr matrix, const long k)
/* return pointer to string containing 1 instance of each character state in
 * column k of matrix, or NULL if more than MAXSTATES states are
 * found; ignores the special meaning of ambiguity codes or gaps, so can return
 * NULL if these are present;
 * N.B. string is static and will be overwritten by later calls */
{
    static char statev[MAXSTATES + 1];	/* array of states */
    long statec;			/* number of states */
    long i;				/* loop counter */

    /* clear record of states */
    statev[0] = '\0';
    statec = 0;

    /* update record of states for column k */
    for (i = 0; i < matrix->n; ++i){
		if (strchr(statev, (int) matrix->row[i][k]) == NULL){	/* new state */
			statev[statec++] = matrix->row[i][k];
			if (statec > MAXSTATES) return NULL;
			statev[statec] = '\0';	/* for strchr() */
		}
    }
	
    return statev;

} /* end getstatev() */

long getminlen(const Dataptr matrix)
/* return minimum length of any tree based on matrix; FIXME not quite right
 * with ambiguity or gaps */
{
    long minlen = 0;	/* return value */
    char *statev;	/* list of states in current character */
    long k;		/* loop counter */

    for (k = 0; k < matrix->m; ++k) {
    	statev = getstatev(matrix, k);
    	if (statev == NULL) minlen += MAXSTATES;
    	else minlen += strlen(statev) - 1;
    }
    return minlen;

} /* end getminlen() */

/**********

=head1 dna_makebin - CONVERT DNA TEXT MATRIX TO BINARY STATESET MATRIX

=head2 SYNOPSIS

    void dna_makebin(const Dataptr mat, Lvb_bool fifthstate,
     unsigned char **enc_mat);

=head2 DESCRIPTION

Converts a matrix of sequence strings to a matrix of binary-encoded
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

void dna_makebin(const Dataptr mat, Lvb_bool fifthstate, unsigned char **enc_mat)
{
    long i;			/* loop counter */
    long j;			/* loop counter */
    char base;			/* current base as text character */
    unsigned char sset = 0U;	/* binary-encoded state set */

    for (i = 0; i < mat->n; i++){
        for (j = 0; j < mat->m; j++) {
			base = mat->row[i][j];

			/* unambiguous bases */
			if (base == 'A') sset = A_BIT;
			else if (base == 'C') sset = C_BIT;
			else if (base == 'G') sset = G_BIT;
			else if (base == 'T') sset = T_BIT;
			else if (base == 'U')	/* treat the same as 'U' */
				sset = T_BIT;

			/* ambiguous bases */
			else if (base == 'Y') sset = C_BIT | T_BIT;
			else if (base == 'R') sset = A_BIT | G_BIT;
			else if (base == 'W') sset = A_BIT | T_BIT;
			else if (base == 'S') sset = C_BIT | G_BIT;
			else if (base == 'K') sset = T_BIT | G_BIT;
			else if (base == 'M') sset = C_BIT | A_BIT;
			else if (base == 'B') sset = C_BIT | G_BIT | T_BIT;
			else if (base == 'D') sset = A_BIT | G_BIT | T_BIT;
			else if (base == 'H') sset = A_BIT | C_BIT | T_BIT;
			else if (base == 'V') sset = A_BIT | C_BIT | G_BIT;
			else if (base == 'N') sset = A_BIT | C_BIT | G_BIT | T_BIT;
			else if (base == 'X') sset = A_BIT | C_BIT | G_BIT | T_BIT;

			/* total ambiguity */
			else if (base == '?') sset = A_BIT | C_BIT | G_BIT | T_BIT | O_BIT;

			/* deletion */
			else if (base == 'O') sset = O_BIT;
			else if (base == '-') {
				if (fifthstate == LVB_TRUE) {
					sset = O_BIT;
				}
				else {
					sset = A_BIT | C_BIT | G_BIT | T_BIT | O_BIT;
				}
			}
			lvb_assert(sset != 0U);
			enc_mat[i][j] = sset;
		}
    }
} /* end dna_makebin() */

void rowfree(Dataptr matrix)
/* free memory used for row strings and array of row strings in matrix,
 * and make the array of row title strings NULL;
 * or, if the array of row title strings is already NULL, do nothing */
{
    long i;	/* loop counter */

    if (matrix->row != NULL) {
    	for(i = 0; i < matrix->n; ++i) free(matrix->row[i]);
    	free(matrix->row);
    	matrix->row = NULL;
    }

} /* end rowfree() */


static long constchar(const Dataptr matrix, Lvb_bool *const togo, const Lvb_bool verbose)
/* Make sure matrix->m-element array togo is LVB_TRUE where matrix column
 * contains only one character state;
 * log details of new columns to ignore if verbose is LVB_TRUE.
 * scratch must point to the first element of an array of at least matrix->m
 * elements or arbitrary (even uninitialised) contents. It will be left with
 * arbitrary contents on return. */
{
    long k;		/* loop counter */
    long i;		/* loop counter */
    long n_columns = 0;

    /* discover variable columns */
    for (k = 0; k < matrix->m; ++k){
    	for (i = 1; i < matrix->n; ++i){
			if (matrix->row[i][k] != matrix->row[0][k]){
				togo[k] = LVB_TRUE;
				n_columns += 1;
				break;
			}
		}
    }

    if (verbose == LVB_TRUE){
    	printf("Ignoring constant columns\n");
    	if (n_columns == 0) printf("... none found.\n");
    	else logcut(togo, matrix->m);
    }
    return n_columns;
} /* end constchar() */

void matchange(Dataptr matrix, const Params rcstruct)
/* change and remove columns in matrix, partly in response to rcstruct,
 * verbosely or not according to value of verbose */
{
    Lvb_bool *togo;	/* LVB_TRUE where column must go */
    int n_columns_to_change = 0;
    /* Allocate memory: this will be free'd just before we return.
     * Dynamic allocation is used because otherwise, each of the
     * arrays would have to have MAX_M elements. That would either
     * limit the program too much, or causes massive waste of
     * address space. */

    togo = alloc(matrix->m * sizeof(Lvb_bool), "'togo' array");

    /* initialize all elements to LVB_FALSE ('don't ignore') */
    memset(togo, LVB_FALSE, matrix->m);

    n_columns_to_change = constchar(matrix, togo, rcstruct.verbose);	/* compuslory cut */

    /* N.B. a function to mark autapomorphic characters for cutting
     * could be called at this point. The effect would be more noticable
     * with unrealistically small test matrices than with real data */

    /* cut the cols as indicated, and crash verbosely if too few remain */
    if (n_columns_to_change != matrix->m) cutcols(matrix, togo, n_columns_to_change);	/* make changes to matrix */
    if (matrix->m < MIN_M)
    	crash("after constant columns are ignored, data matrix has\n"
    			"%ld columns, which is less than LVB's lower limit of\n"
    			"%ld columns.\n", matrix->m, MIN_M);
    else{
    	if (rcstruct.verbose == LVB_TRUE) printf("A total of %ld columns will be ignored\n", matrix->original_m - matrix->m);
    }

    /* free "local" dynamic heap memory */
    free(togo);

} /* end matchange() */

static void cutcols(Dataptr matrix, const Lvb_bool *const tocut, long n_columns_to_change)
/* remove columns in matrix for which the corresponding element of
matrix->m-element array tocut is LVB_TRUE, and update matrix->m;
return the number of columns cut */
{
    char **newrow;			/* rows of reduced matrix */
    long i;				/* loop counter */
    long k;				/* loop counter */
    long newk;				/* current column of reduced matrix */

    /* memory for new matrix row array */
    newrow = (char **) malloc(matrix->n * sizeof(char *));
    for (i = 0; i < matrix->n; ++i) newrow[i] =  (char*) malloc(sizeof(char) * (n_columns_to_change + 1));

    newk = 0;
    for (k = 0; k < matrix->m; ++k){	/* for every column */
		if (tocut[k] == LVB_TRUE){	/* keep this column */
			for (i = 0; i < matrix->n; ++i)	/* fill for ea. row */
				newrow[i][newk] = matrix->row[i][k];
			++newk;	/* fill next row next time */
		}
    }

    /* trap impossible condition */
    lvb_assert(newk == n_columns_to_change);

    /* terminate new row strings */
    for (i = 0; i < matrix->n; ++i) newrow[i][newk] = '\0';

    /* update matrix structure */
    rowfree(matrix);
    matrix->row = newrow;
    matrix->m = n_columns_to_change;

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

    memset(weight_arr, 0, m);

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
