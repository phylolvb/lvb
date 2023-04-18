/* LVB

(c) Copyright 2003-2012 by Daniel Barker.  
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.  
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.  
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.  
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.  
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== RandomNumberGenerator.c - random number generation ========== */

/* exactly as suppplied by EPCC, except rcsid added, printf() then
 * exit() replaced with crash() (and braces made redundant removed),
 * float replaced with double, declaration for rstart added, "LVB.h"
 * included, global variables and rstart() made static, comment after
 * declaration of uni_u shortened to fit 80 columns after the addition
 * of the word 'static', mechanism added for uni() to crash with an
 * assertion failure if rinit() has not yet been called, and check on
 * DBL_MANT_DIG added to rinit(). This last has to be done at run-time
 * for portability, since DBL_MANT_DIG (from float.h) might not be a
 * constant.
 */

/*
 *	C version of Marsaglia's UNI random number generator
 *	More or less transliterated from the Fortran -- with 1 bug fix
 *	Hence horrible style
 *
 *	Features:
 *		ANSI C
 *		not callable from Fortran (yet)
 */

#include "RandomNumberGenerator.h"

static void rstart(int i, int j, int k, int l);

/*
 *	Global variables for rstart & uni
 */

/*const int NUMBER_UNI_AVAILABLE = 10000;
const int NUMBER_UNI_BUFFER = 2000;
static double uni_values[NUMBER_UNI_AVAILABLE + NUMBER_UNI_BUFFER];*/
static double uni_u[98]; /* Was U(97) in Fortran version */
static double uni_c;
const double uni_cd = 7654321.0 / 16777216.0;
const double uni_cm = 16777213.0 / 16777216.0;
static int uni_ui, uni_uj;
static Lvb_bool rinit_called = LVB_FALSE; /* added - DB */

double uni(void)
{
	double luni; /* local variable for uni */

	lvb_assert(rinit_called != LVB_FALSE); /* added - DB */
	luni = uni_u[uni_ui] - uni_u[uni_uj];
	if (luni < 0.0)
		luni += 1.0;
	uni_u[uni_ui] = luni;
	if (--uni_ui == 0)
		uni_ui = 97;
	if (--uni_uj == 0)
		uni_uj = 97;
	if ((uni_c -= uni_cd) < 0.0)
		uni_c += uni_cm;
	if ((luni -= uni_c) < 0.0)
		luni += 1.0;
	return (double)luni;
}

static void rstart(int i, int j, int k, int l)
{
	int ii, jj, m;
	double s, t;

	for (ii = 1; ii <= 97; ii++)
	{
		s = 0.0;
		t = 0.5;
		for (jj = 1; jj <= 24; jj++)
		{
			m = ((i * j % 179) * k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53 * l + 1) % 169;
			if (l * m % 64 >= 32)
				s += t;
			t *= 0.5;
		}
		uni_u[ii] = s;
	}
	uni_c = 362436.0 / 16777216.0;
	uni_ui = 97; /*  There is a bug in the original Fortran version */
	uni_uj = 33; /*  of UNI -- i and j should be SAVEd in UNI()     */
}

/* ~rinit: this takes a single integer in the range
		0 <= ijkl <= 900 000 000
	and produces the four smaller integers needed for rstart. It is
 *	based on the ideas contained in the RMARIN subroutine in
 *		F. James, "A Review of Pseudorandom Number Generators",
 *			Comp. Phys. Commun. Oct 1990, p.340
 *	To reduce the modifications to the existing code, rinit now
 *	takes the role of a preprocessor for rstart.
 *
 *	This is useful for the parallel version of the code as James
 *	states that any integer ijkl will produce a statistically
 *	independent sequence of random numbers.
 *
 *     Very funny. If that statement was worth anything he would have provided
 *     a proof to go with it. spb 12/12/90
 */

void rinit(int ijkl)
{
	int i, j, k, l, ij, kl;

	/* tell later calls to uni() that rinit() has been called - DB */
	rinit_called = LVB_TRUE;

	/* check double type is suitable */
	if (DBL_MANT_DIG < 24)
		crash("FP type unsuitable for uni() random no. generator\n"); /* too small */

	/* check ijkl is within range */
	if ((ijkl < 0) || (ijkl > 900000000))
		crash("rinit: ijkl = %d -- out of range\n\n", ijkl);

	/*        printf("rinit: seed_ijkl = %d\n", ijkl); */

	/* decompose the long integer into the the equivalent four
	 * integers for rstart. This should be a 1-1 mapping
	 *	ijkl <--> (i, j, k, l)
	 * though not quite all of the possible sets of (i, j, k, l)
	 * can be produced.
	 */

	ij = ijkl / 30082;
	kl = ijkl - (30082 * ij);

	i = ((ij / 177) % 177) + 2;
	j = (ij % 177) + 2;
	k = ((kl / 169) % 178) + 1;
	l = kl % 169;

	if ((i <= 0) || (i > 178))
		crash("rinit: i = %d -- out of range\n\n", i);
	if ((j <= 0) || (j > 178))
		crash("rinit: j = %d -- out of range\n\n", j);
	if ((k <= 0) || (k > 178))
		crash("rinit: k = %d -- out of range\n\n", k);
	if ((l < 0) || (l > 168))
		crash("rinit: l = %d -- out of range\n\n", l);
	if (i == 1 && j == 1 && k == 1)
		crash("rinit: 1 1 1 not allowed for 1st 3 seeds\n\n");

	/*        printf("rinit: initialising RNG via rstart(%d, %d, %d, %d)\n",
					i, j, k, l); */

	rstart(i, j, k, l);
}

/**********

=head1 randpint - GET RANDOM POSITIVE INTEGER

=head2 SYNOPSIS

	long randpint(long upper);

=head2 DESCRIPTION

Returns a positive integer with a user-specified upper limit.

Before calling C<randpint()>, ensure the random number generator has
been initialized by calling C<rinit()>.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item upper

The random number will be in the range [0..C<upper>] inclusive.

=back

=head3 RETURN

A random integer in the interval [0..C<upper>].

=cut

**********/

long randpint(const long upper)
{
	double frand;  /* random real */
	double fupper; /* upper limit */
	long rand;	   /* return value */

	lvb_assert(upper >= 0);

	fupper = (double)upper;
	frand = uni();
	frand = frand * fupper;		/* scale to right range */
	rand = (long)(frand + 0.5); /* round to nearest integer */

	/* guard against arithmetic inaccuracy */
	if (rand < 0)
		rand = 0;
	else if (rand > upper)
		rand = upper;

	return rand;

} /* end randpint() */
