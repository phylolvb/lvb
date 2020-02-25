/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== Random_Number_Generator.c - random number generation ========== */

/* exactly as suppplied by EPCC, except rcsid added, printf() then
 * exit() replaced with crash() (and braces made redundant removed),
 * float replaced with double, declaration for rstart added, "lvb.h"
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

#include "lvb.h"

#ifndef NP_Implementation
#include "Parallel_Checkpointing.h"
#endif

static void rstart(int i, int j, int k, int l);

/*
 *	Global variables for rstart & uni
 */

/*const int NUMBER_UNI_AVAILABLE = 10000;
const int NUMBER_UNI_BUFFER = 2000;
static double uni_values[NUMBER_UNI_AVAILABLE + NUMBER_UNI_BUFFER];*/
#define NUMBER_MAX_UNI	98
static double uni_u[NUMBER_MAX_UNI];
static double uni_c;
const double uni_cd = 7654321.0 / 16777216.0;
const double uni_cm = 16777213.0 / 16777216.0;
static int uni_ui, uni_uj;
static Lvb_bool rinit_called = LVB_FALSE;	/* added - DB */


#ifndef NP_Implementation
/* each structure has the number of bytes to read in the first position in the file */
/* the last one is a checksum unsigned long, it not summed in the structure */
unsigned long checkpoint_uni(FILE *fp)
{
	unsigned long n_bytes_to_write = sizeof(uni_u) + sizeof(double) + sizeof(int) + sizeof(int) + sizeof(Lvb_bool) + sizeof(unsigned short);
	unsigned long checksum = 0;
	unsigned short type_block = STATE_BLOCK_UNI;
	fwrite(&n_bytes_to_write, sizeof(n_bytes_to_write), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_write), (unsigned char *) &n_bytes_to_write, checksum);
	fwrite(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
	fwrite(uni_u, sizeof(uni_u), 1, fp);
    for (int i = 0; i < NUMBER_MAX_UNI; i ++) checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) &uni_u[i], checksum);
    fwrite(&uni_c, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(uni_c), (unsigned char *) &uni_c, checksum);
    fwrite(&uni_ui, sizeof(int), 1, fp); checksum = CalculateBlockCRC32(sizeof(uni_ui), (unsigned char *) &uni_ui, checksum);
    fwrite(&uni_uj, sizeof(int), 1, fp); checksum = CalculateBlockCRC32(sizeof(uni_uj), (unsigned char *) &uni_uj, checksum);
    fwrite(&rinit_called, sizeof(Lvb_bool), 1, fp); checksum = CalculateBlockCRC32(sizeof(rinit_called), (unsigned char *) &rinit_called, checksum);
    fwrite(&checksum, sizeof(unsigned long), 1, fp);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(fflush(fp) == 0);
//   print_information_checkpoint("Save data uni", n_bytes_to_write, checksum);
    return checksum;
}


unsigned long restore_uni(FILE *fp)
{
	unsigned long n_bytes_to_write = sizeof(uni_u) + sizeof(double) + sizeof(int) + sizeof(int) + sizeof(Lvb_bool) + sizeof(unsigned short), n_bytes_to_read = 0;
	unsigned long checksum = 0, checksum_read, n_read_values;
	unsigned short type_block;
	n_read_values = fread(&n_bytes_to_read, sizeof(n_bytes_to_read), 1, fp); checksum = CalculateBlockCRC32(sizeof(n_bytes_to_read), (unsigned char *) &n_bytes_to_read, checksum);
	n_read_values = fread(&type_block, sizeof(type_block), 1, fp); checksum = CalculateBlockCRC32(sizeof(type_block), (unsigned char *) &type_block, checksum);
	n_read_values = fread(uni_u, sizeof(uni_u), 1, fp);
    for (int i = 0; i < NUMBER_MAX_UNI; i ++) checksum = CalculateBlockCRC32(sizeof(double), (unsigned char *) &uni_u[i], checksum);
    n_read_values = fread(&uni_c, sizeof(double), 1, fp); checksum = CalculateBlockCRC32(sizeof(uni_c), (unsigned char *) &uni_c, checksum);
    n_read_values = fread(&uni_ui, sizeof(int), 1, fp); checksum = CalculateBlockCRC32(sizeof(uni_ui), (unsigned char *) &uni_ui, checksum);
    n_read_values = fread(&uni_uj, sizeof(int), 1, fp); checksum = CalculateBlockCRC32(sizeof(uni_uj), (unsigned char *) &uni_uj, checksum);
    n_read_values = fread(&rinit_called, sizeof(Lvb_bool), 1, fp); checksum = CalculateBlockCRC32(sizeof(rinit_called), (unsigned char *) &rinit_called, checksum);
    n_read_values = fread(&checksum_read, sizeof(unsigned long), 1, fp);
    lvb_assert(n_read_values == 1);
    lvb_assert(ferror(fp) == 0);
    lvb_assert(n_bytes_to_read == n_bytes_to_write);
    lvb_assert(checksum_read == checksum);
//    print_information_checkpoint("Read data uni", n_bytes_to_write, checksum);
    return checksum;
}

#endif

double uni(void)
{
	double luni;			/* local variable for uni */

	lvb_assert(rinit_called != LVB_FALSE);	/* added - DB */
	luni = uni_u[uni_ui] - uni_u[uni_uj];
	if (luni < 0.0) luni += 1.0;
	uni_u[uni_ui] = luni;
	if (--uni_ui == 0) uni_ui = NUMBER_MAX_UNI - 1;
	if (--uni_uj == 0) uni_uj = NUMBER_MAX_UNI - 1;
	if ((uni_c -= uni_cd) < 0.0) uni_c += uni_cm;
	if ((luni -= uni_c) < 0.0) luni += 1.0;
	return (double) luni;
}

static void rstart(int i, int j, int k, int l)
{
	int ii, jj, m;
	double s, t;

	for (ii = 1; ii < NUMBER_MAX_UNI; ii++)
	{
		s = 0.0;
		t = 0.5;
		for (jj = 1; jj <= 24; jj++) {
			m = ((i*j % 179) * k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l+1) % 169;
			if (l*m % 64 >= 32)
				s += t;
			t *= 0.5;
		}
		uni_u[ii] = s;
	}
	uni_c  = 362436.0   / 16777216.0;
	uni_ui = 97;	/*  There is a bug in the original Fortran version */
	uni_uj = 33;	/*  of UNI -- i and j should be SAVEd in UNI()     */
}

void rinit(int ijkl)
{
	int i, j, k, l, ij, kl;

	/* tell later calls to uni() that rinit() has been called - DB */
	rinit_called = LVB_TRUE;

	/* check double type is suitable */
	if (DBL_MANT_DIG < 24) crash("FP type unsuitable for uni() random no. generator\n"); /* too small */

	/* check ijkl is within range */
	if( (ijkl < 0) || (ijkl > 900000000) ) crash("rinit: ijkl = %d -- out of range\n\n", ijkl);

/*        printf("rinit: seed_ijkl = %d\n", ijkl); */

	/* decompose the long integer into the the equivalent four
	 * integers for rstart. This should be a 1-1 mapping
 	 *	ijkl <--> (i, j, k, l)
	 * though not quite all of the possible sets of (i, j, k, l)
	 * can be produced.
	 */

	ij = ijkl/30082;
	kl = ijkl - (30082 * ij);

	i = ((ij/177) % 177) + 2;
	j = (ij % 177) + 2;
	k = ((kl/169) % 178) + 1;
	l = kl % 169;

	if( (i <= 0) || (i > 178) ) crash("rinit: i = %d -- out of range\n\n", i);
	if( (j <= 0) || (j > 178) ) crash("rinit: j = %d -- out of range\n\n", j);
	if( (k <= 0) || (k > 178) ) crash("rinit: k = %d -- out of range\n\n", k);
	if( (l < 0) || (l > 168) ) crash("rinit: l = %d -- out of range\n\n", l);
	if (i == 1 && j == 1 && k == 1) crash("rinit: 1 1 1 not allowed for 1st 3 seeds\n\n");

/*        printf("rinit: initialising RNG via rstart(%d, %d, %d, %d)\n", i, j, k, l); */
	rstart(i, j, k, l);
}

long randpint(const long upper)
{
    double frand;	/* random real */
    double fupper;	/* upper limit */
    long rand;		/* return value */

    lvb_assert(upper >= 0);

    fupper = (double) upper;
    frand = uni();
    frand = frand * fupper;		/* scale to right range */
    rand = (long) (frand + 0.5);	/* round to nearest integer */

    /* guard against arithmetic inaccuracy */
    if (rand < 0)
	rand = 0;
    else if (rand > upper)
	rand = upper;

    return rand;

} /* end randpint() */
