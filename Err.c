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

/* ========== err.c - error functions ========== */

#include "Lvb.h"

void crash(const char *const fmt, ...)
{
	const char *const warning = "\nFATAL ERROR";	/* dire warning */
	va_list args;					/* arguments */

	va_start(args, fmt);
	printf("%s: ", warning);
	vprintf(fmt, args);
	va_end(args);
	printf("\n");

  #ifdef NP_Implementation
	cleanup();
  #endif // NP_Implementation #54

  #ifdef MPI_Implementation
  int n_error_code = 1;
	MPI_Abort(MPI_COMM_WORLD, n_error_code);
  #endif // MPI_Implementation #58

	exit(1);
}	/* end crash() */
void lvb_assertion_fail(const char *test, const char *file, int line)
/* Log dire warning followed by message of form "assertion failed at
 * '<file>' line <line>: <test>", and exit abnormally. This function
 * should only be called through the lvb_assert() macro. */
{
    crash("assertion failed at '%s' line %d: %s", file, line, test);
}

void scream(const char *const format, ...)
/* log a dire warning, partly composed of vprintf-acceptable user-supplied
 * message */
{
	const char *const warning = "ERROR";
	va_list args;		/* supplied message */

	va_start(args, format);
	printf("%s: ", warning);
	vprintf(format, args);
	va_end(args);
	printf("\n");

	/* flush standard output so the warning is immediately visible */
	if (fflush(stdout) == EOF)
		crash("write error on log");	/* may not work! */

}	/* end scream() */
