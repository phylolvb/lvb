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

/* ========== fops.c - file operations ========== */

#include "lvb.h"

/**********

=head1 CheckFileOpening - OPEN FILE WITH ERROR CHECKING

=head2 SYNOPSIS

    FILE *CheckFileOpening(const char *const nam, const char *const mod);

=head2 DESCRIPTION

Open a file, or crash verbosely on failure.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item nam

Pointer to first byte of a string giving the name of the file to open.

=item mod

Mode to open the file with. Values are as for the standard C library
function C<fopen>.

=back

=head2 RETURN

Returns a pointer to the newly opened file structure.

=cut

**********/

FILE *CheckFileOpening(const char *const nam, const char *const mod)
{
    FILE *fp;	/* file */

    fp = fopen(nam, mod);
    if (fp == NULL)
    {
	if (strcmp(mod, "w") == 0)
	    CrashVerbosely("cannot create file '%s'", nam);
	else if (strcmp(mod, "r") == 0)
	    CrashVerbosely("cannot open file '%s' for reading", nam);
	else if (strcmp(mod, "a") == 0)
	    CrashVerbosely("cannot open file '%s' for appending to", nam);
	else	/* rare mode */
	    CrashVerbosely("cannot open file '%s' with mode '%s'", nam, mod);
    }

    return fp;

} /* end CheckFileOpening() */

/**********

=head1 CheckFileClosure - CLOSE FILE WITH ERROR CHECKING

=head2 SYNOPSIS

    void CheckFileClosure(FILE *const fp, const char *const fnam);

=head2 DESCRIPTION

Close a file, crashing verbosely on file error or failure. Or, if passed
a NULL file pointer, do nothing.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item fnam

Pointer to first byte of a string giving the name of the file to close.

=back

=head3 INOUT

=over 4

=item fp

Pointer to the file to check for errors and close, or NULL.

=back

=head2 RETURN

Returns a pointer to the newly opened file structure.

=cut

**********/

void CheckFileClosure(FILE *const fp, const char *const fnam)
{
    if (fp != NULL)
    {
	if (ferror(fp) != 0)
	    CrashVerbosely("file error on file '%s'", fnam);
	if (fclose(fp) != 0)
	    CrashVerbosely("cannot close file '%s'", fnam);
    }

} /* end CheckFileClosure() */
