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

/* ========== FileOperations.c - file operations ========== */

#include "FileOperations.h"

/**********

=head1 file_exists - CHECK A FILE EXISTS

=head2 SYNOPSIS

    Lvb_bool file_exists(const char *const nam);

=head2 DESCRIPTION

Checks whether a file may be opened for reading. This actually tests
both the presence of the file and its permissions.

=head2 PARAMETERS

=head3 INPUT

=over 4     

=item nam

Pointer to first byte of a string giving the name of the file to check.

=back

=head2 RETURN

Returns LVB_TRUE if the file exists and may be opened for reading, LVB_FALSE
otherwise.

=cut

**********/

Lvb_bool file_exists(const char *const nam)
{
    Lvb_bool val;		/* return value */
    FILE *fp = fopen(nam, "r");	/* store file ptr so can close it */
	
    if (fp == NULL)
	val = LVB_FALSE;
    else
    {
	fclose(fp);
	val = LVB_TRUE;
    }

    return val;

} /* end */

/**********

=head1 clnopen - OPEN FILE WITH ERROR CHECKING

=head2 SYNOPSIS

    FILE *clnopen(const char *const nam, const char *const mod);

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

FILE *clnopen(const char *const nam, const char *const mod)
{
    FILE *fp;	/* file */

    fp = fopen(nam, mod);
    if (fp == NULL)
    {
	if (strcmp(mod, "w") == 0)
	    crash("cannot create file '%s'", nam);
	else if (strcmp(mod, "r") == 0)
	    crash("cannot open file '%s' for reading", nam);
	else if (strcmp(mod, "a") == 0)
	    crash("cannot open file '%s' for appending to", nam);
	else	/* rare mode */
	    crash("cannot open file '%s' with mode '%s'", nam, mod);
    }

    return fp;

} /* end clnopen() */

/**********

=head1 clnclose - CLOSE FILE WITH ERROR CHECKING

=head2 SYNOPSIS

    void clnclose(FILE *const fp, const char *const fnam);

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

void clnclose(FILE *const fp, const char *const fnam)
{
    if (fp != NULL)
    {
	if (ferror(fp) != 0)
	    crash("file error on file '%s'", fnam);
	if (fclose(fp) != 0)
	    crash("cannot close file '%s'", fnam);
    }

} /* end clnclose() */

/**********

=head1 clnremove - DELETE FILE WITH ERROR CHECKING

=head2 SYNOPSIS

    void clnremove(const char *const fnam);

=head2 DESCRIPTION

Delete a file, logging a strong warning on failure.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item fnam

Pointer to first byte of a string giving the name of the file to delete.

=back

=head2 RETURN

None.

=head2 BUGS

For portability, the file should not be open on calling C<clnremove>.

=cut

**********/

void clnremove(const char *const fnam)
{
    if (remove(fnam) != 0)
	scream("cannot delete file '%s'", fnam);

} /* end clnremove */

/**********

=head1 f2str - READ FILE INTO STRING

=head2 SYNOPSIS

    char *f2str(FILE *const stream);

=head2 DESCRIPTION

Read contents of a file into a new string. The text read is terminated by
newline if not already present, and finally by '\0'. The file is treated
as text. One may deallocate memory for the new string using the standard
library function C<free>.

=head2 PARAMETERS

=head3 INOUT

=over 4

=item stream

Pointer to file, which must be open in a readable mode and positioned at
the start.

=back

=head2 RETURN

Pointer to first byte in string containing the contents of the file when
read as text.

=head2 BUGS

Does not know how to interpret newlines if the file has a different
representation of newline to that used by the current system.

Will not work correctly if the file size alters during the call to C<f2str>.

Does not work correctly if the stream is positioned anywhere other than
at the start of the file on calling C<f2str>.

=cut

**********/

char *f2str(FILE *const stream)
{
    char *input;				/* input string */
    unsigned long inbytes;			/* bytes for string */
    unsigned long off;				/* position in string */
    unsigned long offmax = 0UL;			/* length of string */
    const unsigned long maxom = ULONG_MAX - 3;	/* maximum initial offmax */

    /* calculate file size and allocate appropriately */
    while (getc(stream) != EOF)
    {
	offmax++;
	if (offmax >= maxom)	/* crash while value has meaning */
	    crash("input is too long");
    }
    if (ferror(stream) != 0)
	crash("file error on reading file");
    inbytes = offmax + 2UL;	/* '\0', possible '\n' */
    input = (char *) alloc(inbytes, "input");

    /* get string */
    rewind(stream);
    for (off = 0; off < offmax; off++)
	input[off] = (char) getc(stream);
    if (ferror(stream) != 0)
	crash("file error on reading file");

    /* terminate string, also adding newline at end if not present */
    if (input[off-1] != '\n')
    {
	input[off] = '\n';
	off++;
    }
    input[off] = '\0';

    return input;

} /* end f2str() */
