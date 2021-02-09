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

/* ========== Sops.c - string manipulation functions ========== */

#include "LVB.h"

enum { SAME, DIFFERENT };	/* strings same or different */

/**********

=head1 cistrcmp - CASE-INDEPENDENT STRING COMPARISON

=head2 SYNOPSIS

    long cistrcmp(const char *s1, const char *s2);

=head2 DESCRIPTION

Compares two strings in a case-independent manner. The strings may have
the same or different length. (Strings of nonequal length are always
considered different from each other.)

The strings may overlap in memory or occupy the same memory.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item s1

Pointer to the first text character of one of the strings to be
compared.

=item s2

Pointer to the first text character of the other string to be compared.

=back

=head2 RETURN

Returns 0 if the string are the same or nonzero if different.

=cut

**********/

long cistrcmp(const char *const s1, const char *const s2)
{
    size_t i;		/* loop counter */
    size_t len1;	/* length of s1 */
    size_t len2;	/* length of s2 */
    int character_1;	/* current character of s1 */
    int character_2;	/* current character of s2 */

    len1 = strlen(s1);
    len2 = strlen(s2);

    if (len1 != len2)	/* can't be identical */
	return DIFFERENT;

    for (i = 0; i < len1; i++)
    {
        character_1 = tolower(s1[i]);
        character_2 = tolower(s2[i]);
	if (character_1 != character_2)
	    return DIFFERENT;
    }
    return SAME;

} /* end cistrcmp() */

/**********

=head1 nextnonwspc - FIND NEXT NON-WHITE-SPACE CHARACTER

=head2 SYNOPSIS

    char *nextnonwspc(const char *string);

=head2 DESCRIPTION

Finds the next non-white-space text character in a string, if any.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item string

Pointer to the first character in the string in question.

=back

=head2 RETURN

Pointer to the first non-white-space character in the string, or
C<NULL> if all characters in the string are white space.

=cut

**********/

char *nextnonwspc(const char *string)
{
    while (isspace(*string))
	string++;
    if (*string)
	return (char *) string;
    else
	return NULL;

} /* end nextnonwspc() */


/**********

=head1 supper - CONVERT TO UPPER CASE

=head2 SYNOPSIS

    char *supper(char *s);

=head2 DESCRIPTION

Convert a string to upper case.

=head2 PARAMETERS

=head3 INPUT/OUTPUT

=over 4

=item s

Pointer to first text character in the string to be converted.

=back

=head2 RETURN

Returns C<s>.

=cut

**********/

char *supper(char *const s)
{
    int character;		/* current character in uppercase */
    char *elementptr = s;	/* pointer to current character */

    while (*elementptr)
    {
        character = toupper(*elementptr);
	*elementptr = (char) character;
	elementptr++;
    }

    return s;

} /* end supper() */