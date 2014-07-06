/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

sops.c - string manipulation functions

=cut

**********/

#include "lvb.h"

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

=head1 salloc - STRING ALLOCATION

=head2 SYNOPSIS

    char *salloc(long len, const char *msg);

=head2 DESCRIPTION

Allocates dynamic heap memory for a string. One may later free this
memory using the standard library function C<free()>. The new memory is
not initialized.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item len

Length of the string for which memory is being allocated, excluding the
terminating C<'\0'>. C<len> must be positive or zero.

=item msg

Pointer to the first character of a string describing the object being
allocated for. On failure, C<salloc()> will crash with message 'out of
memory: cannot allocate for ' followed by this string.

=back

=head2 RETURN

Pointer to the first character of the newly-allocated memory.

=cut

**********/

char *salloc(const long len, const char *const msg)
{
    char *s;	/* new space */

    lvb_assert(len >= 0);
    s = alloc(len + 1, msg);
    return s;
	
} /* end salloc() */

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
