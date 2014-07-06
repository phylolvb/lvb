/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/**********

=head1 NAME

err.c - error functions

=cut

**********/

#include "lvb.h"

/**********

=head1 crash - CRASH VERBOSELY

=head2 SYNOPSIS

    void crash(const char *fmt, ...);

=head2 DESCRIPTION

Prints a message consisting of 'FATAL ERROR: ' followed by a
user-supplied, C<vfprintf>-acceptable message. Then calls C<cleanup()>
and exits abnormally.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item fmt

Pointer to the first character in a format string acceptable to
C<vfprintf()>.

=item ...

Any parameters whose presence is dictated by C<fmt>.

=back

=cut

**********/

void crash(const char *const fmt, ...)
{
	const char *const warning = "\nFATAL ERROR";	/* dire warning */
	va_list args;					/* arguments */

	va_start(args, fmt);
	printf("%s: ", warning);
	vprintf(fmt, args);
	va_end(args);
	printf("\n");

	cleanup();
	exit(EXIT_FAILURE);

}	/* end crash() */

void lvb_assertion_fail(const char *test, const char *file, int line)
/* Log dire warning followed by message of form "assertion failed at
 * '<file>' line <line>: <test>", and exit abnormally. This function
 * should only be called through the lvb_assert() macro. */
{
    crash("assertion failed at '%s' line %d: %s", file, line, test);
}

/**********

=head1 scream - log dire warning

=head2 SYNOPSIS

    void scream(const char *format, ...);

=head2 DESCRIPTION

Prints a message consisting of 'ERROR: ' followed by a user-supplied,
C<vfprintf>-acceptable message.

=head2 PARAMETERS

=head3 INPUT

=over 4

=item format

Pointer to the first character in a format string acceptable to
C<vfprintf()>.

=item ...

Any parameters whose presence is dictated by C<format>.

=back

=cut

**********/

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
