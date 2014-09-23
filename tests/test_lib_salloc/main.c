/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include "lvb.h"

/* test of alloc() */

#define S1_LENGTH 120

static const char *s2 = "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	"@@@@@@@";	/* length 120 */


int main(void)
{
    char *s1;	/* test string */
    char *s3;	/* test string */
    long i;	/* loop counter */
    long j;	/* loop counter */

    lvb_initialize();

    for (i = 0; i < 50000; i++)	/* repeat to check heap seems OK */
    {
	s1 = alloc(S1_LENGTH + 1, "test string s1");
	for (j = 0; j < S1_LENGTH; j++) s1[j] = '@';
	s1[S1_LENGTH] = 0;
	lvb_assert(strcmp(s1, s2) == 0);
	s3 = alloc(strlen(s1)+1, "test string s3");
	strcpy(s3, s1);
	free(s1);
	lvb_assert(strcmp(s3, s2) == 0);
	free(s3);
    }

    printf("test passed\n");
    return 0;
}
