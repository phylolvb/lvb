/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include <lvb.h>

/* positive test for cistrcmp() */

int main(void)
{
    /* strings that are the same, considered case-insensitively */
    static char s1[] = "`1234567890-=!\"$%^&*("
     ")_+qwertyuiopasdfghjklzxcvbnm[]{};'#:@~,./<>? \t\n\r";
    static char s2[] = "`1234567890-=!\"$%^&*("
     ")_+QWERTYUIOPASDFGHJKLZXCVBNM[]{};'#:@~,./<>? \t\n\r";

    /* string that is almost the same but actually different */
    static char s3[] = "`1234567890-=!\"$%^&*("
     ")_+QWERTYUIOPASDFGHJKLZXCVBNN[]{};'#:@~,./<>? \t\n\r";

    lvb_initialize();

    lvb_assert(cistrcmp(s1, s1) == 0);
    lvb_assert(cistrcmp(s3 + 5, s3) != 0);
    lvb_assert(cistrcmp(s1, s2) == 0);
    lvb_assert(cistrcmp(s2, s3) != 0);

    printf("test passed\n");
    return 0;
}
