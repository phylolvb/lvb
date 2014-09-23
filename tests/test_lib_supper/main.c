/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include "lvb.h"

/* test for supper() */

int main(void)
{
    /* test strings */
    char s[] = "ASDFGHJZJKL\t \n!,;:@!-X";
    char t[] = "asdfgHjzjkl\t \n!,;:@!-x";
    char u[] = "";

    char *copy;		/* copy of string */
    char *value;	/* copy of pointer value */

    lvb_initialize();

    copy = alloc(strlen(s), "copy of the string");
    strcpy(copy, s);
    value = supper(s);
    lvb_assert(value == s);
    lvb_assert(strcmp(copy, s) == 0);
    free(copy);
    
    value = supper(t);
    lvb_assert(value == t);
    lvb_assert(strcmp(s, t) == 0);

    value = supper(u);
    lvb_assert(value == u);
    lvb_assert(strcmp(u, "") == 0);

    printf("test passed\n");
    return 0;
}
