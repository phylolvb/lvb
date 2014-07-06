/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

#include "lvb.h"

/* test for nextnonwspc() */

int main(void)
{
    const char *p1 = "";
    const char *p2 = "\b\a!\"$%^&*()_+NO_SPACE+IN$HERE!@#~'";
    const char *p3 = "\v\t\r\n \f\aHello Goodbye \n";

    lvb_initialize();

    lvb_assert(nextnonwspc(p1) == NULL);
    lvb_assert(nextnonwspc(p2) == p2);
    lvb_assert(nextnonwspc(p3) == p3 + 6);

    printf("test passed\n");
    return 0;
}
