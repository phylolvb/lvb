/* LVB
 * (c) Copyright 2003-2012 by Daniel Barker.
 * (c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
 * Permission is granted to copy and use this program provided that no fee is
 * charged for it and provided that this copyright notice is not removed. */

/* ********** getparam.c - get and set configurable parameters ********** */

#include "lvb.h"

static int get_default_seed(void)
/* return a default integer in the interval [0..MAX_SEED], obtained from the
 * system clock, or exit with an error message if the system time is
 * unavailable */
{
    time_t tim;			/* system time */
    unsigned long ul_seed;	/* seed value obtained from system time */

    tim = time(NULL);
    lvb_assert(tim != -1);
    ul_seed = (unsigned long) tim;
    ul_seed = ul_seed % (1UL + (unsigned long) MAX_SEED);
    lvb_assert(ul_seed <= MAX_SEED);
    return (int) ul_seed;

} /* end get_default_seed() */

static void read_line(char *buffer)
/* Read a line from standard input into memory starting at buffer.
 * Input will be truncated if its storage requirement exceeds
 * LVB_INPUTSTRING_SIZE. The buffer must have at least this allocation. */
{
    fgets(buffer, LVB_INPUTSTRING_SIZE, stdin);
    if (feof(stdin))
    {
	crash("end-of-file on standard input");
    }

} /* end read_line() */

void defaults_params(Params *const prms)
/* set seed in *prms to unacceptable value, and other parameters to their
 * defaults_params from lvb.h */
{
    prms->bootstraps = 0;	/* sensible default */

    /* meaningful value that is not user-configurable */
    prms->verbose = VERBOSE_OUTPUT;

    /* default value that will usually be used */
    prms->seed = get_default_seed();

    prms->p_file_name = "infile";

} /* end defaults_params() */

static void user_adjust(Params *prms)
/* Interact with the user and set the parameters they want */
{
    static char buffer[LVB_INPUTSTRING_SIZE];	/* input string */
    long lval;		/* input integer */

    /* matrix format */
  /*
//    printf("The data matrix must be in a file named 'infile', in PHYLIP 3.6 "
//     "format.\n"
//     "Are the aligned sequences INTERLEAVED or SEQUENTIAL?\n");
//    do
//    {
//        printf("Enter I for INTERLEAVED or S for SEQUENTIAL:\n");
//	read_line(buffer);
//    } while ((cistrcmp(buffer, "S\n") != 0) && (cistrcmp(buffer, "I\n") != 0));
//    switch (toupper(buffer[0]))
//    {
//    case 'I':
//	prms->interleaved = LVB_TRUE;
//        break;
//    case 'S':
//	prms->interleaved = LVB_FALSE;
//        break;
//    }
    prms->interleaved = LVB_FALSE;*/
    
    /* treatment of gaps */
    printf("\nPlease choose treatment of gaps represented by '-' in the data "
         "matrix.\n"
         "These should usually be treated as UNKNOWN.\n");
    do
    {
        printf("Enter U for UNKNOWN or F for FIFTH STATE:\n");
        read_line(buffer);
    } while ((cistrcmp(buffer, "U\n") != 0) && (cistrcmp(buffer, "F\n") != 0));
    switch (toupper(buffer[0]))
    {
    case 'U':
	prms->fifthstate = LVB_FALSE;
	break;
    case 'F':
	prms->fifthstate = LVB_TRUE;
	break;
    }

    /* cooling schedule */
    printf("\nPlease select a cooling schedule. The schedule chosen\n" 
     "will affect the quality and speed of the simulated annealing\n"
     "search. The GEOMETRIC schedule will take significantly less time,\n"
     "but may produce lower quality results. The LINEAR schedule may\n"
     "produce higher quality results, at the cost of increased runtime.\n"
     "Currently, the DEFAULT is the GEOMETRIC schedule.\n");
    do
    {
        printf("Enter G for GEOMETRIC or L for LINEAR\n"
            "or press RETURN for default:\n");
        read_line(buffer);
        if ((strcmp(buffer, "\n") == 0))
        {
            prms->cooling_schedule = 0;
            break;
        }
    } while ((cistrcmp(buffer, "G\n") != 0) && (cistrcmp(buffer, "L\n") != 0));
    switch (toupper(buffer[0]))
    {
    case 'G':
	prms->cooling_schedule = 0;
	break;
    case 'L':
	prms->cooling_schedule = 1;
	break;
    }

    /* seed */
    printf("\nPlease specify a random number seed, or use default.\n"
     "The default varies (it is taken from the system clock) and is "
     "recommended.\n");
    do
    {
        printf("Enter an integer in the range 0 to %ld inclusive,\n"
         "or press RETURN for default:\n", (long) MAX_SEED);
        read_line(buffer);
	if ((strcmp(buffer, "\n") == 0))
	    lval = prms->seed;
	else
	    lval = strtol(buffer, NULL, 10);
    } while ((lval < 0L) || (lval > MAX_SEED));
    prms->seed = (int) lval;

    /* bootstraps */
    printf("\nPlease indicate whether you require bootstrapping.");
    do
    {
	printf("\nEnter the number of bootstrap replicates required\n");
        printf("as an integer in the range 1 to %ld inclusive,\n"
	    "or press RETURN for no bootstrapping:\n", (long) MAX_BOOTSTRAPS);
	read_line(buffer);
	if ((strcmp(buffer, "\n") == 0))
	    lval = 0;
	else
	    lval = strtol(buffer, NULL, 10);
    } while ((lval < 0L) || (lval > MAX_BOOTSTRAPS));
    prms->bootstraps = lval;

    printf("\n");

}

void getparam(Params *prms)
/* Get configuration parameters. This function fills *prms with
 * run-time configuration parameters */
{
	defaults_params(prms);
    user_adjust(prms);

} /* end getparam() */
