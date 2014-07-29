

/* matrix and associated information */
typedef struct data
{
    char **row;			/* array of row strings */
    char **rowtitle;	/* array of row title strings */
    long m;				/* number of columns */
    long original_m;	/* number of columns read from matrix*/
    long n;				/* number of rows */
    long nbranches; 	/* number of possible braches */
} *Dataptr, DataStructure;

