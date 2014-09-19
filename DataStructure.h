

/* matrix and associated information */
typedef struct data
{
    char **row;			/* array of row strings */
    char **rowtitle;	/* array of row title strings */
    long m;				/* number of columns */
    long original_m;	/* number of columns read from matrix*/
    long n;				/* number of rows */
    long nbranches; 	/* number of possible braches */
    int n_threads_getplen;  /* number of possible threads in getplen function */
    int n_threads_try_combination;  /* number of possible threads try combinations trees */
    int n_slice_size_getplen;  /* slice size in getplen function, usually m/n_threads_getplen  */
} *Dataptr, DataStructure;

