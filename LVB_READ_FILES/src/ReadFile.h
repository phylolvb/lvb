/*
 * ReadFile.h
 *
 *  Created on: Jun 24, 2014
 *      Author: mmp
 */
#include "CReadFiles.h"
#include <stdio.h>
#include <string.h>
using namespace std;

#ifndef READFILE_H_
#define READFILE_H_

/* matrix and associated information */
typedef struct data
{
    char **row;		/* array of row strings */
    long m;		/* number of columns */
    long n;		/* number of rows */
    char **rowtitle;	/* array of row title strings */
} *Dataptr, DataStructure;

// nm -D libLVB_READ_FILES_LD.so | grep " T "
/// to export
extern "C" void read_file(char *file_name, DataStructure *p_lvbmat);
extern "C" void phylip_mat_dims_in_external(char *file_name, long *species_ptr, long *sites_ptr);


void read_file(char *file_name, DataStructure *p_lvbmat);
void phylip_mat_dims_in_external(char *file_name, long *species_ptr, long *sites_ptr);
void free_lvbmat_structure(DataStructure *p_lvbmat);


#endif /* READFILE_H_ */
