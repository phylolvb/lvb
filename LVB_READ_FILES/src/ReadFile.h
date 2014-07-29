/*
 * ReadFile.h
 *
 *  Created on: Jun 24, 2014
 *      Author: mmp
 */
#include "CReadFiles.h"
#include "../../DataStructure.h"
#include <stdio.h>
#include <string.h>
using namespace std;

#ifndef READFILE_H_
#define READFILE_H_

extern "C" void read_file(char *file_name, int n_file_type, Dataptr p_lvbmat);
extern "C" void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);


void read_file(char *file_name, int n_file_type, DataStructure *p_lvbmat);
void phylip_mat_dims_in_external(char *file_name, int n_file_type, long *species_ptr, long *sites_ptr, int *max_length_name);
void free_lvbmat_structure(DataStructure *p_lvbmat);
long brcnt(long n) { return (n << 1) - 3; }; /* return number of branches in unrooted binary tree structure containing n tips */

#endif /* READFILE_H_ */
