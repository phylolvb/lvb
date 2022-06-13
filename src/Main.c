/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and
Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2020 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/* ========== Main.c ========== */

#include "Admin.h"
#include "Clock.h"
#include "DataOperations.h"
#include "Hash.h"
#include "Print.h"
#include "LogFile.h"
#include "LVB.h"
#include "MemoryOperations.h"
#include "SearchParameters.h"
#include "Solve.c"
#include "Solve.h"
#include "Verbose.h"

int main(int argc, char **argv)
{
	Dataptr MSA;	/* data MSA */
	int val;			/* return value */
	Parameters rcstruct;		/* configurable parameters */
	long iter;			/* iterations of annealing algorithm */
	long trees_output_total = 0L;	/* number of trees output, overall */
	long trees_output;		/* number of trees output for current rep. */
	long final_length;		/* length of shortest tree(s) found */
	FILE *outtreefp;		/* best trees found overall */
	outtreefp = (FILE *) alloc (sizeof(FILE), "alloc FILE");
	Lvb_bool log_progress;	/* whether or not to log Anneal search */

	#ifdef LVB_MAPREDUCE

	MPI_Init(&argc,&argv);

	MISC misc;

	MPI_Comm_rank(MPI_COMM_WORLD,&misc.rank);
	MPI_Comm_size(MPI_COMM_WORLD,&misc.nprocs);

	MapReduce *mrTreeStack = new MapReduce(MPI_COMM_WORLD);
	mrTreeStack->memsize = 1024;
	mrTreeStack->verbosity = 0;
	mrTreeStack->timer = 0;

	MapReduce *mrBuffer = new MapReduce(MPI_COMM_WORLD);
	mrBuffer->memsize = 1024;
	mrBuffer->verbosity = 0;
	mrBuffer->timer = 0;

	#endif

    /* entitle standard output */
    PrintLVBCopyright();
	PrintLVBInfo();

    /* start timer */
    clock_t Start, End;
    double overall_time_taken;

    Start = clock();
    lvb_initialize();

    getparam(&rcstruct, argc, argv);
    StartTime();

    /* read and alloc space to the data structure */
	MSA = (data *) alloc(sizeof(DataStructure), "alloc data structure");
    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, MSA);

    /* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
	treestack = CreateNewTreestack();
	if(rcstruct.algorithm_selection ==2)
    stack_treevo = CreateNewTreestack();

    matchange(MSA, rcstruct);	/* cut columns */
	#ifdef LVB_MAPREDUCE
    writeinf(rcstruct, MSA, argc, argv, misc.nprocs);
	#else
	writeinf(rcstruct, MSA, argc, argv);
	#endif
    calc_distribution_processors(MSA, rcstruct);

    if (rcstruct.verbose == LVB_TRUE) {
		printf("MinimumTreeLength: %ld\n\n", MinimumTreeLength(MSA));
    }
    rinit(rcstruct.seed);
	log_progress = LVB_TRUE;

    outtreefp = clnopen(rcstruct.file_name_out, "w");
    FILE * treEvo;
	treEvo = (FILE *) alloc(sizeof(FILE), "alloc FILE");
    if(rcstruct.algorithm_selection ==2)
    treEvo = fopen ("treEvo.tre","w");
	iter = 0;
	#ifdef LVB_MAPREDUCE
	final_length = GetSoln(MSA, rcstruct, &iter, log_progress, &misc, mrTreeStack, mrBuffer);
	if (misc.rank == 0) {
		trees_output = PrintTreestack(MSA, &treestack, outtreefp, LVB_FALSE);
	}

	#else
	final_length = GetSoln(MSA, rcstruct, &iter, log_progress);
	trees_output = PrintTreestack(MSA, &treestack, outtreefp, LVB_FALSE);

	#endif

	trees_output_total += trees_output;
    if(rcstruct.algorithm_selection ==2)
		PrintTreestack(MSA, &stack_treevo, treEvo, LVB_FALSE);
    ClearTreestack(&treestack);
	printf("--------------------------------------------------------\n");
	#ifdef LVB_MAPREDUCE
	/* clean the TreeStack and buffer */
	mrTreeStack->map( mrTreeStack, map_clean, NULL );
	mrBuffer->map( mrBuffer, map_clean, NULL );
	/* END clean the TreeStack and buffer */
	#endif

	if(rcstruct.algorithm_selection ==2)
    fclose(treEvo);
	
	clnclose(outtreefp, rcstruct.file_name_out);

    End = clock();

	overall_time_taken = ((double) (End - Start)) /CLOCKS_PER_SEC;

	PrintLogFile(iter, trees_output_total, final_length, overall_time_taken);

	double consistency_index = MinimumTreeLength(MSA);
	double homoplasy_index = 0;

	consistency_index = consistency_index/final_length;
	homoplasy_index = 1 - consistency_index;

	PrintOutput(iter, trees_output_total, final_length, consistency_index, homoplasy_index, overall_time_taken, rcstruct.file_name_out);

	/* "file-local" dynamic heap memory */
    if (rcstruct.algorithm_selection ==2)
    FreeTreestackMemory(MSA, &stack_treevo);
	FreeTreestackMemory(MSA, &treestack);
    rowfree(MSA);
    free(MSA);

    if (cleanup() == LVB_TRUE) val = EXIT_FAILURE;
    else val = EXIT_SUCCESS;

	#ifdef LVB_MAPREDUCE
	FreeTreestackMemory(MSA, &treestack);
	MPI_Barrier(MPI_COMM_WORLD);

	delete mrTreeStack;
	delete mrBuffer;

	MPI_Finalize();
	#endif

    return val;

} /* end main() */
