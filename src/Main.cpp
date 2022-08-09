/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl.
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.
(c) Copyright 2022 by Joseph Guscott and Daniel Barker.

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

/* ========== Main.cpp ========== */

#ifdef LVB_MPI
	#include <mpi.h>
	#include <iostream>
	#include <vector>
	#include <algorithm>

using namespace std;
#endif

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

	#ifdef LVB_MPI
		/* Parallel properties */

		int clusterSize;
		int rank;
		int rootProcess = 0;
		int errorMPI;
		char fileNameMPI[3000];

		MPI_Status status;

		errorMPI = MPI_Init(&argc, &argv);

		errorMPI += MPI_Comm_size(MPI_COMM_WORLD, & clusterSize);
		errorMPI += MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		/* Quick checks */

		if(errorMPI != 0) {
			printf("MPI failed to initialise\n");
			exit(1);
		}

		if(clusterSize <= 1) {
			printf("Number of MPI processes must be greater than one\n");
			MPI_Finalize();
			exit(1);
		}
	#endif



    /* entitle standard output */

	#ifdef LVB_MPI
		if(rank == 0) {
    		PrintLVBCopyright();
			PrintLVBInfo();
		}
	#else
		PrintLVBCopyright();
		PrintLVBInfo();
	#endif

    /* start timer */
    clock_t Start, End;
    double overall_time_taken;

    Start = clock();
    lvb_initialize();

    getparam(&rcstruct, argc, argv);

	#ifdef LVB_MPI
		if(rank == 0)
    		StartTime();
	#else
		StartTime();
	#endif

    /* read and alloc space to the data structure */
	MSA = (data *) alloc(sizeof(DataStructure), "alloc data structure");
    phylip_dna_matrin(rcstruct.file_name_in, rcstruct.n_file_format, MSA);

    /* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
	treestack = CreateNewTreestack();
	if(rcstruct.algorithm_selection ==2)
    stack_treevo = CreateNewTreestack();

    matchange(MSA, rcstruct);	/* cut columns */

	#ifdef LVB_MPI
		if(rank == 0)
			writeinf(rcstruct, MSA, argc, argv, clusterSize);
	#else
		writeinf(rcstruct, MSA, argc, argv);
	#endif

    calc_distribution_processors(MSA, rcstruct);

    if (rcstruct.verbose == LVB_TRUE) {
		printf("MinimumTreeLength: %ld\n\n", MinimumTreeLength(MSA));
    }

	#ifdef LVB_MPI
		int seedMPI = rcstruct.seed;
		vector<int> MPIVect;

		if(rank == 0) {
			for(int i = 0; i < clusterSize; i++)
				MPIVect.push_back(seedMPI + i);
		}

		/* for(auto i: MPIVect)
			cout << i << ' ';
		cout << endl; */

		MPI_Scatter(MPIVect.data(), 1, MPI_INT, &seedMPI, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	rinit(seedMPI);
	#else
		rinit(rcstruct.seed);
	#endif
    
	log_progress = LVB_TRUE;

    outtreefp = clnopen(rcstruct.file_name_out, "w");
    FILE * treEvo;
	treEvo = (FILE *) alloc(sizeof(FILE), "alloc FILE");
    if(rcstruct.algorithm_selection ==2)
    treEvo = fopen ("treEvo.tre","w");
	iter = 0;

	#ifdef LVB_MPI
		final_length = GetSoln(MSA, rcstruct, &iter, log_progress, rank);
	#else
		final_length = GetSoln(MSA, rcstruct, &iter, log_progress);
	#endif

	trees_output = PrintTreestack(MSA, &treestack, outtreefp, LVB_FALSE);

	#ifdef LVB_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		for(int i = 0; i < clusterSize; i++) {
			if(rank == i) {
				sprintf(fileNameMPI, "%s_%d", rcstruct.file_name_out, i);
				outtreefp = clnopen(fileNameMPI, "w");
				trees_output = PrintMPITreestack(MSA, &treestack, outtreefp, rank, LVB_FALSE);
			}
		}

		MPI_Gather(&final_length, 1, MPI_INT, MPIVect.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

		/* for(auto i: MPIVect)
			cout << i << ' ';
		cout << endl; */

		MPI_Barrier(MPI_COMM_WORLD);

		int position = 0;

		 vector<int>::iterator minValue = min_element(MPIVect.begin(), MPIVect.end());
		position = distance(MPIVect.begin(), minValue);
		/* cout << "Position, rank: " << position << ", " << rank << endl; */

		MPI_Barrier(MPI_COMM_WORLD);
	#endif

	trees_output_total += trees_output;
    if(rcstruct.algorithm_selection ==2)
		PrintTreestack(MSA, &stack_treevo, treEvo, LVB_FALSE);
    ClearTreestack(&treestack);

	#ifdef LVB_MPI
		if(rank == 0) {
			printf("--------------------------------------------------------\n");
			printf("================================================================================\n");
			printf("\nMPI Search Results: \n\n");
			printf("----------------------------------------------------------------------\n");
			printf (" File:        Seed:          Topologies:     Score:         Runtime-:\n");
			printf("----------------------------------------------------------------------\n");
		}
	#else
		printf("--------------------------------------------------------\n");
	#endif

	if(rcstruct.algorithm_selection ==2)
    fclose(treEvo);
	
	#ifdef LVB_MPI
		clnclose(outtreefp, fileNameMPI);
	#else
		clnclose(outtreefp, rcstruct.file_name_out);
	#endif

    End = clock();

	overall_time_taken = ((double) (End - Start)) /CLOCKS_PER_SEC;

	#ifdef LVB_MPI
		PrintMPILogFile(iter, trees_output_total, final_length, overall_time_taken, seedMPI);
	#else
		PrintLogFile(iter, trees_output_total, final_length, overall_time_taken, rcstruct.seed);
	#endif

	double consistency_index = MinimumTreeLength(MSA);
	double homoplasy_index = 0;

	consistency_index = consistency_index/final_length;
	homoplasy_index = 1 - consistency_index;

	#ifdef LVB_MPI
		/* if(rank == 0) */
		PrintMPIOutput(iter, trees_output_total, final_length, consistency_index, homoplasy_index, overall_time_taken, fileNameMPI, seedMPI, rank);
	#else
		PrintOutput(iter, trees_output_total, final_length, consistency_index, homoplasy_index, overall_time_taken, rcstruct.file_name_out);
	#endif

	/* "file-local" dynamic heap memory */
    if (rcstruct.algorithm_selection ==2)
    FreeTreestackMemory(MSA, &stack_treevo);
	FreeTreestackMemory(MSA, &treestack);
    rowfree(MSA);
    free(MSA);

    if (cleanup() == LVB_TRUE) val = EXIT_FAILURE;
    else val = EXIT_SUCCESS;

	#ifdef LVB_MPI
		MPI_Finalize();
		if(rank == 0) {
			printf("----------------------------------------------------------------------");
			printf("\n================================================================================\n");
		}
	#endif

    return val;

} /* end main() */
