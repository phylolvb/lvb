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
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== main.c ========== */

#include "Admin.h"
#include "Clock.h"
#include "DataOperations.h"
#include "Hash.h"
#include "Print.h"
#include "LogFile.h"
#include "LVB.h"
#include "MemoryOperations.h"
#include "arguments.h"
#include "Solve.c"
#include "Solve.h"
#include "Verbose.h"

int main(int argc, char **argv)
{
	/* entitle standard output */
	PrintCopyright();
	PrintVersion();

	/* start timer */
	clock_t start = clock();
	StartTime();

	SystemChecks();

	Arguments args;		  /* configurable arguments */
	GetArguments(&args, argc, argv);

	/* read and alloc space to the data structure */
	Dataptr MSA = (data *)alloc(sizeof(DataStructure), "alloc data structure");
	phylip_dna_matrin(args.file_name_in, args.n_file_format, MSA);

	/* "file-local" dynamic heap memory: set up best tree stacks, need to be by thread */
	treestack = CreateNewTreestack();
	if (args.algorithm_selection == 2)
		stack_treevo = CreateNewTreestack();

	matchange(MSA, args); /* cut columns */
	SetNumThreads(MSA, args);
	PrintArguments(args, MSA, argc, argv);

	if (args.verbose == true)
	{
		printf("MinimumTreeLength: %ld\n\n", MinimumTreeLength(MSA));
	}
	rinit(args.seed);
	Lvb_bool log_progress = LVB_TRUE;

	FILE *outtreefp;			  /* best trees found overall */
	outtreefp = (FILE *)alloc(sizeof(FILE), "alloc FILE");
	outtreefp = clnopen(args.file_name_out, "w");
	FILE *treEvo;
	treEvo = (FILE *)alloc(sizeof(FILE), "alloc FILE");
	if (args.algorithm_selection == 2)
		treEvo = fopen("treEvo.tre", "w");
	long iter = 0; /* iterations of annealing algorithm */
	int tree_length = GetSoln(MSA, args, &iter, log_progress);
	int number_trees_output = PrintTreestack(MSA, &treestack, outtreefp, LVB_FALSE);

	if (args.algorithm_selection == 2)
		PrintTreestack(MSA, &stack_treevo, treEvo, LVB_FALSE);
	ClearTreestack(&treestack);
	printf("--------------------------------------------------------\n");

	if (args.algorithm_selection == 2)
		fclose(treEvo);

	clnclose(outtreefp, args.file_name_out);

	clock_t end = clock();

	double overall_time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

	PrintLogFile(iter, number_trees_output, tree_length, overall_time_taken);

	double consistency_index = MinimumTreeLength(MSA);
	double homoplasy_index = 0;

	consistency_index = consistency_index / tree_length;
	homoplasy_index = 1 - consistency_index;

	PrintOutput(iter, number_trees_output, tree_length, consistency_index, homoplasy_index, overall_time_taken, args.file_name_out);

	/* "file-local" dynamic heap memory */
	if (args.algorithm_selection == 2)
		FreeTreestackMemory(MSA, &stack_treevo);
	FreeTreestackMemory(MSA, &treestack);
	rowfree(MSA);
	free(MSA);

	if (cleanup() == true)
		return EXIT_FAILURE;

	return EXIT_SUCCESS;

} /* end main() */
