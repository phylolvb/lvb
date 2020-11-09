/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2020 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
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

#include "Hash.h"

#define CLADESEP ","	/* clade separator for trees */

void TopologyHashing(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
{
	long current_topology_hash = 0;
	long hashstack_count = 0;

	FILE *printalltopologies = fopen("PrintAllTopologies", "a+");
	FILE *printcurrenttopologyforhash = fopen("PrintCurrentTopologyForHash", "w");

	CallPrintHashTree(MSA, printalltopologies, BranchArray, root);
	CallPrintHashTree(MSA, printcurrenttopologyforhash, BranchArray, root);

	fclose(printalltopologies);
	fclose(printcurrenttopologyforhash);

	current_topology_hash = HashCurrentTree();

	// printf("Current Topology Hash = %lu \n", current_topology_hash);

	hashstack_count = CountHashesInFile();

	// printf("HashStackCount = %ld \n", hashstack_count);

	// CompareHashStringToStringStack(hashstack_count, current_topology_hash);

}

void CallPrintHashTree (Dataptr MSA, FILE *const stream, const TREESTACK_TREE_BRANCH *const BranchArray, const long root)
	/* print tree in BranchArray (of root root) in bracketed text form to stream stream,
	 * in unrooted form */
	{
		PrintHashTree(MSA, stream, BranchArray, root);
	} /* end lvb_treeprint() */


void PrintHashTree(Dataptr MSA, FILE *const stream, const TREESTACK_TREE_BRANCH *const BranchArray, const long root) {
    /* send tree in BranchArray, of root root, to file pointed to by stream in
	 * unrooted form */
	    long obj;					/* current object */
	    static Lvb_bool doneabsroot = LVB_FALSE;	/* have output root */
	    static Lvb_bool usecomma;			/* output clade sep. */
	    char *tmp_title;				/* temporary string */

	    obj = root;

	    if (doneabsroot == LVB_FALSE)	/* print whole tree */
	    {
			/* start tree */
			tmp_title = (char *) alloc(strlen(MSA->rowtitle[obj]) + 1, "temp. title");
			strcpy(tmp_title, MSA->rowtitle[obj]);

			while(tmp_title[strlen(tmp_title) - 1] == ' '){
				tmp_title[strlen(tmp_title) - 1] = '\0';
			}
			fprintf(stream, "(%s", tmp_title);
			free(tmp_title);    /* VERY LOCAL dynamic heap memory */
			usecomma = LVB_TRUE;
			doneabsroot = LVB_TRUE;

			PrintHashTree(MSA, stream, BranchArray, BranchArray[root].left);
			PrintHashTree(MSA, stream, BranchArray, BranchArray[root].right);
			/* end tree */
			fprintf(stream, ");\n");
			if (ferror(stream))
				crash("file error when writing unrooted tree");

			/* clean up for next call */
			usecomma = LVB_FALSE;
			doneabsroot = LVB_FALSE;
	    }
	    else	/* print remainder of tree */
	    {
			if (usecomma == LVB_TRUE) fprintf(stream, "%s", CLADESEP);
			if (root < MSA->n)	/* leaf */
			{
				tmp_title = (char *) alloc(strlen(MSA->rowtitle[obj]) + 1, "temp. title");
				strcpy(tmp_title, MSA->rowtitle[obj]);
				while(tmp_title[strlen(tmp_title) - 1] == ' '){
					tmp_title[strlen(tmp_title) - 1] = '\0';
				}
				fprintf(stream, "%s", tmp_title);
				free(tmp_title);	/* VERY LOCAL dynamic heap memory */
				usecomma = LVB_TRUE;
			}
			else
			{
				fprintf(stream, "(");
				usecomma = LVB_FALSE;
				PrintHashTree(MSA, stream, BranchArray, BranchArray[root].left);
				PrintHashTree(MSA, stream, BranchArray, BranchArray[root].right);
				fputc(')', stream);
				usecomma = LVB_TRUE;
			}
	    }

}

long HashCurrentTree() /* complete */
{
string line;
ifstream myfile ("PrintCurrentTopologyForHash");

FILE *printhashvalue = fopen("PrintAllHashValues", "a+");

string str;
long str_hash = 0;

if (myfile.is_open())
{
    while ( getline (myfile,line))
    {
        str = line;
        str_hash = hash<string>{}(str);

		fprintf(printhashvalue, "%lu\n", str_hash);
		//printf("Current hash value %lu\n", str_hash);
    }
    myfile.close();
}

else cout << "Unable to open file";

fclose(printhashvalue);

return str_hash;
}

long CountHashesInFile() /* complete */
{
	unsigned int number_of_hashes = 0;
	int ch;

	FILE *hashfile = fopen("PrintAllHashValues", "r");

	while (EOF != (ch=getc(hashfile)))
		if ('\n' == ch)
			++number_of_hashes;
	return number_of_hashes;

	fclose(hashfile);
}

/* for number of lines {
	if compare == True, break
	else add
} */

long CompareHashStringToStringStack(long hashstack_count, long current_topology_hash)
{
	int i; /* loop counter */
	FILE *hashfile = fopen("PrintHashTree", "r");

	/* for hashstack_count iterations compare propsed hash to each line, if match break, else add */

	for (i = 0; i < hashstack_count; i++) {
		// compare current_hash to hash at line n
	}
	
	printf("Current_Hash = %ld \n", current_topology_hash);
	// printf("Hashstack_count_in_comparison = %ld \n", hashstack_count);
	// printf("Hashstack_count_in_myfile = %ld \n", hashfile);

	fclose(hashfile);

	return 0;
}