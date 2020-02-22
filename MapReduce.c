/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2019 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
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

#ifdef MAP_REDUCE_SINGLE

	#include "lvb.h"

	void map_clean(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr)
	{
	  kv->add(NULL, 0, NULL, 0);
	}

	void reduce_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
	{
	   char *value;
	   MISC *misc = (MISC *) ptr;
	   int check;
	   int ID;

	   uint64_t totalnvalues;
	//   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
	//   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

		int macro_nblocks = 1;
		totalnvalues = nvalues;
	//	cerr << "totalnvalues 1 = " << totalnvalues << endl;
		MapReduce *macro_mr = NULL;
		if (!(multivalue)) {
			macro_mr = (MapReduce *) (valuebytes);
			totalnvalues = macro_mr->multivalue_blocks(macro_nblocks);
		}
	//	cerr << "totalnvalues 2 = " << totalnvalues << endl;

		for (int macro_iblock = 0; macro_iblock < macro_nblocks; macro_iblock++) {
			if (macro_mr)
				(nvalues) = macro_mr->multivalue_block(macro_iblock, &(multivalue),&(valuebytes));

			check = 0;
		   value = multivalue;
		//   cerr << "value 2 0 = " << *value << " " << nvalues << endl;
		   for (int i=0; i<nvalues; i++) {
				ID = *(int *) value;
				if(ID == 0) {
					check = 1;
					break;
				}
				value += valuebytes[i];
				ID = *(int *) value;
		//		cerr << "value 2 = " << i << " " << ID << endl;
		   }

		   if (check == 1) {
				value = multivalue;
		//		cerr << "value 3 0 = " << *value << endl;
				for (int i=0; i<nvalues; i++) {
					ID = *(int *) value;
					misc->count[ID]++;
					value += valuebytes[i];
					ID = *(int *) value;
		//			cerr << "value 3 = " << i << " " << ID << " " << valuebytes[i] << endl;
				}
		   }

		}

	// if(misc->rank ==0) {
	//   long *set;
	//   set = (long *) key;
	//   int n = (int) (keybytes / sizeof(long));
	//   for (int i=0; i<n; i++) {
	//	cerr << "\t" << set[i];
	//   }
	//   cerr << endl << " -------------- " << endl;
	// }

	}

	void reduce_filter(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr)
	{
	   char *value_i, *value_j;
	   int ID_i, ID_j;
	   int check;
	   uint64_t nvalues_total;

		if(nvalues > 1) {

		   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
		   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

			 value_i = multivalue;
			 for (int i=0; i<(nvalues-1); i++) {
			check = 0;
			ID_i = *(int *) value_i;
				value_j = value_i;
				for(int j=(i+1); j<nvalues; j++) {
				ID_j = *(int *) value_j;
				if(ID_i == ID_j) {
					check = 1;
					break;
				}
						value_j += valuebytes[j];
				}
			if(check == 0) kv->add(key, keybytes, (char *) &ID_i, sizeof(int));
				value_i += valuebytes[i];
			 }

		   END_BLOCK_LOOP

		} else if (nvalues == 1) {

			 CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
			 BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

			 value_i = multivalue;
			 kv->add(key, keybytes, value_i, valuebytes[0]);

			 END_BLOCK_LOOP
		}

	}
#endif
