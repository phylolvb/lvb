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

/* ========== Hash.cpp - hashing functions ========== */

#include "Hash.h"

long CompareHashTreeToHashstack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
{
    long i = 0, new_root = 0;
    static TREESTACK_TREE_BRANCH *copy_2 = NULL;			/* possibly re-rooted tree 2 */
    Lvb_bool b_First = LVB_TRUE;
    string current_site_states;
    unsigned long current_site_states_hash = 0;
    static vector<unsigned long> hashstackvector;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
    treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0){
    	lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

    if (sp->next == 0){
     	current_site_states = MakeHashSet(MSA, copy_2, new_root /* always root zero */);
      hashstackvector.clear();
      current_site_states_hash = HashSiteSet(current_site_states);
    } else{
            for (i = sp->next - 1; i >= 0; i--) {
            if (TopologicalHashComparison(MSA, hashstackvector.at(i), copy_2, b_First, current_site_states, current_site_states_hash) == 0) return 0;
                b_First = LVB_FALSE;
              }
          }
    hashstackvector.push_back(current_site_states_hash);

    /* topology is new so must be pushed */
    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);
    return 1;

} /* end CompareHashTreeToHashstack() */

//if hash != return 1, else return 0
long TopologicalHashComparison(Dataptr MSA, unsigned long stored_hash, const TREESTACK_TREE_BRANCH *const tree_2, Lvb_bool b_First, string current_site_states, unsigned long& current_site_states_hash) {
  if (b_First == LVB_TRUE) {
    current_site_states = MakeHashSet(MSA, tree_2, 0 /* always root zero */);
    current_site_states_hash = HashSiteSet(current_site_states);
  }
  return HashComparison(stored_hash, current_site_states_hash);
}

long HashComparison(unsigned long stored_hash, unsigned long current_site_states_hash) {
    if (stored_hash != current_site_states_hash) return 1;
    return 0;
}

string ConvertSiteSetToString(Dataptr MSA, Objset *oset_1)
{
  ostringstream os;
	for (int i = 0; i < MSA->nsets; i++){
		os << i << "    " << oset_1[i].cnt << "    ";
		for (int x = 0; x < oset_1[i].cnt; x++) 
    os << oset_1[i].set[x] << "   ";
		os << endl;
	}
  
  string sitesetstr(os.str());
  return sitesetstr;
}

unsigned long HashSiteSet(string currentsiteset)
{
  unsigned long str_hash = hash<string>{}(currentsiteset);
  return str_hash;
}

unsigned long Hashcurrent_site_states()
{
  ifstream file;
  file.open("PrintObjectset");

  stringstream strStream;
  strStream << file.rdbuf();
  string str = strStream.str();

  unsigned long str_hash = hash<string>{}(str);

  FILE *printallhash = fopen("PrintAllHashes", "a+");
    fprintf(printallhash, "%lu \n", str_hash);
  fclose(printallhash);

  return str_hash;
}