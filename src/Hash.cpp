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
    std::string current_site_states;
    unsigned long current_site_states_hash = 0;
    static std::vector<unsigned long> hashstackvector;
    
    std::string key;
    std::string index; 
    static std::vector<std::vector<std::string>> hashtable;

    PushToHashTable("100","200", hashtable);

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
    treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0){
    	lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

  /* if treestack is empty, add current config */
    if (sp->next == 0){
     	current_site_states = MakeHashSet(MSA, copy_2, new_root); /* make sset and return sset as string */
      hashstackvector.clear();
      current_site_states_hash = HashSiteSet(current_site_states);

      key = current_site_states;
      index = current_site_states_hash;
    } else{
            for (i = sp->next - 1; i >= 0; i--) {
            if (TopologicalHashComparison(MSA, hashstackvector.at(i), copy_2, b_First, current_site_states, current_site_states_hash) == 0) return 0; /* if current hash matches stored hash, exit */
                b_First = LVB_FALSE;
              }
          }
    /* add new hash to vector */
    hashstackvector.push_back(current_site_states_hash);

    /* topology is new so must be pushed */
    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);
    return 1;

} /* end CompareHashTreeToHashstack() */

 /*if hash != return 1, else return 0 */
long TopologicalHashComparison(Dataptr MSA, unsigned long stored_hash, const TREESTACK_TREE_BRANCH *const tree_2, Lvb_bool b_First, std::string current_site_states, unsigned long& current_site_states_hash) {
  if (b_First == LVB_TRUE) {
    current_site_states = MakeHashSet(MSA, tree_2, 0); /* make sset and return sset as string */
    current_site_states_hash = HashSiteSet(current_site_states); /* hash sset string */
  }
  return HashComparison(stored_hash, current_site_states_hash); 
}

long HashComparison(unsigned long stored_hash, unsigned long current_site_states_hash) {
    if (stored_hash != current_site_states_hash) return 1;
    return 0;
}

unsigned long HashSiteSet(std::string currentsiteset)
{
  unsigned long str_hash = std::hash<std::string>{}(currentsiteset);
  return str_hash;
}

void PushToHashTable(std::string index, std::string key, std::vector<std::vector<std::string>> hashtable) {
    
    FILE * printhashtable = fopen ("src/HashTable","w");   
    
    for (int i = 0; i < 3; i++) {
      std::vector<std::string> index_key;
      for (int j = 0; j < 1; j++) {
        index_key.push_back(index);      // index (hash)
        index_key.push_back(key);      // key(tree)
      }
      hashtable.push_back(index_key);
    }

    for(auto i = hashtable.begin(); i != hashtable.end(); i++)
      fprintf(printhashtable, "%s\t%s\n", index.c_str(), key.c_str());

    /* print hash table */
    for (int i = 0; i < hashtable.size(); i++) {
      for (int j = 0; j < hashtable[i].size(); j++) {
        std::cout << hashtable[i][j];
      }
      std::cout << std::endl;
    }

    fclose(printhashtable);
}
   