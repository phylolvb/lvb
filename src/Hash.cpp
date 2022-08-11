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

/* ========== Hash.cpp - hashing functions ========== */

#include "Hash.h"

long CompareHashTreeToHashstack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray, const long root, Lvb_bool b_with_sitestate, Parameters rcstruct) {

  if(rcstruct.searchSelection == 0) {
    /* LINEAR SEARCH */
    long i = 0, new_root = 0;
    static TREESTACK_TREE_NODES *copy_2 = NULL;
    Lvb_bool b_First = LVB_TRUE;
    std::string current_site_states;
    unsigned long long current_site_states_hash = 0;
    static std::vector<unsigned long long> hashstackvector;

    /* allocate "local" static heap memory - static - do not free! */
    if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
      treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0) {
      lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

    /* if treestack is empty, add current config */
    if (sp->next == 0) {
      current_site_states = MakeHashSet(MSA, copy_2, new_root);
      hashstackvector.clear();
      current_site_states_hash = HashSiteSet(current_site_states);
    } else {    
      for (i = sp->next - 1; i >= 0; i--) {
        if (TopologicalHashComparison(MSA, hashstackvector.at(i), copy_2, b_First, current_site_states, current_site_states_hash) == 0) {
          return 0; /* if current hash matches stored hash, exit */
        }
        b_First = LVB_FALSE;
      }
    
    }
    hashstackvector.push_back(current_site_states_hash);

    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);

    return 1;
  }

  if(rcstruct.searchSelection == 1) {
    /* BINARY SEARCH */
    long i = 0, new_root = 0;
    static TREESTACK_TREE_NODES *copy_2 = NULL;
    Lvb_bool b_First = LVB_TRUE;
    std::string current_site_states;
    unsigned long long current_site_states_hash = 0;
    static std::vector<unsigned long long> hashstackvector;

    /* allocate "local" static heap memory - static - do not free! */
    if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
      treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0) {
      lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

    /* if treestack is empty, add current config */
    if (sp->next == 0) {
      current_site_states = MakeHashSet(MSA, copy_2, new_root);
      hashstackvector.clear();
      current_site_states_hash = HashSiteSet(current_site_states);
    } else {    
      current_site_states = MakeHashSet(MSA, copy_2, 0);
      current_site_states_hash = HashSiteSet(current_site_states);

      if(std::binary_search(hashstackvector.begin(), hashstackvector.end(), current_site_states_hash)){
        return 0;
      }
    }
    hashstackvector.push_back(current_site_states_hash);
    std::sort(hashstackvector.begin(), hashstackvector.end());

    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);

    return 1;
  }

  if(rcstruct.searchSelection == 2) {
      /* SET SEARCH */
    long i = 0, new_root = 0;
    static TREESTACK_TREE_NODES *copy_2 = NULL;
    Lvb_bool b_First = LVB_TRUE;
    std::string current_site_states;
    unsigned long long current_site_states_hash = 0;
    static std::vector<unsigned long long> hashstackvector;
    static std::unordered_set <unsigned long long> hashSet;
    unsigned long long HashKey = 0;

    /* allocate "local" static heap memory - static - do not free! */
    if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
      treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0) {
      lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

    /* if treestack is empty, add current config */
    if (sp->next == 0) {
      current_site_states = MakeHashSet(MSA, copy_2, new_root);
      hashSet.clear();
      HashKey = HashSiteSet(current_site_states);
    } else {
      current_site_states = MakeHashSet(MSA, copy_2, 0);
      HashKey = HashSiteSet(current_site_states);
    
      if(hashSet.find(HashKey) != hashSet.end()) 
        return 0;
    }    
  
    hashSet.insert(HashKey);
    
    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);

    return 1;
  }
} /* end CompareHashTreeToHashstack() */

long TopologicalHashComparison(Dataptr MSA, unsigned long long stored_hash, const TREESTACK_TREE_NODES *const tree_2, Lvb_bool b_First,
                                std::string current_site_states, unsigned long long& current_site_states_hash) {
  if (b_First == LVB_TRUE) {
    current_site_states = MakeHashSet(MSA, tree_2, 0);
    current_site_states_hash = HashSiteSet(current_site_states);
  }
  return HashComparison(stored_hash, current_site_states_hash);
} /* end TopologicalHashComparison() */

long HashComparison(unsigned long long stored_hash, unsigned long long current_site_states_hash) {
  if (stored_hash == current_site_states_hash) return 0;

  return 1;
} /* end HashComparison() */

unsigned long long HashSiteSet(std::string currentsiteset) {
  unsigned long long str_hash = std::hash<std::string>{}(currentsiteset);
  return str_hash;
} /* end HashSiteSet() */

/*
long CollisionResolution(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const copy_2, Lvb_bool b_First) {
  b_First = LVB_FALSE;

  for (long i = sp->next - 1; i >= 0; i--) {
    if (TopologyComparison(MSA, sp->stack[i].p_sitestate, copy_2, b_First) == 0) return 0;
  }
} end CollisionResolution() */

int linearHashSearch(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray, const long root, Lvb_bool b_with_sitestate) {

  
}

int binaryHashSearch(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray, const long root, Lvb_bool b_with_sitestate) {

}

int setHashSearch(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_NODES *const BranchArray, const long root, Lvb_bool b_with_sitestate) {

}