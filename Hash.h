#ifndef LVB_HASH_H_
#define LVB_HASH_H_

using namespace std;
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <functional>
#include <unordered_set>
#include <vector>
#include <iterator>
#include <bits/stdc++.h>

#include "LVB.h"

long TopologicalHashComparison(Dataptr restrict, unsigned long, const TREESTACK_TREE_BRANCH *const, Lvb_bool b_first, unsigned long&);
unsigned long HashCurrentSiteStates();
long HashComparison(unsigned long, unsigned long);
long CompareHashTreeToHashstack(Dataptr, TREESTACK *, const TREESTACK_TREE_BRANCH *const, const long, Lvb_bool b_with_sitestate);

#endif