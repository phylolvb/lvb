#include "LVB.h"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <functional>
#include <unordered_set>
#include <vector>
#include <iterator>
#include <bits/stdc++.h>

using namespace std;

static void PrintHashTree(Dataptr restrict, FILE *const stream, const TREESTACK_TREE_BRANCH *const BranchArray, const long root);
long CompareHashToHashStack(unsigned long, long);
long CountHashesInFile();