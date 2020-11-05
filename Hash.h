#include "LVB.h"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <functional>
#include <unordered_set>

using namespace std;

static void PrintCurrentTree(Dataptr restrict, FILE *const stream, const TREESTACK_TREE_BRANCH *const BranchArray, const long root);
