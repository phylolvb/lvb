#include "LVB.h"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <functional>
#include <unordered_set>

using namespace std;

static void PrintHashTree(Dataptr restrict, FILE *const stream, const TREESTACK_TREE_BRANCH *const BranchArray, const long root);
long CompareHashToHashStack(unsigned long, long);
long CountHashesInFile();
long ConvertHashStackToArray(long *, long);