#include "Hash.h"

long CompareHashTreeToHashstack(Dataptr MSA, TREESTACK *sp, const TREESTACK_TREE_BRANCH *const BranchArray, const long root, Lvb_bool b_with_sitestate)
{
    long i = 0, new_root = 0;
    static TREESTACK_TREE_BRANCH *copy_2 = NULL;			/* possibly re-rooted tree 2 */
    Lvb_bool b_First = LVB_TRUE;
    unsigned long current_hash = 0;
    static vector<unsigned long> hashstackvector;

	/* allocate "local" static heap memory - static - do not free! */
	if (copy_2 == NULL) copy_2 = treealloc(MSA, b_with_sitestate);
    treecopy(MSA, copy_2, BranchArray, b_with_sitestate);
    if (root != 0){
    	lvb_reroot(MSA, copy_2, root, new_root, b_with_sitestate);
    }

    if (sp->next == 0){
     	makesets(MSA, copy_2, new_root /* always root zero */);
        hashstackvector.clear();
        current_hash = HashCurrentSiteStates();
    } else{
            for (i = sp->next - 1; i >= 0; i--) {
            if (TopologicalHashComparison(MSA, hashstackvector.at(i), copy_2, b_First, current_hash) == 0) return 0;
                b_First = LVB_FALSE;
              }
          }
    hashstackvector.push_back(current_hash);

    /* topology is new so must be pushed */
    lvb_assert(root < MSA->n);
    PushCurrentTreeToStack(MSA, sp, BranchArray, root, b_with_sitestate);

    return 1;

} /* end CompareHashTreeToHashstack() */

unsigned long HashCurrentSiteStates()
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

//if hash != return 1, else return 0
long TopologicalHashComparison(Dataptr MSA, unsigned long stored_hash, const TREESTACK_TREE_BRANCH *const tree_2, Lvb_bool b_First, unsigned long& current_hash) {
  if (b_First == LVB_TRUE) {
    makesets(MSA, tree_2, 0 /* always root zero */);
    current_hash = HashCurrentSiteStates();

    // cout << stored_hash << " VS " << current_hash << endl;
  }
  return HashComparison(stored_hash, current_hash);
}

long HashComparison(unsigned long stored_hash, unsigned long current_hash) {
    if (stored_hash != current_hash) return 1;
    return 0;
}
