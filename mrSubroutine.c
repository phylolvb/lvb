
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
 
   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   check = 0;
   value = multivalue;
   for (int i=0; i<nvalues; i++) {
	ID = *(int *) value;
	if(ID == 0) {
		check = 1;
		break;
	}
	value += valuebytes[i];
   }

   if (check == 1) {
    value = multivalue;
    for (int i=0; i<nvalues; i++) {
        ID = *(int *) value;
	misc->count[ID]++; 
	value += valuebytes[i];
    }
   }

   END_BLOCK_LOOP

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

