
#include <iostream>

using namespace std;

#define PROGNAM "LVB "			
#define LVB_VERSION "4.0"	
#define LVB_MAP_REDUCE "MapReduce Multicore version "
#define LVB_RELEASE_DATE "February 2019" 
#define LVB_WIKI "github.com/phylolvb/lvb"

#ifdef NP_Implementation
extern "C" void print_LVB_COPYRIGHT();
extern "C" void print_LVB_INFO();
#endif

void print_LVB_COPYRIGHT();
void print_LVB_INFO();