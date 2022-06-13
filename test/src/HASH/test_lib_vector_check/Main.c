/* LVB

(c) Copyright 2003-2012 by Daniel Barker.
(c) Copyright 2013, 2014 by Daniel Barker and
Maximilian Strobl.
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2022 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Chang Sik Kim, Fernando Guntoro, Maximilian Strobl, Chris Wood
and Martyn Winn.

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

#include "src/LVB.h"

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

/* positive test for cistrcmp() */

int main(void)
{

    const unsigned long test_hash_1 = 179349889;
    const unsigned long test_hash_2 = 782364598;
    bool hash_found_1 = false;
    bool hash_found_2 = false;

    static vector<unsigned long> testvector = {357435487, 158374413, 873543581, 179349889, 357431354, 574318349};

    for(int i = 0; i < testvector.size(); i++) {
        if(test_hash_1 == testvector.at(i)) hash_found_1 = true; /* test_hash_1 located at position 3*/
    }

    for(int i = 0; i < testvector.size(); i++) {
        if(test_hash_2 == testvector.at(i)) hash_found_2 = true; /* test_hash_2 not located in test_vector */
    }

    if(hash_found_1 == true && hash_found_2 == false) {
        printf("test passed\n");
    } else {
        printf("test failed\n");
    }

    return 0;
}
