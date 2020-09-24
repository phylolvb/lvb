#ifdef LVB_NP

/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
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

#include "info.h"

void PrintLVBInfo() {
  cout << "==============================================="
          "=================================\n\n";
  cout << "LVB v." << LVB_VERSION << " ";
  #ifdef LVB_MAPREDUCE  
    cout << LVB_MAPREDUCE_VERSION;
  #endif
  cout << "built for Linux 64-bit \n";

  cout << "Released: " << LVB_RELEASE_DATE " by the Barker Lab\n"
  "Developed by: Joseph Guscott and Daniel Barker\n"
  "For help, see the GitHub Wiki page at: " << LVB_WIKI "\n"
  "Please send any questions to joseph.guscott@ed.ac.uk"
  " or daniel.barker@ed.ac.uk\n\n";
  cout << "==============================================="
          "=================================\n\n";
}

void PrintLVBCopyright() {
  cout << "(c) Copyright 2003-2012 by Daniel Barker\n"
  "(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl\n"
  "(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro\n"
  "and Maximilian Strobl\n"
  "(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro,\n"
  "Maximilian Strobl and Chris Wood.\n"
  "(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,\n"
  "Fernando Guntoro, Maximilian Strobl and Chris Wood.\n"
  "(c) Copyright 2019 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,\n"
  "Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, "
  "Martyn Winn and Chris Wood.\n"
  "All rights reserved.\n"
  "\n"
  "Redistribution and use in source and binary forms, with or without\n"
  "modification, are permitted provided that the following conditions\n"
  "are met:\n"
  "\n"
  "1. Redistributions of source code must retain the above copyright\n"
  "notice, this list of conditions and the following disclaimer.\n"
  "\n"
  "2. Redistributions in binary form must reproduce the above\n"
  "copyright notice, this list of conditions and the following\n"
  "disclaimer in the documentation and/or other materials provided\n"
  "with the distribution.\n"
  "\n"
  "3. Neither the name of the copyright holder nor the names of its\n"
  "contributors may be used to endorse or promote products derived\n"
  "from this software without specific prior written permission.\n"
  "\n"
  "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
  "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
  "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS\n"
  "FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
  "COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,\n"
  "INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"
  "INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n"
  "SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)\n"
  "HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,\n"
  "STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)\n"
  "ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF\n"
  "ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n";
}

#elif LVB_PARALLEL_SEARCH



#endif