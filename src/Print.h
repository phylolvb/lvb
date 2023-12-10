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
(c) Copyright 2023 by Joseph Guscott and Daniel Barker.

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

/* ========== Print.h - interface for Print.cpp ========== */

#ifndef LVB_PRINT_H_
#define LVB_PRINT_H_

#include <iostream>

#include "LVB.h"
#include "Clock.h"

#define PROGNAM "LVB "
#define LVB_VERSION "4.2"
#define LVB_RELEASE_DATE "April 2023"
#define LVB_WIKI "github.com/phylolvb/lvb"

void PrintCopyright();
void PrintVersion();
void PrintOutput(long, long, long, double, double, double, char *);
void PrintArguments(Arguments, Dataptr, int, char **);

#endif // LVB_INFO_H_
 