/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl
and Chris Wood
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Chang Sik Kim,
Maximilian Strobl and Martyn Winn
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

// macros to simplify use of multi-page KMVs
// Karen Devine, March 2010

#ifdef MAP_REDUCE_SINGLE

	#ifndef _BLOCKMACROS_HPP
	#define _BLOCKMACROS_HPP

	// macros to loop over blocks when reduce multivalues may span more than 1 block
	// use CHECK_FOR_BLOCKS initially to get # of blocks in the multivalue
	// enclose code for each block between BEGIN_BLOCK_LOOP and END_BLOCK_LOOP
	// NOTE: DO NOT put a semicolon afer these macros

	#define CHECK_FOR_BLOCKS(multivalue, valuebytes, nvalues, totalnvalues)  \
	  int macro_nblocks = 1; \
	  totalnvalues = nvalues; \
	  MapReduce *macro_mr = NULL; \
	  if (!(multivalue)) { \
		macro_mr = (MapReduce *) (valuebytes); \
		totalnvalues = macro_mr->multivalue_blocks(macro_nblocks); \
	  }

	#define BEGIN_BLOCK_LOOP(multivalue, valuebytes, nvalues)  \
	  for (int macro_iblock = 0; macro_iblock < macro_nblocks; macro_iblock++) { \
		if (macro_mr)  \
		  (nvalues) = macro_mr->multivalue_block(macro_iblock, \
												 &(multivalue),&(valuebytes));

	#define END_BLOCK_LOOP }

	#endif

#endif