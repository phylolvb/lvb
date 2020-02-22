/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2019 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
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

/* ========== Mathematical_Wrapper.c - wrappers for standard maths functions ========== */

/*

Provides wrappers for mathematical functions in the standard library.
These wrappers include error handling in the sense that they crash
verbosely if the standard library function signifies an error through
C<errno> or its return value.

**********/

#include "lvb.h"

double exp_wrapper(double x)
/*@globals errno@*/ /*@modifies nothing@*/
{
    double val;	/* return value */

    val = exp(x);
    if (errno == EDOM)
        crash("internal error detected in function exp_wrapper():\n"
         "domain error. x is %g", x);
    else if ((val == HUGE_VAL) || (val == - HUGE_VAL))
        crash("internal error detected in function exp_wrapper():\n"
         "range error. x is %g", x);
    else if (val == 0.0)
        crash("internal error detected in function exp_wrapper():\n"
         "underflow. x is %g", x);
    return val;

}	/* end exp_wrapper() */

double log_wrapper(double x)
/*@globals errno@*/ /*@modifies nothing@*/
{
    double val;	/* return value */

    if (x <= 0.0)
    {
        crash("internal error detected in function log_wrapper():\n"
         "domain error. x is %g\n", x);
    }
    val = log(x);
    if (errno == EDOM)
        crash("internal error detected in function log_wrapper():\n"
         "domain error. x is %g\n", x);
    return val;

}	/* end log_wrapper() */

double pow_wrapper(double x, double y)
/*@globals errno@*/ /*@modifies nothing@*/
{
    double val = -1;	/* return value, initialized to stop warnings */

    if (x < 0)
    {
        if ((ceil(y) != y) || (floor(y) != y))
        {
            crash("internal error detected in function pow_wrapper():\n"
            #ifndef NP_Implementation
             "domain error. x is %g, y is %g, ceil(y) is %g, floor(y) is %g",
		        x, y, ceil(y), floor(y));
            #else
             "domain error. x is %g, y is %g", x, y);
            #endif
        }
    }
    else if (x == 0.0) 
    {
        if (y > 0.0)
            val = 0.0;
        else if (y == 0.0)
            crash("internal error detected in function pow_wrapper():\n"
             "domain error. x is %g, y is %g", x, y);
    }
    else
    {
        val = pow(x, y);
        if (errno == EDOM)
            crash("internal error detected in function pow_wrapper():\n"
             "domain error. x is %g, y is %g", x, y);
        else if ((val == HUGE_VAL) || (val == - HUGE_VAL))
            crash("internal error detected in function pow_wrapper():\n"
             "range error. x is %g, y is %g", x, y);
        else if (val == 0.0)
            crash("internal error detected in function pow_wrapper():\n"
             "underflow. x is %g, y is %g", x, y);
    }
    return val;

}	/* end pow_wrapper() */