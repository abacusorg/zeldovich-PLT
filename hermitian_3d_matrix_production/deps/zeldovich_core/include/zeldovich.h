#pragma once

#include <complex>

/* v1.0 -- Initial version

v1.1 -- Fixed critical bug in Box-Muller implementation.
Fixed minor bug in spline_function initialization that was causing
crashes under Linux.

v1.2 -- Added feature to output only one XY slab (at a chosen z)
so as to keep the files smaller for debugging or plotting.

v1.3-- Changed to use ParseHeader to handle input files

v1.4-- Changed output to match input specification for abacus

v1.5-- Changed standard random call to mersene twister from the GSL

v1.6-- Support for "oversampled" simulations (same modes at different PPD) via the
k_cutoff option

v1.7-- Support for PLT eigenmodes and rescaling

v2.0-- Support for oversampling with and without new power.
       N.B. The generation of modes from RNG has changed!
       Use "ZD_Version = 1" to get the old phases
       (but beware version 1 phases depend on ZD_NumBlock).
*/
#define VERSION "zeldovich_v2.0"

using Complx = std::complex<double>;

const int64_t MAX_PPD = 65536;
