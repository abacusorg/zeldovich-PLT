# Random Number Generation for zeldovich-PLT
Author: Lehman Garrison

zeldovich-PLT requires fast, parallel random number generation (RNG)
for realizing the Gaussian random field upon which the initial conditions
are based.  This document describes the meta-RNG (and meta-meta-RNG)
strategy used by this code to support parallelism and oversampling.

# Meta RNG

Generating initial conditions requires generating a Fourier cube
of complex numbers, with amplitudes drawn from Gaussians with
mean and width determined by the power spectrum and phases drawn
from a uniform random distribution.  So, we need to fill a ppd^3
cube with (complex) random numbers (ppd = "particle per dimension").
We want these "random" numbers to be reproducible, so we use a deterministic
pseudo-RNG like the Mersenne Twister instead of true randomness.

But ppd^3 is potentially large (100s of billions), and we also
need to do some processing of each complex number in addition 
to the RNG call(s).  So we naturally seek to parallelize.  But
most RNGs are sequential, so we can't have multiple threads all
querying the same RNG but filling different spots in the Fourier
cube.  This would not be reproducible (or at least slow if
synchronized appropriately).

It's natural to consider giving each thread an RNG, but then
the results would depend on the number of threads.  In zeldovich-PLT,
we parallelize over y-planes in Fourier space, so a natural
solution is to give each of the ppd planes an RNG.  That way,
it doesn't matter which thread is doing the work.

But we also want to support oversampling, in which we can generate
phase-matched ICs with increasing particle number.  In Fourier
space, this means that the boundaries of the cube are growing,
so the challenge becomes generating the same random numbers
for the part of the cube that is unchanged, while generating
new random numbers for the new part that has appeared.

If one fills out to the new k_max = k_Nyquist for a Fourier
pencil (i.e. all kx for a given ky,kz), the the RNG-per-plane
scheme would mean that the RNG would be in a different spot
when starting the next pencil in the plane.  So that does not
support oversampling.

So then we arrive at one RNG per pencil, which is essentially
what we do in zeldovich-PLT.  But we actually need two per pencil,
because we need to fill all interior modes (i.e. modes that don't
change when oversampling) before moving onto exterior ones.
And the FFTW packing of FFT modes for a pencil is ordered:
```
[ 0, 1, ..., N/2-1, -N/2, ..., -1]*(2*pi/L)
```
so we need to fill inwards from the left and right edges (i.e.
lowest wavenumbers) before proceeding to the center.  We could
alternative left and right from a single RNG, but that would
probably incur some code duplication in the inner loop of the
cube-filling.  It's also possible an alternating scheme would
be slower due to requiring (pre)fetching out of multiple locations
in memory and lower reuse of cache lines.

Regardless of 1 or 2 RNG per pencil, ppd^2 pencils is a
non-trivial number of RNG to initialize.  Using GSL RNG,
each initialization requires a memory allocation and a
bit of math to initialize the internal state from the seed.
Furthermore, we want to warm up each RNG with a few calls
(for reasons that will be explained below), so even if we were
to write our own "mass initializer", we would want to
parallelize the process.

So the challenge becomes generating 2\*ppd^2 seeds in a
parallel fashion from a base/master seed (this is the seed
specified in the simulation parameter file).  We don't make each
seed a simple increment over its predecessor, because we sometimes
run suites of simulations whose seeds we simply increment by 1.
So we can consider seeding each RNG with the output of a
"meta RNG" seeded with the base seed.  But how do we parallelize
querying the meta RNG?  We can parallelize over rows, much as
we parallelized over planes in the 3D problem, but now we're back
to our initial problem of parallelizing the RNG!  And we can solve
it the same way: with a "meta-meta RNG" that gives the seeds
for the ppd meta RNGs that give the seeds for the 2\*ppd^2 RNGs.
The meta-meta RNG is seeded with the base seed.

The meta RNG needs to be seeded in the same "outside in" order
as the primary RNG to synchronize Fourier cubes of different sizes.
But here we just operate in inward-marching pairs and duplicate the
code in the inner loop, since it's quite minimal.  For the primary
RNG, we instead changed the loop logic to switch direction halfway
through.


# Warm up

An RNG may take a few rounds to each a high-entropy state if the
initial seed had low entropy (e.g. 1, 2, 3).  Since we commonly
choose low-entropy seeds for convenience, it makes sense to do
a few rounds of warm-up.  This is mostly true for the meta-meta RNG,
since the meta RNG seeds will have high entropy if we initialize
the meta-meta RNG properly.  But we do a few rounds of warm up
at all three RNG levels since it's cheap.  Many GSL RNGs seem
to do a tiny bit of warm up internally, but we know that we can
afford to do more.

Our meta-seeding is further motivation to warm up.  Some RNGs
yield the seed as the first value or the value before the first.
This would correlate the levels of the meta RNG hierarchy, or at
least pass low entropy seeds to the bottom.

See, e.g., http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf

# RNG algorithm and memory

2\*ppd^2 RNGs is not a trivial number.  For ppd=7000, that's 100 million.
The memory requirements on the RNG state alone may be large if we aren't
careful.  For example, GSL's default RNG (mt19937) uses 5K of state per
RNG, a huge amount!  This would result in 500 GB of memory usage from
the RNGs alone.  Fortunately, there is at least one simulation-quality,
low footprint RNG in GSL: the taus2 "Tausworthe" generator, which uses
24 bytes of state.  This is 2.4 GB for ppd 7000, and is thus unlikely to
be a limiting factor.

Having less state, the taus2 RNG necessarily has a lower period than
the mersenne twister.  But it's still huge: 2^88, or 10^26, far more
than we'll ever use.  taus2 also passes just about all of the diehard
randomness purity tests, performing only a shade worse than mt19937.
It's more than sufficient for generating a Gaussian random field.

We use the GSL library because it's portable and well-tested.  But it
means that we're tied to their initialization routines, which involve
a separate malloc for each RNG.  So 100 million multi-threaded allocations
of 24 bytes each cries out for tcmalloc or a similar thread-aware allocator.
So we now link against tcmalloc when it's available; this is many times
faster than using the built-in allocator.

It's probably possible to write a bulk RNG allocator that still hooks
into GSL, but it only matters for the largest allocation sizes, where
we're already using tcmalloc for Abacus.

Finally, it's worth noting that 

RNG Purity vs complexity
