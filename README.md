# zeldovich-PLT

## Authors
Daniel Eisenstein & Lehman Garrison  
lgarrison@cfa.harvard.edu  
https://github.com/lgarrison/zeldovich-PLT


## Overview
This code generates Zel'dovich approximation (ZA) initial conditions (i.e. first-order Lagrangian perturbation theory) for cosmological N-body simulations, optionally applying particle linear theory (PLT)
corrections.  This code does not provide second-order ICs (2LPT), but one can use these ICs with the
config-space 2LPT detailed in [Garrison et al. (2016)](https://arxiv.org/abs/1605.02333).

If you do not intend to use the config-space 2LPT, then it's better to use a Fourier-space 2LPT
code (e.g. [2LPTic](http://cosmo.nyu.edu/roman/2LPT/)) than to rely on ZA, even with PLT corrections.

This code supports two types of PLT corrections: (1) PLT eigenmodes, and (2) rescaling.  Both corrections are most important on small scales (approaching, e.g., k<sub>Nyquist</sub>/4).

1. PLT eigenmodes (`ZD_qPLT`) initializes the simulation in the correct eigenmodes for a particle lattice.  This eliminates transients that arise due to the common assumption that the particle system obeys the continuum modes.

2. Rescaling (`ZD_qPLT_rescale` and `ZD_PLT_target_z`) increases the amplitude of initial power on small scales to exactly cancel out future (unavoidable) under-growth.  This requires (1).

This code does not presently support glass initial conditions, only particle lattices.

The code uses double precision internally, but you can set output format to single or double precision (see the `ICFormat` option).

The code can run small problems in memory, but it can also operate out-of-core to support large problems. Set the `-DDISK` option in the Makefile to turn this on.

## Usage
Build with `make`, and run with `./zeldovich <param_file>`.  An example parameter file (`example.par`) is provided, and all of the options are listed in `parameter.cpp`.  See the "Parameter file options" section below for detailed descriptions of the options.

### Dependencies
Zeldovich-PLT needs FFTW 3 and GSL, and the ParseHeader library needs flex and Bison >= 3.0.  The code has been tested with g++, but it should work with the Intel compilers as well.

### Convergence testing
This code supports testing N-body simulation convergence at linear order by increasing particle density ("oversampling")
for a fixed set of modes.  One should keep fixed the starting redshift, volume, softening length, RNG seed, etc. when doing this kind of convergence testing.

To generate oversampled initial conditions (for example, 128<sup>3</sup> initial conditions that sample the same modes as 64<sup>3</sup> initial conditions), invoke the code twice: once with NP = 64<sup>3</sup> to generate the fiducial simulation, then again with NP = 128<sup>3</sup> and `ZD_k_cutoff = 2` to generate the oversampled.  This truncates the modes in the 128<sup>3</sup> sim at k<sub>Nyquist</sub>/2, and ensures the RNG is appropriately synchronized.  **Important**: do not change `ZD_NumBlock` between invocations!

## Citation
If you use this code, please cite [Garrison et al. (2016)](https://arxiv.org/abs/1605.02333).

## Technical details
We're doing four big `[z][y][x]` transform.  But we don't store the
data that way.  Instead, we block the z & y directions.  So for a 
`PPD^3` problem, we have `NB` blocks of `P=PPD/NB` two-d information.  Each
block contains the full `PPD` x-dimension for all four problems.
We require that `NB` divides `PPD` evenly, just for sanity.

This is done so that we first operate on a y row of blocks, so we
can do the x-z FFT.  Then we operate on a z column of blocks, so that 
we can do the y FFT and output ordered in z.

We also store the four transforms packed into two complex transforms.
We store these two arrays interleaved at the block level.

We order the full array as:

```
double complex data[zblock=0..NB-1][yblock=0..NB-1][arr=0..1][zresidual=0..P-1][yresidual=0..P-1][x=0..PPD-1]
```

where
```
z = zblock*P+zresidual
zblock = floor(z/P)
zresidual = z-P*floor(z/P)
```
and the same for y.

Each block is `4*P*P*PPD = 4*PPD^3/NB^2` numbers.  For 4096 and 64, that's
a total array of 256 Gdoubles = 2 TB.  Each slab is 32 GB and each 
block is 512 MB.  Each skewer of X's is 4096 complex doubles, which is 
64 KB, so memory movements are efficient.

The concept is that `2/NB` of the problem has to fit into memory at one time.
The coefficient is 2 because we have to construct the conjugate wavenumber
simultaneously so as to enforce the correct (anti-)Hermitian structure.

One can refer to the array as `[z][y][x]`, but of course the intent is 
that one operates the outer loop by block, so that one can load a big
chunk from disk.

### FFT packing and mode generation
We're going to load the 4 real FFTs into 2 complex FFTs.  This will
allow us to do the iFFTs as a simple cubic approach.  We'll have
```
Re A = density
Im A = q_x
Re B = q_y
Im B = q_z
```
So `Ahat = D + iF, and Bhat = G + iH`, where
```
D(k) = D(-k)^* = delta(k)
F(k) = F(-k)^* = delta(k)*i k_x/k^2
G(k) = G(-k)^* = delta(k)*i k_y/k^2
H(k) = H(-k)^* = delta(k)*i k_z/k^2
```

In detail, we treat `k_x/k^2` as `j/n^2`, where `n^2 = (j^2+l^2+m^2)*(2*pi/L)`
for integral triples `(j,l,m)`.  We must wrap `0..N-1` to `-N/2+1..N/2` for this
computation.  In this language, `-k` means `(j,l,m) -> (N-j,N-l,N-m)`.

`delta(k)` is a complex Gaussian with variance in `Re(delta)` and `Im(delta)`
of `P_theory/Volume_box`.  This is the correct normalization if the iFFT
is in the FFTW convention (no prefactor of `1/N`).

We will force `delta(k)` to be 0 at `k=0` and also if `j=N/2`, `l=N/2`, or `m=N/2`.
The latter avoids bookkeeping about Nyquist aliasing and can't matter 
physically.

In regards to Nyquist aliasing, if we set `k` and `-k` each time we
generate a `delta`, it doesn't matter if we do an element more than
once.  It just overwrites both elements with a different random
number.  That avoids some Nyquist bookkeeping.

We must load `k` and `-k` at the same time (or play horrid tricks to 
reset and resyncronize the random number generator).  There is no
avoiding having two slabs in memory for this.  Given that our first
sweep is in Y slabs, it is probably easiest to split the half-space
on `y=0`.

## PLT eigenmodes
The PLT eigenmode features of this code are developed and tested in
[Garrison et al. (2016)](https://arxiv.org/abs/1605.02333) based on the work of Marcos, et al. (2006).

Nominally, this code only produces ZA displacements, from which velocities can 
be computed in config space later.  However, using the PLT eigenmodes requires computing the
velocities in Fourier space, so we have to do that here.  Thus, we require four
complex FFTs instead of two.

We provide a precompted set of 128<sup>3</sup> numerical eigenmodes with this code.
The code does linear interpolation if a finer FFT mesh is being used.

## Parameter file options
`ZD_Seed`: *integer*  
The random number seed.  Both this variable and `ZD_NumBlock` affect the output phases.

`ZD_NumBlock`: *integer*  
This is the number of blocks to break the FFT
into, per linear dimension.  This must be an even number; moreover,
it must divide `PPD = NP^(1/3)` evenly.  The default is 2, but you
may need a higher number.

This is a key tuning parameter for the code.  The full problem
requires `32*NP` bytes (or `64*NP` if `ZD_qPLT` is being used), which may exceed the amount of RAM.
The zeldovich code holds `2/NumBlock` of the full volume in memory,
by splitting the problem in 2 dimensions into `NumBlock^2` parts.
Each block therefore is `32*NP/NumBlock^2` bytes.  It
is important that these blocks be larger than the latency of the
disk, so sizes of order 100 MB are useful.  We are holding
`2*NumBlock` such blocks in memory.

Hence, for a computer with `M` bytes of available memory and a
problem of `NP` particles, we need `NumBlock > 64*NP/M`
(preferably to be the next larger even number that divides evenly
into `NP^(1/3)`) and we prefer that `32*NP/NumBlock^2` is 
larger than the latency.

For example, for a `4096^3` simulation, `32*NP` is 2 TB.  If we use
`NumBlock` of 128, then we will need 32 GB of RAM and each block saved
to disk will be 128 MB.  For a `2048^3` simulation, `32*NP` is 256 GB
and `NumBlock` of 16 will require 32 GB of RAM and each block will
be 1024 MB.  For a `8192^3` simulation, `32*NP` is 16 TB and `NumBlock`
of `256` will require 128 GB of RAM and a block size of 256 MB.

If `ZD_k_cutoff != 1`, then the actual `ZD_NumBlock` will be `ZD_NumBlock*ZD_k_cutoff`.
See `ZD_k_cutoff` for details.

Both this variable and `ZD_Seed` affect the output phases.

`ZD_Pk_filename`: *string*  
The file name of the input power spectrum.
This can be a CAMB power spectrum.

`ZD_Pk_scale`: *double*  
This is the quantity by which to multiply the 
wavenumbers in the input file to place them in the units that will be used
in the zeldovich code, in which the fundamental wavenumber is `2*pi` divided
by `BoxSize`.  Default value is `1.0`.

As a common example, one might need to convert between `Mpc^-1` and
`h Mpc^-1` units.  The zeldovich code does not use the `hMpc` value,
so it doesn't know what the units of `BoxSize` are.  If `BoxSize` is in 
`h^-1 Mpc` units, so that `hMpc=1`, and if the `P(k)` file had `k`
in `h Mpc^-1` units, then all is well: use the value of `1.0`.

However, if `BoxSize` were in Mpc units and the input power were in 
`h Mpc^-1` units, then we want to convert the wavenumbers to 
`Mpc^-1` units.  That means multiplying by `h`, so we should use
`ZD_Pk_scale = h`.

`ZD_Pk_norm`: *double*  
The scale at which to normalize the input
`P(k)`, in the same units as `BoxSize`.  For example, if `BoxSize` is
given in `h^-1 Mpc`, then one might choose 8 to select `sigma_8`.
If this value is 0, then the power spectrum will not be renormalized
(but `ZD_Pk_scale` will be applied to the wavevectors, so beware that
the power isn't thrown off, as it does have units of volume).  The
default is 0, but we recommend controlling the normalization.

`ZD_Pk_sigma`: *double*  
The amplitude to use to normalize the
fluctuations of the density field, with a tophat of radius `ZD_Pk_norm`.
This must be scaled to the initial redshift by the growth function; 
the zeldovich code does not know about cosmology.  The default is 0,
but one almost certainly wants to change this!

Note that using this parameter means that the choice of unit of power 
in the input power spectrum file is irrelevant (but we do care about
the unit of wavenumber, see above).

`ZD_Pk_smooth`: *double*  
The length scale by which to smooth the input 
power spectrum before generating the density field.  This is applied as
a Gaussian smoothing as `exp(-r^2/2a^2)` on the density field, which is
 `exp(-k^2 a^2)` on the power spectrum.  Smoothing the power spectrum 
is useful for testing, as it reduces grid artifacts.  The smooth occurs
after the power spectrum has been normalized.  Default is 0.

`ZD_qPk_fix_to_mean`: *integer*  
Fix the amplitude of the modes to sqrt(P(k)).  The phases are unchanged.
That is, running with this option off then on will produce ICs with identical phases
but different mode amplitudes.  This is useful for producing "paired and fixed" sims,
as suggested by [Angulo & Pontzen (2016)](http://arxiv.org/abs/1603.05253). Default is 0.

`ZD_qoneslab`: *integer*  
If `> 0`, output only one PPD slab.  For debugging only.
The default is `-1`.

`ZD_qonemode`: *integer*  
If `> 0`, zero out all modes except the one with the wavevector specified in `ZD_one_mode`.

`ZD_one_mode`: *three ints*  
This is the one wavevector that will be inserted into the box if `ZD_qonemode > 0`.
This is useful for automatically iterating through a series of wavevectors, for examining isotropy or Nyquist effects, for example.
Each component can be an integer in the range `[-ppd/2,ppd/2]`.

`ZD_qPLT`: *integer*  
If `> 0`, turn on particle linear theory corrections.
This tweaks the displacements and velocities, mostly near `k_Nyquist`, to ensure everything starts in the growing mode.
The output format should include velocities if you turn this on, either in the `RVZel` or `RVdoubleZel` format.

`ZD_PLT_filename`: *string*  
The file containing the PLT eigenmodes; i.e. the true growing modes for the grid.
This file usually contains something like a 128<sup>3</sup> grid,
and the code linearly interpolates the eigenmodes and eigenvalues to finer meshes as needed.

`ZD_qPLT_rescale`: *integer*  
If `> 0`, increase the amplitude of the displacements on small scales (near `k_Nyquist`)
to preemptively compensate for future undergrowth that we know happens on a grid.

`ZD_PLT_target_z`: *double*  
If `ZD_qPLT_rescale > 0`, then increase the initial displacements such that they will match the linear theory prediction
at this redshift.  Recall that modes on the grid (mostly) grow more slowly than linear theory, which is why we (mostly) increase the initial displacements.
This redshift should be in a quasi-linear regime where linear theory is still mostly valid, e.g. `z~5`.

`ZD_k_cutoff`: *double*  
The wavenumber above which not to input any power, expressed such that `k_max = k_Nyquist / k_cutoff`, e.g. `ZD_k_cutoff = 2` means we null out modes above half-Nyquist.  Non-whole numbers like 1.5 are allowed.  This is useful for doing convergence tests, e.g. run once with `PPD=64` and `ZD_k_cutoff = 1`, and again with `PPD=128` and `ZD_k_cutoff = 2`.  This will produce two boxes with the exact same modes (although the PLT corrections will be slightly different), but the second box's modes are oversampled by a factor of two.  To keep the random number generation synchronized between the two boxes (fixed number of particle planes per block), `ZD_NumBlock` is increased by a factor of `ZD_k_cutoff`.

`BoxSize`: *double*  
This is the box size, probably in Mpc or h<sup>-1</sup>Mpc.  The zeldovich code only cares about the units to the extent that they should match the units in the power spectrum file.  See `ZD_Pk_scale` for further discussion.

`NP`: *long long int*  
Number of particles.  This must be a perfect cube for the zeldovich code to work.

`CPD`: *integer*  
Cells per dimension. Zeldovich will output CPD slabs, each with CPD^2 cells.  These slabs are only defined by the 
initial grid position; we do not guarantee that the initial displacement may not have
taken the particle out of the slab.  Usually the deviations will be small enough that 
particles move by at most 1 slab.

`InitialConditionsDirectory`: *string*  
The location to write the output.  In addition,
the zeldovich code will use this for the swap space for the block transpose.
This will generate files of the name `zeldovich.%d.%d`, which should be automatically
deleted after the code has finished.

`InitialRedshift`: *double*  
The output redshift.  This is **only used for determining rescaling amplitude**; i.e. this option has no effect if `ZD_qPLT_rescale` is not set.  This code does not compute growth functions; `ZD_Pk_sigma` controls the power spectrum normalization.

`ICFormat`: *string*  
Valid options are: `RVZel`, `RVdoubleZel`, or `Zeldovich`.

One should use one the `RV` options if `ZD_qPLT` is set, because the velocities have been explicitly computed in Fourier space.

All displacements are comoving displacements in the same units as `BoxSize`, and the velocities are comoving redshift-space displacements (same units as `BoxSize`).
To get to physical velocities from comoving redshift-space displacements, multiply by a*H(z).

The comoving positions of the initial lattice (in unit-box units where the domain is [0,1)) are simply given by
```
x = i/PPD
y = j/PPD
z = k/PPD
```
Thus, the global, absolute positions of the particles can be formed by adding the displacements (which this code outputs) to this lattice.

- `RVZel`  
Single-precision output of the displacement and velocity, and particle lattice location.
```
class RVZelParticle {
public:
    unsigned short i,j,k;
    float displ[3];
    float vel[3];
};
```

- `RVdoubleZel`  
Double-precision output of the displacement and velocity, and particle lattice location.
```
class RVdoubleZelParticle {
public:
    unsigned short i,j,k;
    double displ[3];
    double vel[3];
};
```

- `Zeldovich`  
Double-precision output of the displacement and particle lattice location.
```
class ZelParticle {
public:
    unsigned short i,j,k;
    double displ[3];
};
```


## License
[MIT](LICENSE)

If you use this code, please cite [Garrison et al. (2016)](https://arxiv.org/abs/1605.02333).

