## MPI zeldovich-PLT 
Generating particle initial conditions using zeldovich-PLT with MPI.

### Building

- **Dependencies**
  - An MPI C/C++ toolchain (e.g. `mpicc`, `mpicxx` or Intel MPI).
  - FFTW3 built with OpenMP support (single or double precision)

- **Environment**
  - Set `FFTW_DIR` to the root of your FFTW installation (it must contain `include/` and `lib/`).
  - Optionally set `COMPILER_FAMILY` to `gcc` (default) or `intel`.
  - Optional `CFLAGS` are passed through to the compiler; define `-DUSE_DOUBLE_PRECISION` here to link against double-precision FFTW.

- **Build command**

```bash
cd zeldovich-PLT
make FFTW_DIR=/path/to/fftw          # add COMPILER_FAMILY=... and CFLAGS=... as needed
```

This produces the `mpi_zeldovich` executable in this directory.

### Running

Prepare a parameter file (e.g. `param_PPD16_CPD2.par`) and run the code under MPI on chosen number of MPI ranks in z and x (grid decomposition). Output is written under the directory specified by `InitialConditionsDirectory` in the parameter file (e.g. `ic/`).

### Parameter file overview

Most parameters are inherited directly from the upstream `InitialConditions/zeldovich-PLT` code; see that README for detailed cosmology and PLT options. The commonly used options here are:

- **BoxSize**: Size of the simulation box (e.g. in \(h^{-1}\,\mathrm{Mpc}\)). Must be consistent with the power spectrum units.
- **NP**: Total number of particles. Must be a perfect cube; \( \mathrm{PPD} = \mathrm{NP}^{1/3} \).
- **CPD**: Cells per dimension. The code outputs `CPD` x-slabs, each nominally containing `CPD^2` particles.
- **ICFormat**: Output particle format (`RVZel`, `RVdoubleZel`, or `Zeldovich`), as in the upstream code.
- **InitialConditionsDirectory**: Directory where IC files and temporary transpose files are written (e.g. `ic`).
- **InitialRedshift**: Initial redshift; used for PLT rescaling and \(f_{\mathrm{NL}}\) if enabled.

Cosmological and PLT-related quantities (e.g. `ZD_Seed`, `ZD_Pk_*`, `ZD_qPLT`, `ZD_PLT_*`, `ZD_f_NL`, etc.) behave exactly as described in `InitialConditions/zeldovich-PLT/README.md`.

#### New MPI grid parameters

2D MPI process grid in the \(x\)–\(z\) plane:

- **`ZD_grid_x`**: Number of MPI ranks along the global \(x\)-direction.
- **`ZD_grid_z`**: Number of MPI ranks along the global \(z\)-direction.

The total number of MPI ranks you launch should be compatible with this decomposition (and with any implicit decomposition in the remaining direction). These parameters control how slabs are distributed across ranks:

- **1D outputs**: Set `ZD_grid_z = 1`. All slabs share a single \(z\)-rank and the output directory layout is effectively 1D in \(x\).
- **2D outputs**: Require `ZD_grid_z > 1`. Each \(z\)-rank writes its own subset of slabs in separate subdirectories (see below).

### IC Output format

In **Mode 3**, the code writes per–x-slab particle ICs. The output directory structure depends on the \(z\)-grid:

- **`ZD_grid_z == 1`**
  - All slabs are written directly under the IC directory:
    - `ic/ic_{slab:04d}`

- **`ZD_grid_z > 1`**
  - Each MPI rank in \(z\) writes to its own subdirectory:
    - `ic/z{MPI_rank_z:02d}/ic_{slab:04d}_z{MPI_rank_z:02d}`

Here `slab` is the x-slab index, and `MPI_rank_z` is the rank’s index in the \(z\)-direction of the MPI grid.

#### Indices within files

Each particle record in these files contains lattice indices:

- `i` and `j` are **global** grid indices in \(x\) and \(y\).
- `k` is **slab-local** in the rank’s portion of \(z\).

To convert `k` to a global \(z\) index, you must combine it with the slab and the rank’s \(z\)-range (i.e. use the slab index and `MPI_rank_z` to reconstruct the global \(k\)).

Records in a file are written in the following order:

1. Iterate over all global \(z\) values within the rank’s \(z\)-range,
2. For each \(z\), iterate over all global \(y\),
3. For each \((z, y)\), iterate over slab-local `k`.

This ordering should be taken into account when post-processing or converting the ICs.

