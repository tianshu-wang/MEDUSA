## To-do

- Check Defaults
- GW cadence defaults, nstep_log, nstep_analysis

## Branches:

The `master` branch is, and shall always be, *the* production branch.
Everything in `master` should have been tested, and be guaranteed to work.
For the most part, changes made to `master` should always be done via pull
requests, primarily from the `dev` branch.

The `dev` branch is the main developmental branch. This branch should contain
new features that are expected to work but may be untested or only tested in
simplified cases.

All experimental features should be in their own branch that was copied from
either the `master` or `dev` branches.

## Recent changes

1) C-only code. Default: On. To turn it off: set `USE_FORTRAN_CODE=TRUE` in one of
   the make files. Some Fortan file names now have a `_F90` suffix.

- `eos/burrows/collapse.f90` and `eos/burrows/eos_stuff.f90` have been combined into
  `eos_stuff.c`
- `opacities/opacities_table_module.f90` has been replaced by `opacities_table_module.c`


Other files that call the Fortran function are modified correspondingly. Notice that in
`problems/`, only `init_rad_ccsn.c` is modified. For other problems, the user has to
modify `init_*.c` so that all Fortran calls are replaced by C functions with the same
name (C is case sensitive while Fortran is not, so be careful!).

This branch is checked in 1D. The result shows no difference from the Fortran version
until machine precision. In higher dimension, there are some deviations after chaos
developed. Need further check.

2) Interpolate the exponentials. Default: Off. To turn it on: set `USE_EXP_INTERP` to
   `TRUE` in a make file.

3) Kompaneets scheme of neutrino inelastic scattering on nucleons. Default: On.
   To turn it off: in the input file of Fornax, add `use_kom=0`.

4) Higher resolution at pns surface. Default: On. To turn it off: in setup.mk
   set `RCOORD=SINH`.

5) Linear memory allocation. Default: On. To turn it off: set
   `USE_LINEAR_ALLOCATION=FALSE` in a make file.

6) Larger time step method. Default: Off. To turn it on: set `USE_LARGER_STEP=TRUE`.

7) Others: MPI can be turned off for debugging. In `opacities/opac_emis_neutrino.c`, a bug
   in the function `inelastic_opac_emis()` is fixed. In `step.c`, `io.c` and `mpi.c`, some
   sentences use `sim.p[ii]` to get the ii-th element of sim.p (in 1D) or use
   `sim.p[ii][jj]` to get the (ii,jj)-th element of sim.p (in 2D). They are now replaced with macros.

## Pre-run Checklist

After compiling and before starting a new run. One should go through
[CHECKLIST.md](CHECKLIST.md).

```sh
source ./machines/della.env
```

Finally, you can compile with

```sh
make clean; make -j
```
### Common Options:

* Gravity

  In `setups.mk`, set `GRAV=SPHERICAL_MONOPOLE_GRAV` or `GRAV=SPHERICAL_MULTIPOLE_GRAV` for
  either monopole-only or multipole gravity in spherical geometry (note that `GEOM=SPHERICAL`
  must also be defined).  Then, for multipole, set the input file parameter `multipole_lmax=???`
  to set the maximum multipole order (no more than 12 is recommended).  $\Phi$

* Perturbations

  1. Perturbations in the style of [Muller & Janka (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.448.2141M)
  rely on spherical harmonics.  Rather than compute these by hand, we use functions from the GNU Scientific
  Library (GSL).  In `setups.mk`, set `PERTURB=VELOCITY_SPHERICAL_HARMONIC` and ensure that both `libgsl`
  and `libgslcblas` are in your `PATH` (either as `.so` or `.a`), or set the environment variable using
  `GSL_DIR=/path/to/libgsl` (this is done automatically on, say, `cori` by issuing `module load gsl`).
  The `Makefile` will handle the rest.

  2. Random density perturbations in the style of [Summa et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...825....6S)
  can be introduced by setting `PERTURB=DENSITY_RANDOM` in `setups.mk`.  This will cause a small-amplitude
  random density perturbation to be applied 10ms after bounce of the form <a href="https://www.codecogs.com/eqnedit.php?latex=\rho'&space;=&space;\rho&space;\left(1&space;&plus;&space;A&space;r&space;\right&space;)," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho'&space;=&space;\rho&space;\left(1&space;&plus;&space;A&space;r&space;\right&space;)" title="\rho' = \rho \left(1 + A r \right )" /></a>,
  where A is the amplitude (in percentage) and r is a random number between -1 and 1.  Use the input file
  variable `perturb_level=???` to set the amplitude.

### Compiler flags:

* OpenMP:

  To compile for OpenMP, set `USE_OMP=TRUE` in `machine.mk` and, if needed, reset the platform-dependent compiler
  flag to whatever it needs be.  For example, the default flag is `-fopenmp`, but cori requires instead `-qopenmp`, hence
  [cori-intel-knl.mk](machines/cori-intel-knl.mk) resets the flag using `OMPFLAGS = -qopenmp`.  The number of threads
  per MPI task must be set via the environment variable `OMP_NUM_THREADS=???` at run time.  Often this means you
  should set this variable in your batch submission script just prior to running the executable.  The value of this
  variable is echoed in the log file at startup.

* KNL Nodes:

  To compile for Knight's Landing (KNL) nodes, one typically includes a platform-dependent compiler flag to indicate the
  use of the `AVX512` instruction set.  See [`cori-intel-knl.mk`](machines/cori-intel-knl.mk) for an example.  Additionally,
  to achieve per-node performance comparable to or exceeding that of Haswell nodes, compile with OpenMP and use at
  least 2 threads per MPI task.

## Running

### 3-d Runs:

Since the pre-bounce epoch is largely unaffected by multidimensionality, it may help to run in
1-d first, say, to 10ms after bounce.  To do this, first compile with `NDIM=1` (e.g., use the
setup file [`rad_ccsn_1d.mk`](setups/rad_ccsn_1d.mk)), and in the input file, set

```sh
restart_from_1d = 0
dump_hdf = 0
model_name = s15.wh07.fornax
```

This will cause Fornax to create 1-d ASCII dump files.  If `USE_AUX_SRC=TRUE` is used in `setup.mk`,
then [`init_rad_ccsn.c`](problems/init_rad_ccsn.c) will compute the time of bounce.  You can then
take an ASCII dump at, say, 10ms after bounce, copy it to the working directory, and rename it
something informative.  Then, recompile with `NDIM=2` or `NDIM=3`
(e.g., use the setup file [`rad_ccsn_3d.mk`](setups/rad_ccsn_3d.mk)), and in the input file, set instead

```sh
restart_from_1d = 1
dump_hdf = 1
model_name = s15_wh07_n608_20g_10ms_1d
```

This will cause Fornax to restart from the 1-d ASCII dump file and create multi-d HDF5 dump files henceforth.
Optionally, the 1-d evolved model can be perturbed before resuming the calculation in multi-d.

Note. After the first segment of a 3-d run it is necessary to set `restart_from_1d` to `0` in order
to be able to recover from a 3-d restart file.

### MPI Domain Decomposition:

The multi-d domain decomposition is done using a type of greedy algorithm.  For large MPI task counts,
this can be somewhat of a slow process.  Fornax can read/write a computed decomposition to a shareable
ASCII file to avoid having to do this more than once for a given MPI task count, resolution, and
dimensionality.  To load a decomposition file, set the input file parameter `decomp_from_file=1`
and place the file in the working directory.  When not loading a decomposition file (default),
Fornax will (re)compute the decomposition and attempt to write it to a file.  However, if the file
already exists in the working directory, overwriting will be prohibited and you will be forced to
explicitly delete the old file first.

Currently, the "efficiency" of the decomposition measures how well it balances the computational
load across MPI tasks, defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon&space;=&space;1&space;-&space;\left(&space;\frac{N_\mathrm{worst}&space;-&space;N_\mathrm{ideal}}{N_\mathrm{ideal}}&space;\right&space;)." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon&space;=&space;1&space;-&space;\left(&space;\frac{N_\mathrm{worst}&space;-&space;N_\mathrm{ideal}}{N_\mathrm{ideal}}&space;\right&space;)." title="\epsilon = 1 - \left( \frac{N_\mathrm{worst} - N_\mathrm{ideal}}{N_\mathrm{ideal}} \right )." /></a>

The efficiency will be 1 for a perfectly-balanced decomposition, but it could be 0 or even negative
for an inefficient decomposition.  Moreover, the efficiency of Fornax's decompositions is neither an
analytic nor monotonic funciton of the number of MPI tasks, nor is it necessarily optimal.  Finally,
the decomposition does not currently account for the connectivity of the MPI task domains, so in 3-d,
domains along the polar axes may share a communication boundary with many other domains compared to,
say, the global average.

The decomposition is done in serial.  To compute several decompositions simultaneously, you can set the following
input file variables:

```sh
decomp_only=1
decomp_link=1
decomp_npmin=10000
decomp_npmax=100000
decomp_npskip=1000
```

In this example, Fornax will _only_ compute the decompositions for the number of MPI tasks in the range from
10,000 to 100,000 in increments of 1000.  If you are only interested in the efficiencies, set `decomp_link=0`
to avoid the costly computation of neighbor adjacency, and the efficiencies will be written to the file
`efficiency.txt` in the working directory.  Otherwise, if `decomp_link=1`, the full decomposition will be computed
and the results saved to a file.  If `USE_MPI==TRUE` is set in `machine.mk`, then each MPI task will compute a
separate decomposition in the sequence.

### Parameter File:

It is highly recommended to look at the `parse_input()` function in [`io.c`](io.c) for a complete list of input file parameter
options and their default values.

In addition to the Fornax executable, to run you will need a parameter file, a
progenitor model and, if you are using HPC resources, a batch script. Examples
can be found in the `./problems`, `./progenitors`, and `./machines` folders.

The recommended way to run is with the
[batchtools](https://bitbucket.org/dradice/batchtools) package. You should
download batchtools and add it to your PATH

```sh
export PATH=/path/to/batchtools/bin:$PATH
```

The recommended workflow with batchtools is

1. Create a simulation directory

```sh
mkdir /path/to/mysim && cd /path/to/mysim
```

2. Initialize the simulation using batchtools

```sh
batchtools init --exe /path/to/fornax/fornax \
                --parfile /path/to/fornax/problems/input.rad_ccsn \
                --batch /path/to/fornax/machines/della.sub \
                --include /path/to/fornax/setup.mk \
                --include /path/to/fornax/machine.mk
```

3. Edit the ```BATCH/CONFIG``` file to configure your run

4. Create a segment and run with

```sh
batchtools makesegment && batchtools submit
```

### Timing

It is possible to collect detailed timing statistics during the run.

Timing needs to be activated in the setup.mk

```sh
DO_TIMING=TRUE
```

With default settings timing statistics are only saved at the end of the run,
but it is possible to enable more frequent output by setting the input file
option `nstep_timing`.

Timing of each routine is preferentially enabled using the TIMER_START and
TIMER_STOP macros. If you need to enable more than one timer in a routine see
the example in the `main` function.

The timing routines are thread safe, but they will only collect timing for the
root thread. Per-thread timing is not supported at this time.
