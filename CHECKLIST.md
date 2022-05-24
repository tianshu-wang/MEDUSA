# Pre-run Checklist

This is the checklist that should be followed before starting a new simulation
from scratch. This file does not include help regarding compilation; for that,
see the [README.md](README.md) file, but please make sure you compile correctly
(e.g. ensure that 1D/2D/3D, multipole/monopole gravity are properly set).

## Run directory

You need to setup the directory that you plan to run your new simulation in.
If you are in an HPC environment then you want to have this directory on
"scratch" or equivalent the high IO efficiency filesystem. You may also need to
stripe the directory appropriately.  The stripping of the run directory on the
following clusters (along with their recommended stripping) have been found to
be particularly important for IO efficiency.

- Theta: `lfs setstripe -c 32 <directory>`
- Frontera: `lfs setstripe -c 16 <directory>`

Can someone check the above?

## Input File

To start a simulation run you must provide Fornax with an input file, which
often have a name similar to `input.rad_ccsne`
(see [input.rad_ccsne.example](input.rad_ccsne.example) for an example).
The sections following this one outline specific choices that need to be made
which are set in the input file.

## Restarts and Initial Conditions

The below settings are used to specify the initial conditions of a run and/or
how to restart a run.

    # Prepare initialization/restart data
    restart_from_1d   = 0   # flag to restart from 1-d ASCII file
    restart_from_3d   = 0   # flag to restart from a 3-d HDF5 dump
    restart_from_last = 0   # flag to ignore input file and use the most recent restart dump instead
    restart_from_hdf  = 0   # flag to restart from HDF5 dump
    #restart_delay = 0.010  # post-bounce 1-d ASCII restart file delay [s]
    model_name = <path_to_progenitors>/s9.0.swbj15.fornax  # 1D progenitor/IC data file

## EOS and Opacities

EOS data is located
[here](https://www.astro.princeton.edu/~burrows/josh/radice/EOS),
and Opacity data is located
[here](https://www.astro.princeton.edu/~burrows/josh/radice).

Currently (as of 2021/09/15) the "goto" EOS is
[SFHo](https://www.astro.princeton.edu/~burrows/josh/radice/EOS/table30030050.stellarcollapse.lr2_15.5.Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.dat) ([Steiner et al. 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...774...17S/abstract))
but using other EOSes may also be useful. Note that currently the parameters
for the EOS table (e.g. density range) are hard coded into
[collapse.f90](eos/burrows/collapse.f90) and [eos_stuff.c](eos/burrows/eos_stuff.c) and needs to be updated when changing
the EOS file.

There are two opacity files that are needed: the parameter file
(typically of the form `opacbin.<name>.param`) and the binary data file
(`opacbin.<name>.bin` or `opacity.<name>.bin`). Note that some of the parameter files can be used for multi binary data files.

Currently (as of 2021/09/15) the "goto" opacity data is
[opacity.SFHo.juo.horo.brem1.extendedT.bin](https://www.astro.princeton.edu/~burrows/josh/radice/opacity.SFHo.juo.horo.brem1.extendedT.bin)
which uses the "extended T" (T=temperature) parameter file
[opacbin.extendedT.param](https://www.astro.princeton.edu/~burrows/josh/radice/opacbin.extendedT.param).
While changing the opacity can be useful, the "extended T" files are preferred.

The filenames for the opacity data typically describe the physical effects
accounted for within the file:
- `SFHo` (EOS; should be the same as the above EOS)
- `juo` ([Juodagalvis et al. 2010](https://ui.adsabs.harvard.edu/abs/2010NuPhA.848..454J/abstract) electron capture rates)
- `horo` ([Horowitz et al. 2017](https://ui.adsabs.harvard.edu/abs/2017PhRvC..95b5801H/abstract) manybody correction)
- `brem1` (1 times standard bremsstrahlung)


Example settings for the input file:
    # opacity and EOS files
    opac_param_file = <path_to_data>/opacbin.extendedT.param
    opac_file = <path_to_data>/opacity.SFHo.juo.horo.brem1.extendedT.bin
    eos_file = <path_to_data>/table30030050.stellarcollapse.lr2_15.5.Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.dat
    include_inelastic = 1  # include inelastic scattering
    inelastic_root = <path_to_data>/300.mixed/  # directory of inelastic scattering data

    num_comp = 1  # number of compositions to track (default=1 for Ye)


## Grid and Resolution

Below are the standard settings for 3D grids:

    # grid/resolution settings
    n1  = 1024            # radial gridzones
    n2  = 128             # theta gridzones
    n3  = 256             # phi gridzones
    outer_radius = 3.0e9  # max radius [cm]
    ccsn_dr_min = 5.0e4   # smallest dr (in core) [cm], best to leave 5e4 [cm] = 0.5 [km]
    grid_units = km       # used for some outputs

## Perturbations

Velocity perturbations have been shown to make a difference in explodability of
CCSN models (and I'm too lazy to find a reference for that). It is important to
specify physically motivated perturbations (e.g. around the Si/O interface)
that are useful (e.g. do not get accreted onto the PNS very early). The lmax
used can make a big difference in outcome with lower lmax leading to more
explodable models. Below is an example of how to set the perturbations in one
region, however up to three different perturbations can be specified by
replacing the `1` in the keyword with the correct region number. Please note
that the radial range in this example may not be valid for your run.

    # Initial perturbations settings (up to 3 regions)
    perturb_r1m = 200e5   # velocity perturbation, region 1, rmin [cm]
    perturb_r1p = 1000e5  # velocity perturbation, region 1, rmax [cm]
    perturb_l1  = 10      # velocity perturbation lmax, region 1
    perturb_m1  = 10      # velocity perturbation mmax, region 1
    perturb_n1  = 4       # velocity perturbation nmax, region 1
    perturb_dv1 = 1.0e7   # velocity perturbation amplitude [cm/s], region 1

## Radiation

Below are the standard settings for radiation. Note that opacity and inelastic
scattering data exists for only a few possible choices for the below settings.
Do not modify these settings unless you know there is corresponding opacity
and inelastic scattering data, and you have properly specified them in the
above step.

    # radiation settings
    # Do not modify this unless you know what your doing and have the appropriate opacity data
    freq_type = mev  # unit for energy bins
    nr1 = 12         # number of energy bins for species 1 (electron neutrinos)
    nr2 = 12         # number of energy bins for species 2 (anti electron neutrinos)
    nr3 = 12         # number of energy bins for species 3 (other neutrinos)
    emin1 = 1.       # smallest energy for species 1 (electron neutrinos)
    emax1 = 300.     # largest energy for species 1 (electron neutrinos)
    emin2 = 1.       # smallest energy for species 2 (anti electron neutrinos)
    emax2 = 100.     # largest energy for species 2 (anti electron neutrinos)
    emin3 = 1.       # smallest energy for species 3 (other neutrinos)
    emax3 = 100.     # largest energy for species 3 (other neutrinos)

force_flux_limit = 1  # force the energy density to be >= 0 and the flux to be <= cE
use_chat = 0  # Flag to use a reduced speed of light (best to leave this zero)

## Decomposition

Check the directory of decomposition files for 3D simulations, and calculate
new ones (separately) IF necessary.  Any new gridding/resolution requires new
decomposition files. These settings can be specified as below:

    # decomposition file settings
    decomp_only = 0       # flag to test decomposition only without allocating or running
    decomp_from_file = 0  # flag to load decomposition proc info from file ("proc_info.txt")
    #decomp_path = <directory_path>  # directory that "proc_info.txt" is located in

## Tracers

To add tracers to a run set `dt_pdump` (dt between tracer dumps),
`n_tracer_target` (number of tracers you want), mass_inside_tracers (tracers
will fill the volume outside of this interior mass, in solar masses). E.g.:

    # tracer settings
    dt_pdump = 0.001            # tracer dump frequency, 1 ms
    mass_inside_tracers = 1.35  # place tracers exterior to this mass cutoff
    n_tracer_target = 1000000   # place this many tracers (defaults to zero)

## Rotation

One can enable initial rotation of the progenitor with the below settings.
The parameter `rotate_Omega0` sets the core/central rotation rate in radians
per second, and `rotate_A` sets the (cylindrical) falloff radius in km for the
rotation profile; far outside this radius there is no rotation.

    # initial rotation settings
    rotate_Omega0 = 0.1  # core rotation rate [rad s^-1] (defaults to zero)
    rotate_A = 10000     # rotation falloff radius [km]

## Other Settings

Multipole Gravity (make sure the code is compiled with this feature if you want
to use it):

    # gravity max multipole
    multipole_lmax = 12

## Backwards Compatibility

This section contains information on migrating some data/setting from old
versions of Fornax, and can be skipped if not dealing with Fornax data made in
this decade (2020s).

If you want to use the old restart format for reading/writing, you need to open
io.c, go to the functions restart_dump and restart_read and change, e.g.,
restart_read_v2 (newformat) to restart_read_v1 (old format). I would suggest to
do this for the first time when migrating a run from the old restart format to
the new one.

## Analysis (Post-processing)

Check `config.py` to ensure `tbounce.prescribed` is set correctly.
