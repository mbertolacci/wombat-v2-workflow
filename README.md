# Inferring changes to the global carbon cycle with WOMBAT v2.0, a hierarchical flux-inversion framework

This repository contains code to reproduce the results in the paper [Inferring changes to the global carbon cycle with WOMBAT v2.0, a hierarchical flux-inversion framework](https://arxiv.org/abs/2210.10479). This README assumes familiarity with the paper.

The code implements the WOMBAT v2 framework in the context of the CO2 application described in the paper. With some work, this can also be modified to estimate fluxes for another trace gas.

Unless stated otherwise, all commands are to be run in the root directory of the repository.

# Accessing the estimated fluxes and parameters

The outputs of WOMBAT v2.0 for the application to CO2 are [available on Zenodo](https://doi.org/10.5281/zenodo.8015839). These comprise samples from the posterior distribution of the model parameters, bottom-up estimates and samples from the posterior distribution of the fluxes, and bottom-up estimates and samples from the posterior distribution of the trend/seasonality terms.

# Installation/setting up an environment

This workflow requires R version 4+, along with a variety of dependency packages in both languages. The easiest way to set up an environment in which to run this code is to use [Anaconda](https://www.anaconda.com/). Instructions for setting up an appropriate conda environment are provided below. If, for some reason, you don't want to do that, you can adapt the instructions below to a local installation using a local copy of R - everything should work okay. For what remains, we assume you have Anaconda installed.

Create a conda environment and set it to use [conda-forge](https://conda-forge.org/):

```
conda create --yes --prefix .conda_env
conda activate ./.conda_env
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
```

Install conda packages

```
conda install r-base pkg-config udunits2 libgdal netcdf4 cdo nco
```

Set your CRAN mirror. WOMBAT was made in Australia, so we choose the mirror run by CSIRO, but you can pick your favourite:

```
cat > .conda_env/lib/R/etc/Rprofile.site <<- EOF
local({
  r <- getOption('repos')
  r['CRAN'] <- 'https://cran.csiro.au'
  options(repos = r)
})
EOF
```

Install R dependencies (this might take awhile):

```
Rscript -e "renv::restore()"
```

Check out the required version of GEOS-Chem:

```
(cd external/GEOS_Chem && git clone --branch 12.3.2 https://github.com/geoschem/geos-chem.git Code.v12.3.2)
```

# Getting data

All input datasets go into the `data` directory. There are a few files already there, but the rest will need to be retrieved from their primary sources.

## Data required for all steps

In addition to those data sets already provided, the following data are needed to run the inversions:

- Available from the [OCO-2 v10 MIP website](https://gml.noaa.gov/ccgg/OCO2_v10mip/download.php):
  - 10s averages of the OCO-2 B10 retreivals, in a file named `OCO2_b10c_10sec_GOOD_r5.nc4`
  - The ObsPack measurements used by the v10 MIP, `obspack_co2_1_OCO2MIP_v3.2_2021-05-20`. These should be untarred into the data directory.
  - The TCCON retrievals in the file `downloaded_20211217.tgz`. These were not actually used in the paper, but are required for the matching step to work.
- The Landschutzer ocean fluxes in the file `spco2_MPI-SOM_FFN_v2020.nc`, available from this [NOAA website](https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0160558/MPI_SOM-FFN_v2020/)
- Three additional data sets are needed that are not generally available at this time. These can be provided upon request, and we will endeavour to make them freely available. These are:
  - GFED4.1s fire emissions, preprocessed for ingestion by HEMCO
  - The SiB4 bottom-up estimates of GPP and respiration
  - Fossil-fuel emissions based on ODIAC postprocessed for ingestion by HEMCO.
  - Data for the Lauder CO2 collection site in NZ (not used in the project, but required for the matching step).

GEOS-Chem requires meteorological fields and CO2 emission to run. These go into the `data/GEOS_Chem` directory (if you already have some of these, you could symlink them in). There are instructions for how to download these files on the [GEOS-Chem wiki](http://wiki.seas.harvard.edu/geos-chem/index.php/Downloading_data_from_Compute_Canada).

The directories you need to download are:

```
ExtData/GEOS_2x2.5/MERRA2
ExtData/CHEM_INPUTS/MODIS_LAI_201204
ExtData/CHEM_INPUTS/Olson_Land_Map_201203
ExtData/HEMCO/CO2/v2015-04/BIOFUEL
```

## Intermediate files to reproduce just the inversion

The most computationally expensive part of the workflow are steps 1 and 2 below, where the basis-function runs are computed and post-processed. To ease reproduction of the inversion results (steps 3 and 4 below), we can provide the necessarily post-processed outputs of steps 1 and 2 sufficient to run the inversion and generate the results. These can be provided by the author upon request.

Once you have the archive, extract the files into the root directory of this repository with

```
tar xzf ~/path/to/WOMBAT_v2_CO2_intermediates.tar.gz
```

Then you can run the workflow from step 3 onwards as described below.

# Running the workflow

The workflow of this repository is split into four steps:

1. `1_transport`: Creates GEOS-Chem basis function runs. This includes setting up the inventories, creating the run directories, and setting up the configuration files. You will then need to find a way to run GEOS-Chem for each run.
2. `2_matching`: Postprocesses the run output from the previous step by extracting the modelled mole-fraction values for each basis function at each observation time and location.
3. `3_inversion`: Performs the inversions.
4. `4_results`: Summarise the results as a series of plots and other outputs. This reproduces all the figures in the paper.

## Step 1: running transport

In this step many runs of GEOS-Chem need to be completed, and is by far the most computationally intensive part of the workflow. The code creates several run directories corresponding to basis functions, and then GEOS-Chem can be run in parallel for each directory.

To set up the basis-function runs, run the command:

```
WOMBAT_LOG_LEVEL=debug make -j4 1_transport_targets
```

This will set up GEOS-Chem run directories in `1_transport/intermediates/runs` and `1_transport/intermediates/runs-r10-r15-rNZ`. In each folder, these is a base run called `base`. You need to compile the base run code compiled by changing into the directory and calling make as follows:

```
(cd 1_transport/intermediates/runs && BPCH_DIAG=n make -j24 mpbuild)
(cd 1_transport/intermediates/runs-r10-r15-rNZ && BPCH_DIAG=n make -j24 mpbuild)
```

The code for the other directories all symlink back to this directory to access the built GEOS-Chem executable.

Now you need to run all the GEOS-Chem runs in these directories, including the base run. Their outputs will require around 10TB of space. To run then, the simplest scheme is to change into a run directory and call geos.mp. For example:

```
(cd 1_transport/intermediates/runs/residual_20210301_part003 && OMP_NUM_THREADS=8 ./geos.mp)
```

After this is completed, you should run the `postprocess-run.sh` script in the same directory:

```
(cd 1_transport/intermediates/runs/residual_20210301_part003 && bash postprocess-run.sh)
```

It may however be easier to run these through a batch system like Slurm. Example scripts to do this are provided for the [Gadi supercomputer system run by the Australian NCI](https://nci.org.au/) in `run-gadi.sh`, and for a local Slurm setup in `run-niasra-hpc-sbatch.sh`. All the runs are independent from each other so they can be done in parallel.

### A note on the transport runs

When we first performed the basis-function runs, there was a bug in the region mask for a few regions. Because the cost of re-running all the runs is prohibitive, and most runs did not need to change, we added a second set of runs that corrects those regions. In this repository we have retained that "two run" structure: the first set of runs are in `1_transport/intermediates/runs`, and the second in `1_transport/intermediates/runs-r10-r15-rNZ`.

If reproducing the results from scratch or adapting the framework to a new problem, the second set of runs should be removed.

## Step 2: postprocessing runs (matching)

Similar to the runs above, the matching step first creates a directory for each run with scripts that need to be run to postprocess the GEOS-Chem outputs. Then, you run each of those scripts, potentially in parallel.

The directories can be created with

```
WOMBAT_LOG_LEVEL=debug make -j4 2_matching_targets
```

This will create a directory structure in `2_matching/intermediates/runs` and `2_matching/intermediates/runs-r10-r15-rNZ`, which parallels the structurs in `1_transport` from the previous step.

Once that's done, there are two steps to run for each directory. The first step aggregates the basis-function fluxes to a monthly resolution. Scripts for running that on the Gadi supercomputer are in `2_matching/intermediates/runs/<RUN>/run-aggregate-flux-gadi.sh`; the script just uses [CDO](https://code.mpimet.mpg.de/projects/cdo/embedded/index.html) and can be adapted to your needs. The second step extracts those portions of the mole fraction outputs of GEOS-Chem that correspond to observations. Scripts for running that on the Gadi supercomputer are in `2_matching/intermediates/runs/<RUN>/run-matching-gadi.sh`. This calls the R script `2_matching/src/match.R` and can again be adapted to your needs.

## Steps 3 and 4: inversion and results

Once steps 1 and 2 are completed (or you've downloaded the intermediate files mentioned earlier), you can run the inversions. The simplest way to do this is to run

```
WOMBAT_LOG_LEVEL=debug OMP_NUM_THREADS=8 make -j4 4_results_targets
```

This should work it's way through the steps, ultimately running the inversion and then generating all the plots and outputs. You can modify the `-j` option and the `OMP_NUM_THREADS` variable to suit your local system.
