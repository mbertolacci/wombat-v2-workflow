NOTE THE STUFF WITH THE EXTRA SET OF RUNS??? PROBABLY SHOULD. OR JUST FIX IT SILENTLY

LAUDER STUFF?

# Installation/setting up an environment

This workflow requires R version 4+, along with a variety of dependency packages in both languages. The easiest way to set up an environment in which to run this code is to use [Anaconda](https://www.anaconda.com/). Instructions for setting up an appropriate conda environment are provided below. If, for some reason, you don't want to do that, you can adapt the instructions below to a local installation using `pip` and your local copy of R - everything should work okay. For what remains, we assume you have Anaconda installed.

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

# Running the inversions

The workflow of this repository is split into four steps:

1. `1_transport`: Creates GEOS-Chem basis function runs. This includes setting up the inventories, creating the run directories, and setting up the configuration files. You will then need to find a way to run GEOS-Chem for each run.
2. `2_matching`: Postprocesses the run output from the previous step by extracting the modelled mole-fraction values for each basis function at each observation time and location.
3. `3_inversion`: Performs the inversions.
4. `4_results`: Summarise the results as a series of plots and other outputs.

# Transport

Note that there was a bug in the region mask for a few regions for the initial set of runs (1_transport/intermediates/runs)

This was corrected in a second set of runs made for those regions only (1_transport/intermediates/runs-r10-r15-rNZ)