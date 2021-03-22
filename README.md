# BATS_Ecomodel
I'm trying to replicate the Luo et al. (2010) ecosystem model with data assimilation, for the Bermuda Atlantic Time-Series.

## Introduction
Many of my first-attempts can be found in `/Draft`. What I did first was to take one of the FORTRAN files (bacterial component) and adapt it to MATLAB. Taking inspiration from the way the FORTRAN was written, I coded everything (parameters, state vars) as `structs`. This was a shortsighted, as MATLAB's ODE solvers only use column vectors. 
If interested in *using* these functions, start with the last-modified set of codes.

## Use of these codes

As it stands, this is a stable forward model using dummy parameters from the original paper.

1. Clone this repo and navigate to the latest update (whatever isn't in `/Draft` or `/sandbox`). 
2. I recommend also grabbing `AlbumMaps.mat` from [here](github.com/germo006/NoahMaps). The plots rely on these colormaps. 
3. Play around with `tester.m`. This is a wrapper for the respective model functions (`ode_mod_ecosys_xxx_vec.m`, `getForcing.m`, `getppbp.m`, `constructor.m`).
  * `tester.m` has an odd section where it takes the input from `constructor.m` and rips the `struct` variables apart, piecing them back together as column vectors. When you get the output and you are curious which variable is now what, look at the `xxx_names` variables. They show the original names of the inputs in the same order as they are indexed in the solver-readable vectors. 
4. If you are interested in changing more than time intervals and plots, other functions need to be changed or run.
  a. Changing the ecosystem model (the differential equations themselves) is best done by editing the original function, written with human-readable variables (`ode_mod_ecosys_xxx.m`).
    * If you have done this, copy the text of `ode_mod_ecosys_xxx.m` into the "Copy of" function and run `vectorizer.m`. This will replace the structure notation with column vectors for parameters and state variables.
    * You will then need to change the first line of the `ode_mod_ecosys_xxx_vec` function so that the name and inputs match the filename and the requirements of the wrapper. 
  b. Changing the parameter values and initial state values can be done in `constructor.m`. Simply change the value you want to alter. 
  c. Changing `getForcing.m` is up to you. I'm still working on making that one more robust to use real data. 

## What Remains

There is a lot to be done yet. As the goal of this is to bring metabolite chemistry into play and adapt an inverse-modeling approach to BATS, this repo is effectively a big testbed at the moment.

1. No data assimilation yet, only some forcing.
2. Have not constructed adjoint model for optimization. 
3. Currently working on extending to 1-D from 0-D.