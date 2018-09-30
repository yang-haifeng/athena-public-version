#ifndef DIFFUSION_DRIVER_HPP_
#define DIFFUSION_DRIVER_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file diffusion_driver.hpp
//  \brief defines the control parameters for hydro and field diffusion modules

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "mesh.hpp"

class DiffusionDriver {
public:
  DiffusionDriver(Mesh *pmy_mesh, ParameterInput *pin);
  ~DiffusionDriver() {};
    
  Mesh *pm;
    
  bool phys_def[NPHYS], operator_split_def[NPHYS];
  bool diffusion_defined, operator_split, super_timestepping;
  int ndiffmax;          // maximum number of diffusion substeps
  Real dt_ratio[21];     // ratio of full to diffusion timesteps

  Real phys_dt[NPHYS];   // timestep constraint for individual physics
  Real dt_sub[20];       // substeps for super timestepping
  Real diff_dt, nu_sts;  // diffusion timestep for super timestepping
  int ndiffstep;         // number of diffusion substeps
    
  Real UpdateDtAndDiffDt(Real dt, Real min_dt[NPHYS+1]);
  void CalculateSubTimeSteps();
    
  Real CalculateNuForSTS(int ndiffstep);
  Real CalculateDtRatioForSTS(int ndiffstep, Real nu);
};
#endif // DIFFUSION_DRIVER_HPP_
