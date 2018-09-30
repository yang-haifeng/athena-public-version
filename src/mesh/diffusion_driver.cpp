//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file diffusion_driver.cpp
//  \brief sets the control parameters for hydro and field diffusion modules

// C/C++ headers
#include <iostream>   // endl
#include <cfloat>     // FLT_MAX
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "mesh.hpp"
#include "diffusion_driver.hpp"

//----------------------------------------------------------------------------------------
//! \fn DiffusionDriver::DiffusionDriver(Mesh *pmy_mesh, ParameterInput *pin)
//  \brief DiffusionDriver constructor
DiffusionDriver::DiffusionDriver(Mesh *pmy_mesh, ParameterInput *pin)
{
  pm = pmy_mesh;
  diffusion_defined  = false;
  operator_split     = false;
  super_timestepping = pin->GetOrAddInteger("time","super_timestepping",0);
    
  for (int i=0; i<NPHYS; ++i){
    phys_def[i]=false;
    operator_split_def[i]=false;
  }
    
  // set up microphysics modules
  if (pin->GetOrAddReal("problem","nuiso",0.0) > 0.0){
    phys_def[ISO_VIS]=true;
  }
  if (pin->GetOrAddReal("problem","kappa",0.0) > 0.0){
    phys_def[ISO_COND]=true;
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    if ((pin->GetOrAddReal("problem","eta_o",0.0) > 0.0) ||
        (pin->GetOrAddReal("problem","eta_a",0.0) > 0.0)) {
      phys_def[OAD]=true;
    }
  }
  // more to be added as more physics implemented
    
  // set up operator split control
  if (pin->GetOrAddInteger("problem","operator_split_viscosity",0) != 0){
    operator_split_def[ISO_VIS]=true;
    operator_split_def[ANI_VIS]=true;
  }
  if (pin->GetOrAddInteger("problem","operator_split_conduction",0) != 0){
    operator_split_def[ISO_COND]=true;
    operator_split_def[ANI_COND]=true;
  }
  if (pin->GetOrAddInteger("problem","operator_split_cooling",0) != 0){
    operator_split_def[COOL]=true;
  }
  if (pin->GetOrAddInteger("problem","operator_split_field_diffusion",0) != 0){
    operator_split_def[OAD]=true;
  }
  if (pin->GetOrAddInteger("problem","operator_split_hall",0) != 0){
    operator_split_def[HALL]=true;
  }
    
  for (int i=0; i<NPHYS; ++i) {
    diffusion_defined = diffusion_defined || phys_def[i];
    operator_split    = operator_split || (phys_def[i] && operator_split_def[i]);
  }
  
  // set up super-timestepping control
  if (operator_split == false) {
    if (super_timestepping == true)
      std::cout<<"Warning in DiffusionDriver: no operator split!" << std::endl
               <<"Reset super timestepping to false."<<std::endl;
      super_timestepping = false;
  }
    
  ndiffmax = pin->GetOrAddInteger("time","ndiffmax",12);
  if (ndiffmax > 20) {
    std::cout<<"Warning: reset ndiffmax to 20."<<std::endl;
    ndiffmax = 20;
  }
  for (int i=0;i<=20;++i) {
    Real nu     = CalculateNuForSTS(i+1);
    dt_ratio[i] = CalculateDtRatioForSTS(i+1, nu);
  }
    
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real DiffusionDriver::UpdateDtAndDiffDt(Real dt, Real min_dt[NPHYS+1])
//  \brief Calculate the diffusion dt and update mhd dt accounting for operator split and
//         super timestepping options
Real DiffusionDriver::UpdateDtAndDiffDt(Real dt, Real min_dt[NPHYS+1])
{
  Real mydt = dt;
  diff_dt = (FLT_MAX);
  for (int i=0;i<NPHYS;++i) {
    phys_dt[i] = min_dt[1+i];
    if ((phys_def[i] && !operator_split_def[i]) || (i==HALL))
      mydt = std::min(mydt,min_dt[1+i]);
    if ((phys_def[i] && operator_split_def[i]) && i!=HALL)
      diff_dt = std::min(diff_dt,phys_dt[i]);
  }
    
  if (operator_split) {
    if(super_timestepping){
      Real myratio = mydt/diff_dt;
      ndiffstep = 1;
      while ((myratio > dt_ratio[ndiffstep-1]) && (ndiffstep <= ndiffmax))
        ndiffstep++;
            
      if (ndiffstep==1) {
        diff_dt = mydt;
        nu_sts  = 0.0;
      } else {
        if (ndiffstep > ndiffmax) {
          ndiffstep = ndiffmax;
          mydt = dt_ratio[ndiffstep-1] * diff_dt;
        }
        else
          diff_dt = dt/dt_ratio[ndiffstep-1];
                
        nu_sts = CalculateNuForSTS(ndiffstep);
      }
    }
    else // SUPER_TIMESTEPPING
    { // subcycling
      ndiffstep = (int)(dt/diff_dt)+1;
      if (ndiffstep == 1)
        diff_dt = mydt;
      else if (ndiffstep > ndiffmax)
        mydt = ndiffmax*diff_dt;
      else
        diff_dt = mydt/ndiffstep;
    }
  }
  else // OPERATOR_SPLIT
  {// embedded to integrator
    ndiffstep = 1;
    dt=std::min(diff_dt,mydt);
    diff_dt = mydt;
  }
    
  return mydt;
}

//----------------------------------------------------------------------------------------
//! \fn void DiffusionDriver::CalculateSubTimeStep(void)
//  \brief At each step, calculate dt for individual substeps
void DiffusionDriver::CalculateSubTimeSteps()
{
  if (super_timestepping) {
    for (int n=0; n<ndiffstep; ++n)
      dt_sub[n] = diff_dt/(1.0+nu_sts-(1.0-nu_sts)*cos(0.5*PI*(2.0*n+1.0)/(Real)(ndiffstep)));
        
      std::cout<<"Super timestepping with "<<ndiffstep<<" substeps."<<std::endl;
  }
  else // subcycling
  {
    for (int n=0; n<ndiffstep; ++n)
      dt_sub[n] = diff_dt;
      std::cout<<"Subcycling with "<<ndiffstep<<" substeps."<<std::endl;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real CalculateNuForSTS(int ndiffstep)
//  \brief Calculate the nu parameter (free parameter between 0 and 1) for super
//         timestepping. In practice, we choose it to be 1/4N^2, but this may be changed.
Real DiffusionDriver::CalculateNuForSTS(int ndiffstep)
{
  return 0.25/SQR((Real)(ndiffstep));
}

//----------------------------------------------------------------------------------------
//! \fn Real DiffusionDriver::CalculateDtRatioForSTS(int ndiffstep, Real nu)
//  \brief Calculate the ratio of super timestep Delta t to diffusion timestep diff_dt
//         for super timestepping
Real DiffusionDriver::CalculateDtRatioForSTS(int ndiffstep, Real nu)
{
  Real nu_sqrt = sqrt(nu);
  return  0.5*((Real)(ndiffstep))/nu_sqrt
          /(pow(1.0+nu_sqrt,2.0*ndiffstep)+pow(1.0-nu_sqrt,2.0*ndiffstep))
          *(pow(1.0+nu_sqrt,2.0*ndiffstep)-pow(1.0-nu_sqrt,2.0*ndiffstep));
}
