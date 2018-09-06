//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file disk.cpp
//  \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//         spherical polar coordinates.  Initial conditions are in vertical hydrostatic
//         equilibrium.
//======================================================================================

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min
#include <cstdlib>    // srand
#include <cfloat>     // FLT_MIN

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../utils/utils.hpp"

// User-defined physical source term
void MySource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &w,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &u);
// User-defined Grid spacing
Real CompressedX2(Real x, RegionSize rs);

// Vector potentials for poloidal fields
static Real AphiOpen(const Real x1, const Real x2, const Real x3);
void AddPoloidalField(MeshBlock *pmb);

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

// problem parameters
static Real GM=0.0,R0=1.0,Rbuf;
static Real rho0,alpha,rho_floor;
static Real HoR0,HoRc,qT,fcool;
static Real mu,Bz0,Am0, lHoH0,res_scale;
static Real taddBp;
static int finest_lev;

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Get parameters for grav terms and radii
  GM = pin->GetOrAddReal("problem","GM", 0.0);
  R0 = pin->GetOrAddReal("problem","R0", 1.0);
  //Rdisk = pin->GetReal("problem","Rdisk");
  Rbuf = pin->GetOrAddReal("problem","Rbuf",2.0);
    
  // Get initial density
  rho0 = pin->GetReal("problem","rho0");
  alpha = pin->GetOrAddReal("problem","alpha",2.25);
  rho_floor = pin->GetReal("problem","rho_floor");
    
  // Get initial temperature (midplane and corona)
  HoR0 = pin->GetOrAddReal("problem","HoR0",0.1);
  HoRc = pin->GetOrAddReal("problem","HoRc",0.5);
  qT    = pin->GetOrAddReal("problem","qT",1.0);
  fcool = pin->GetOrAddReal("problem","fcool",1.0);
    
  // Get B field related parameters
  Bz0 = pin->GetOrAddReal("problem","Bz0",0.0);
  mu  = pin->GetReal("problem","mu");
  Am0 = pin->GetOrAddReal("problem","Am0",1.0);
  lHoH0 = pin->GetOrAddReal("problem","lHoH0",10.0);
  res_scale = pin->GetOrAddReal("problem","res_scale",0.2);
  taddBp = pin->GetOrAddReal("problem","taddBp",0.0);

  finest_lev = pin->GetInteger("problem","finest_lev");

// enroll user-defined grid spacing
  EnrollUserMeshGenerator(X2DIR, CompressedX2);

// Compute the theta-profiles for density for prescribe temperature
  int nx2bin=mesh_size.nx2;
  nx2bin = nx2bin<<finest_lev;

  // allocate user arrays to store temperature/density profiles
  AllocateRealUserMeshDataField(2);
  ruser_mesh_data[0].NewAthenaArray(nx2bin); // temperature
  ruser_mesh_data[1].NewAthenaArray(nx2bin); // density

  // allocate auxiliary arrays
  AthenaArray<Real> thetaf,thetac,lnsinthf,lnsinthc,F;
  thetaf.NewAthenaArray(nx2bin+1);
  thetac.NewAthenaArray(nx2bin);
  lnsinthf.NewAthenaArray(nx2bin+1);
  lnsinthc.NewAthenaArray(nx2bin);
  F.NewAthenaArray(nx2bin+1);

  //std::cout<<mesh_size.x2min<<'\t'<<mesh_size.x2max<<std::endl;
  for (int i=0;i<=nx2bin;++i) {
    Real rx     = (Real)(i)/(Real)(nx2bin);
    thetaf(i)   = MeshGenerator_[X2DIR](rx,mesh_size);
    if ((thetaf(i)<=1.0e-2) || (thetaf(i)>=3.13159265359))
      lnsinthf(i) = log(sin(1.0e-2));
    else
      lnsinthf(i) = log(sin(thetaf(i)));
    //std::cout<<rx<<'\t'<<thetaf(i)<<'\t'<<lnsinthf(i)<<std::endl;
  }
  for (int i=0;i<nx2bin;++i) {
    thetac(i) = ((sin(thetaf(i+1)) - thetaf(i+1)*cos(thetaf(i+1))) -
                 (sin(thetaf(i  )) - thetaf(i  )*cos(thetaf(i  ))))/
                 (cos(thetaf(i  )) - cos(thetaf(i+1)));
    lnsinthc(i) = log(sin(thetac(i)));
  }

  // compute the temperature function g(theta)=SQR(HoR_local)
  Real HoRc0=0.5;
  for (int i=0;i<nx2bin;++i) {
    Real delta = fabs(thetac(i)-0.5*PI);
    Real ext = 3.0;
    Real y   = HoR0*(HoRc0*ext-1.0)/(1.0-ext*HoR0);
    Real myHoR = std::max(HoR0,(y+HoRc0)*delta-y);
    ruser_mesh_data[0](i) = SQR(std::min(myHoR,HoRc0));
    //std::cout<<thetac(i)<<'\t'<<delta<<'\t'<<ruser_mesh_data[0](i)<<std::endl;
  }

  // now solve for F(theta)=g(theta)*f(theta) at grid interfaces
  // By definition, F(theta)=1 at theta=pi/2 (ntheta is dividable by 2)
  // Eventually, density profile f(theta) can be deduced from F and g.
  F(nx2bin/2) = SQR(HoR0);
  for (int i=nx2bin/2-1;i>=0;--i) {
    Real dlnF = (1.0/ruser_mesh_data[0](i)-(alpha+1.0))
               *(lnsinthf(i+1)-lnsinthf(i));
    F(i)      = exp(log(F(i+1))-dlnF);
  }
  for (int i=nx2bin/2+1;i<=nx2bin;++i) {
    Real dlnF = (1.0/ruser_mesh_data[0](i-1)-(alpha+1.0))
               *(lnsinthf(i)-lnsinthf(i-1));
    F(i)      = exp(log(F(i-1))+dlnF);
  }
  for (int i=0;i<nx2bin;++i) {
    ruser_mesh_data[1](i) = sqrt(F(i)*F(i+1))/ruser_mesh_data[0](i);
    //std::cout<<ruser_mesh_data[1](i)<<std::endl;
  }

  thetaf.DeleteAthenaArray();
  thetac.DeleteAthenaArray();
  lnsinthf.DeleteAthenaArray();
  lnsinthc.DeleteAthenaArray();
  F.DeleteAthenaArray();

// enroll user-defined boundary condition
  if(mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, DiskInnerX1);
  }
  if(mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X1, DiskOuterX1);
  }

// enroll user-defined physical source terms
  EnrollUserExplicitSourceFunction(MySource);

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  //AllocateRealUserMeshBlockDataField(n0+nvar);
  AllocateRealUserMeshBlockDataField(2);
// List of User field:
//     0: temperature profile (1d array)
//     1: initial density profile (1d array)

  ruser_meshblock_data[0].NewAthenaArray(block_size.nx2+2*NGHOST);
  ruser_meshblock_data[1].NewAthenaArray(block_size.nx2+2*NGHOST);

// set temperature profile
#pragma simd
  for (int j=js-NGHOST;j<=je+NGHOST;++j) {
    Real delta = fabs(pcoord->x2v(j)-0.5*PI);
    Real ext = 3.0;
    Real y   = HoR0*(HoRc*ext-1.0)/(1.0-ext*HoR0);
    Real myHoR = std::max(HoR0,(y+HoRc)*delta-y);
    ruser_meshblock_data[0](j) = SQR(std::min(myHoR,HoRc));
  }

// set initial density profile 
  int dlev = pmy_mesh->root_level+finest_lev-loc.level;
  int ifac = 1<<dlev;

  for (int j=js;j<=je;++j) {
    int j1=ifac*((j-js)+loc.lx2*block_size.nx2);
    int j2=j1+ifac-1;
    Real val=1.0;
    for (int t=j1; t<=j2; ++t)
      val *= pmy_mesh->ruser_mesh_data[1](t);
    for (int t=0; t<dlev; ++t)
      val = sqrt(val);
    ruser_meshblock_data[1](j) = val;
  }
}

//--------------------------------------------------------------------------------------
//! \fn static Real CompressedX2(Real x, RegionSize rs)
//  \brief Increase the theta grid size towards the pole
Real CompressedX2(Real x, RegionSize rs) 
// This compressedX2 doesn't respect x2min or x2max
{
  //std::cout<<"CompressedX2 called"<<std::endl;
  Real x2mid = 0.5*PI;
  Real lw, rw;

  Real xabs = fabs(x);
  if (xabs > 1.0) xabs = 2.0-xabs;

  if (xabs<1.0e-3) return 0.0;
  if (xabs>0.999)  return 3.141592653589793;
  if (xabs<=0.5) {
    Real ratn=pow(-rs.x2rat,0.5*rs.nx2);
    Real rnx=pow(-rs.x2rat,xabs*rs.nx2);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
    Real val = rs.x2min*lw+x2mid*rw;
    //std::cout<<ratn<<"\t"<<rnx<<"\t"<<std::endl;
    //std::cout<<lw<<"\t"<<rw<<"\t"<<val<<std::endl;
    return x > 0.0 ? val : -val;
  } else {
    Real ratn=pow(1.0/-rs.x2rat,0.5*rs.nx2);
    Real rnx=pow(1.0/-rs.x2rat,(xabs-0.5)*rs.nx2);
    lw=(rnx-ratn)/(1.0-ratn);
    rw=1.0-lw;
    Real val = x2mid*lw+rs.x2max*rw;
    //std::cout<<ratn<<"\t"<<rnx<<"\t"<<std::endl;
    //std::cout<<lw<<"\t"<<rw<<"\t"<<val<<std::endl;
    return x < 1.0 ? val : 2.0*PI-val;
  }
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

// First, set magnetic field
  int nx1 = block_size.nx1+2*NGHOST+1;
  int nx2 = block_size.nx2+2*NGHOST+1;

  AthenaArray<Real> A3,area,len,len_p1;
  A3.NewAthenaArray(nx2,nx1);
  area.NewAthenaArray(nx1);
  len.NewAthenaArray(nx1);
  len_p1.NewAthenaArray(nx1);

  Real betat1 = pin->GetOrAddReal("problem","betat1",0.0);
  int fieldopt = pin->GetOrAddInteger("problem","fieldopt",1);

  if (MAGNETIC_FIELDS_ENABLED)
  {
    // Set values for Aphi (vector potential)
    if (fieldopt == 1) // add open poloidal field
    {
      Real RBmin = 0.;
      Real Phimin = 2.0/(3.0-alpha)*Bz0*pow(RBmin/R0,1.0-0.5*(alpha-1.0));
      std::cout<<Bz0<<std::endl<<std::endl;

#pragma omp for schedule(static)
      for (int j=js; j<=je+1; ++j) {
        Real x2 = pcoord->x2f(j);
        Real sintheta = sin(x2);
        for (int i=is; i<=ie+1; ++i) {
          Real x1 = pcoord->x1f(i);
          Real R  = x1*sintheta;
          A3(j,i) = std::max(AphiOpen(x1,x2,0.0)*R,Phimin)/R;
	  //std::cout<<Phimin<<"\t"<<A3(j,i)<<std::endl;
      }}
    }

    if (fieldopt == 2 || fieldopt == 3) // No field
    {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          A3(j,i) = 0.0;
      }}
    }

    // Obtain poloidal B field
    for (int k=ks; k<=ke; ++k) {
#pragma omp for schedule(static)
    for (int j=js; j<=je+1; ++j) {
      pcoord->Face2Area(k,j,is,ie,area);
      pcoord->Edge3Length(k,j,is,ie+1,len);
#pragma simd
      for (int i=is; i<=ie; ++i) {
        pfield->b.x2f(k,j,i) = -(len(i+1)*A3(j,i+1) - len(i)*A3(j,i))/(area(i)+TINY_NUMBER);
	//std::cout<<pfield->b.x2f(k,j,i)<<std::endl;
      }
    }}

    for (int k=ks; k<=ke; ++k) {
#pragma omp for schedule(static)
    for (int j=js; j<=je; ++j) {
      pcoord->Face1Area(k,j,is,ie+1,area);
      pcoord->Edge3Length(k,j  ,is,ie+1,len);
      pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        pfield->b.x1f(k,j,i) = (len_p1(i)*A3(j+1,i) - len(i)*A3(j,i))/(area(i)+TINY_NUMBER);
      }
    }}

    // Set toroidal field
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real rho = rho0*pow(x1/R0,-alpha)*ruser_meshblock_data[1](j);
        Real pres= rho*(GM/x1)*ruser_meshblock_data[0](j);
        pfield->b.x3f(k,j,i) = sqrt(2.0*pres*betat1);
      }
    }}

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      const Real& b1_i   = pfield->b.x1f(k,j,i  );
      const Real& b1_ip1 = pfield->b.x1f(k,j,i+1);
      const Real& b2_j   = pfield->b.x2f(k,j  ,i);
      const Real& b2_jp1 = pfield->b.x2f(k,j+1,i);
      const Real& b3_k   = pfield->b.x3f(k  ,j,i);
      const Real& b3_kp1 = pfield->b.x3f(k+1,j,i);

      Real& bcc1 = pfield->bcc(IB1,k,j,i);
      Real& bcc2 = pfield->bcc(IB2,k,j,i);
      Real& bcc3 = pfield->bcc(IB3,k,j,i);

      const Real& x1f_i  = pcoord->x1f(i);
      const Real& x1f_ip = pcoord->x1f(i+1);
      const Real& x1v_i  = pcoord->x1v(i);
      const Real& dx1_i  = pcoord->dx1f(i);
      Real lw=(x1f_ip-x1v_i)/dx1_i;
      Real rw=(x1v_i -x1f_i)/dx1_i;
      bcc1 = lw*b1_i + rw*b1_ip1;
      const Real& x2f_j  = pcoord->x2f(j);
      const Real& x2f_jp = pcoord->x2f(j+1);
      const Real& x2v_j  = pcoord->x2v(j);
      const Real& dx2_j  = pcoord->dx2f(j);
      lw=(x2f_jp-x2v_j)/dx2_j;
      rw=(x2v_j -x2f_j)/dx2_j;
      bcc2 = lw*b2_j + rw*b2_jp1;
      const Real& x3f_k  = pcoord->x3f(k);
      const Real& x3f_kp = pcoord->x3f(k+1);
      const Real& x3v_k  = pcoord->x3v(k);
      const Real& dx3_k  = pcoord->dx3f(k);
      lw=(x3f_kp-x3v_k)/dx3_k;
      rw=(x3v_k -x3f_k)/dx3_k;
      bcc3 = lw*b3_k + rw*b3_kp1;

    }}}
  }

// Now, initialize hydro variables
  Real amp = pin->GetOrAddReal("problem","amp",0.05);
  int64_t iseed=gid;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real gc= ruser_meshblock_data[0](j); //temperature
      Real f = ruser_meshblock_data[1](j); //density
#pragma simd
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real myamp = x1 > Rbuf ? amp : 0.0;

        //Real Rcut = 10;
        //Real rho = rho0*pow(x1/R0,-alpha)*f*exp( -x1*x1/Rcut/Rcut ); 
	  // add extra exp( -(R/Rcut)**2 ) term
        Real rho = rho0*pow(x1/R0,-alpha)*f; 
	//std::cout<<f<<std::endl;

        Real cs  = sqrt((GM/x1)*gc);
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho*cs*myamp*(ran2(&iseed)-0.5);
        phydro->u(IM2,k,j,i) = rho*cs*myamp*(ran2(&iseed)-0.5);
        phydro->u(IM3,k,j,i) = rho*sqrt(GM*(1.0-(alpha+1.0)*gc)/x1);
        Real pressure = rho*SQR(cs);
        if(NON_BAROTROPIC_EOS)
          phydro->u(IEN,k,j,i)= pressure/(peos->GetGamma()-1.0)
            + 0.5*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) +
                   SQR(phydro->u(IM3,k,j,i)))/rho;
        if (MAGNETIC_FIELDS_ENABLED)
          phydro->u(IEN,k,j,i) +=  0.5*(SQR(pfield->bcc(IB1,k,j,i))+
            SQR(pfield->bcc(IB2,k,j,i))+SQR(pfield->bcc(IB3,k,j,i)));
      }
    }
  }

  A3.DeleteAthenaArray();
  area.DeleteAthenaArray();
  len.DeleteAthenaArray();
  len_p1.DeleteAthenaArray();

  return;
}

//--------------------------------------------------------------------------------------
//! \fn static Real AphiOpen(const Real x1,const Real x2,const Real x3)
//  \brief Aphi: 3-component of vector potential for open poloidal fields
static Real AphiOpen(const Real x1, const Real x2, const Real x3)
{
  Real R = x1*sin(x2);
  return 2.0/(3.0-alpha)*Bz0*pow(R/R0,-0.5*(alpha-1.0))/pow(1.0+1.0/SQR(tan(x2)*mu),0.625);
}

//--------------------------------------------------------------------------------------
//!\f User source term to set density floor and adjust temperature 
void MySource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &w,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(static)
   for (int j=js; j<=je; ++j){
#pragma simd
     for (int i=is; i<=ie; ++i){
       Real x1   = pmb->pcoord->x1v(i);
       u(IDN,k,j,i) = std::max(u(IDN,k,j,i), rho_floor*pow(x1/R0,-alpha));
  }}}
  return;
}

//--------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
// 

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real x0   = pco->x1v(is);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real rhos = prim(IDN,k,j,is);
//      Real Ts   = (GM/x0)*pmb->ruser_meshblock_data[0](j);
      Real v1s  = prim(IM1,k,j,is);
      Real v2s  = prim(IM2,k,j,is);
      Real v3s  = prim(IM3,k,j,is);
Real Ts = prim(IEN,k,j,is)/prim(IDN,k,j,is);
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1v(is-i); 

        prim(IDN,k,j,is-i) = rhos*pow(x/x0,-alpha);
        prim(IM1,k,j,is-i) = v1s > 0.0 ? 0.0 : v1s;
        prim(IM2,k,j,is-i) = v2s;//*(x/x0);
        prim(IM3,k,j,is-i) = v3s*sqrt(x0/x);
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,is-i) = prim(IDN,k,j,is-i)*Ts*(x0/x);
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    x0   = pco->x1f(is);
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1f(is-i);
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is)*SQR(x0/x);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
      }
    }}

    x0   = pco->x1v(is);
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1v(is-i);
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is)*(x0/x);
      }
    }}
  }

  return;
}

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real x0   = pco->x1v(ie);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real rhoe = prim(IDN,k,j,ie);
//      Real Te   = (GM/x0)*pmb->ruser_meshblock_data[0](j);
      Real v1e  = prim(IM1,k,j,ie);
      Real v2e  = prim(IM2,k,j,ie);
      Real v3e  = prim(IM3,k,j,ie);
Real Te = prim(IEN,k,j,ie)/prim(IDN,k,j,ie);
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1v(ie+i);

        prim(IDN,k,j,ie+i) = rhoe*pow(x/x0,-alpha);
        prim(IM1,k,j,ie+i) = v1e < 0.0 ? 0.0 : v1e;
        prim(IM2,k,j,ie+i) = v2e;
        prim(IM3,k,j,ie+i) = v3e*sqrt(x0/x);
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,ie+i) = prim(IDN,k,j,ie+i)*Te*(x0/x);
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    x0   = pco->x1f(ie+1);
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1f(ie+i+1);
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1))*SQR(x0/x);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
      }
    }}

    x0   = pco->x1v(ie);
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1v(ie+i);
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie)*(x0/x);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//!\f: Userwork in loop: apply external B field and handle user-defined output
//
//void MeshBlock::UserWorkInLoop(ParameterInput *pin)
void MeshBlock::UserWorkInLoop()
{
  Real time = pmy_mesh->time;
  Real dt   = pmy_mesh->dt;

// set temperature profile
#pragma simd
  for (int j=js-NGHOST;j<=je+NGHOST;++j) {
    Real delta = fabs(pcoord->x2v(j)-0.5*PI);
    Real ext = 3.0;
    Real y   = HoR0*(HoRc*ext-1.0)/(1.0-ext*HoR0);
    Real myHoR = std::max(HoR0,(y+HoRc)*delta-y);
    ruser_meshblock_data[0](j) = SQR(std::min(myHoR,HoRc));
  }

// Adjust temperature
  Real mygam = peos->GetGamma();
  Real gam1  = 1.0/(mygam-1.0);

  int dk = NGHOST;
  if (block_size.nx3 == 1)
    dk = 0;

// now adjust temperature
  for (int k=ks-dk; k<=ke+dk; ++k){
#pragma omp for schedule(static)
    for (int j=js-NGHOST; j<=je+NGHOST; ++j){
      Real theta = pcoord->x2v(j);
      if ((pbval->block_bcs[INNER_X2] > 0) && (j<js))
        theta = pcoord->x2v(2*js-j-1);
      if ((pbval->block_bcs[OUTER_X2] > 0) && (j>je))
        theta = pcoord->x2v(2*je-j+1);
      Real sintheta = sin(theta);

#pragma simd
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i){
        Real x1   = pcoord->x1v(i);
        Real R    = x1*sintheta;
        Real OmgK = sqrt(GM/(R*R*R));
        Real odt  = OmgK*dt;
        Real facT = std::max(100.0*exp(-SQR(R-1.0)),1.0);
        Real facR = std::max(std::min(2.0*(Rbuf-R),1.0),0.0);

        // density and velocity adjustment
        Real rho    = phydro->w(IDN,k,j,i);
        Real vr     = phydro->w(IM1,k,j,i);
        Real vtheta = phydro->w(IM2,k,j,i);
        Real vphi   = phydro->w(IM3,k,j,i);
        Real M1     = rho*vr;
        Real M2     = rho*vtheta;//*fmult_v;
        Real M3     = rho*vphi;

        // temperature adjustment
        Real Tw   = phydro->w(IEN,k,j,i)/rho;
        Real myT  = GM/x1*ruser_meshblock_data[0](j);
        Real dT   = (myT-Tw)*(std::min(odt*fcool*facT,1.0));
        Real newT = Tw+dT;
        newT = std::min(std::max(newT, 0.2*myT), 5.0*myT);

        // apply density floor and Alfven limiter
        Real rhonew0 = std::max(rho, rho_floor*pow(x1/R0,-alpha));
        Real b1   = pfield->bcc(IB1,k,j,i);
        Real b2   = pfield->bcc(IB2,k,j,i);
        Real b3   = pfield->bcc(IB3,k,j,i);
        Real bsq  = SQR(b1)+SQR(b2)+SQR(b3);
        Real vA2  = bsq/rhonew0;
        Real enhac= vA2 < SQR(2.0*x1) ? 1.0 : vA2/SQR(2.0*x1);
        Real rhonew = rhonew0*(1.0+(enhac-1.0)*std::min(facR,1.0));

        // update primitive and conserved variables
        phydro->w(IDN,k,j,i) = rhonew;
        phydro->u(IDN,k,j,i) = rhonew;
        phydro->w(IM1,k,j,i) = M1/rhonew;
        phydro->u(IM1,k,j,i) = M1;
        phydro->w(IM2,k,j,i) = M2/rhonew;
        phydro->u(IM2,k,j,i) = M2;
        phydro->w(IM3,k,j,i) = M3/rhonew;
        phydro->u(IM3,k,j,i) = M3;

        Real Ek = 0.5*(SQR(M1)+SQR(M2)+SQR(M3))/rhonew;
        phydro->w(IEN,k,j,i) = rhonew0*newT;
        phydro->u(IEN,k,j,i) = Ek+phydro->w(IEN,k,j,i)*gam1;
        phydro->u(IEN,k,j,i) += 0.5*bsq;
      }

    }
  }

// Apply poloidal field (one shot)
  if (MAGNETIC_FIELDS_ENABLED) {
    if ((time+dt >= taddBp) && (time < taddBp)) {
      AddPoloidalField(this);
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void AddPoloidalField(MeshBlock *pmb)
//  \brief Impose poloidal field in the middle of the simulation
void AddPoloidalField(MeshBlock *pmb)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int nx1 = pmb->block_size.nx1+2*NGHOST+1;
  int nx2 = pmb->block_size.nx2+2*NGHOST+1;
  Hydro *phydro = pmb->phydro;
  Field *pfd = pmb->pfield;

  AthenaArray<Real> A3,area,len,len_p1;
  A3.NewAthenaArray(nx2,nx1);
  area.NewAthenaArray(nx1);
  len.NewAthenaArray(nx1);
  len_p1.NewAthenaArray(nx1);

  int dk = NGHOST;
  if (pmb->block_size.nx3 == 1) dk=0;
  int djm = NGHOST;
  int djp = NGHOST;
  if (pmb->pbval->block_bcs[INNER_X2] > 0) djm = 0;
  if (pmb->pbval->block_bcs[OUTER_X2] > 0) djp = 0;

 if (MAGNETIC_FIELDS_ENABLED)
  {
    Real RBmin = 0.;
    Real Phimin = 2.0/(3.0-alpha)*Bz0*pow(RBmin/R0,1.0-0.5*(alpha-1.0));
#pragma omp for schedule(static)
    for (int j=js-djm; j<=je+djp+1; ++j) {
      Real x2 = pmb->pcoord->x2f(j);
      Real sintheta = sin(x2);
      for (int i=is-NGHOST; i<=ie+NGHOST+1; ++i) {
        Real x1 = pmb->pcoord->x1f(i);
        Real R  = x1*sintheta;
        A3(j,i) = std::max(AphiOpen(x1,x2,0.0)*R-Phimin,0.0)/(fabs(R)+TINY_NUMBER);
      }
    }

    for (int k=ks-dk; k<=ke+dk; ++k) {
#pragma omp for schedule(static)
    for (int j=js-djm; j<=je+djp+1; ++j) {
      pmb->pcoord->Face2Area(k,j,is-NGHOST,ie+NGHOST,area);
      pmb->pcoord->Edge3Length(k,j,is-NGHOST,ie+NGHOST+1,len);
#pragma simd
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        pfd->b.x2f(k,j,i) += -(len(i+1)*A3(j,i+1) - len(i)*A3(j,i))/(area(i)+TINY_NUMBER);
      }
    }}

    for (int k=ks-dk; k<=ke+dk; ++k) {
#pragma omp for schedule(static)
    for (int j=js-djm; j<=je+djp; ++j) {
      pmb->pcoord->Face1Area(k,j,is-NGHOST,ie+NGHOST+1,area);
      pmb->pcoord->Edge3Length(k,j  ,is-NGHOST,ie+NGHOST+1,len);
      pmb->pcoord->Edge3Length(k,j+1,is-NGHOST,ie+NGHOST+1,len_p1);
#pragma simd
      for (int i=is-NGHOST; i<=ie+NGHOST+1; ++i) {
        pfd->b.x1f(k,j,i) += (len_p1(i)*A3(j+1,i) - len(i)*A3(j,i))/(area(i)+TINY_NUMBER);
      }
    }}

    for (int k=ks-dk; k<=ke+dk; k++) {
    for (int j=js-djm; j<=je+djp; j++) {
    for (int i=is-NGHOST; i<=ie+NGHOST; i++) {
      const Real& b1_i   = pfd->b.x1f(k,j,i  );
      const Real& b1_ip1 = pfd->b.x1f(k,j,i+1);
      const Real& b2_j   = pfd->b.x2f(k,j  ,i);
      const Real& b2_jp1 = pfd->b.x2f(k,j+1,i);

      Real& bcc1 = pfd->bcc(IB1,k,j,i);
      Real& bcc2 = pfd->bcc(IB2,k,j,i);
      Real  EB0 = 0.5*(SQR(bcc1)+SQR(bcc2));

      const Real& x1f_i  = pmb->pcoord->x1f(i);
      const Real& x1f_ip = pmb->pcoord->x1f(i+1);
      const Real& x1v_i  = pmb->pcoord->x1v(i);
      const Real& dx1_i  = pmb->pcoord->dx1f(i);
      Real lw=(x1f_ip-x1v_i)/dx1_i;
      Real rw=(x1v_i -x1f_i)/dx1_i;
      bcc1 = lw*b1_i + rw*b1_ip1;

      const Real& x2f_j  = pmb->pcoord->x2f(j);
      const Real& x2f_jp = pmb->pcoord->x2f(j+1);
      const Real& x2v_j  = pmb->pcoord->x2v(j);
      const Real& dx2_j  = pmb->pcoord->dx2f(j);
      lw=(x2f_jp-x2v_j)/dx2_j;
      rw=(x2v_j -x2f_j)/dx2_j;
      bcc2 = lw*b2_j + rw*b2_jp1;
      Real  EB1 = 0.5*(SQR(bcc1)+SQR(bcc2));

      if(NON_BAROTROPIC_EOS)
        phydro->u(IEN,k,j,i) += EB1-EB0;
    }}}
  }

  A3.DeleteAthenaArray();
  area.DeleteAthenaArray();
  len.DeleteAthenaArray();
  len_p1.DeleteAthenaArray();

  return;
}
