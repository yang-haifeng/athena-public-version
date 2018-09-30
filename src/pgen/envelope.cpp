#include <iostream> // cout endl
#include <cmath>
#include <algorithm> // min

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

static Real gm0, rho0, T0, gamma_gas;
static Real rc=1., ind_r=-1., ind_T=-1.;

void Mesh::InitUserMeshData(ParameterInput *pin){
  gm0 = pin->GetOrAddReal("problem","GM", 1.0);
  rho0 = pin->GetReal("problem","rho0");
  T0 = pin->GetReal("problem", "T0");

  if (NON_BAROTROPIC_EOS)
    gamma_gas = pin->GetReal("hydro","gamma");

  EnrollUserBoundaryFunction(INNER_X1, DiskInnerX1);
  EnrollUserBoundaryFunction(OUTER_X1, DiskOuterX1);

  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin){
  for (int k=ks; k<=ke; k++){
    for (int j=js; j<=je; j++){
      for (int i=is; i<ie; i++){
        Real r = pcoord->x1v(i);
	Real rho = rho0*pow(r/rc, ind_r);
	Real T = T0*pow(r/rc, ind_T);
	phydro->u(IDN,k,j,i) = rho;
	phydro->u(IM1,k,j,i) = 0.;
        phydro->u(IM2,k,j,i) = 0.;
        phydro->u(IM3,k,j,i) = 0.;
	if (NON_BAROTROPIC_EOS) {
          Real press = rho*T;
          phydro->u(IEN,k,j,i) = press/(gamma_gas - 1.0);
       }
      }
    }
  }
}

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        Real r = pco->x1v(is-i);
	Real rho = rho0*pow(r/rc, ind_r);
	Real T = T0*pow(r/rc, ind_T);
	prim(IDN,k,j,is-i) = rho;
        prim(IM1,k,j,is-i) = 0.;
        prim(IM2,k,j,is-i) = 0.;
        prim(IM3,k,j,is-i) = 0.;
	if (NON_BAROTROPIC_EOS) {
          Real press = rho*T;
          prim(IEN,k,j,is-i) = press/(gamma_gas - 1.0);
       }
      }
    }
  }
}

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        Real r = pco->x1v(ie+i);
	Real rho = rho0*pow(r/rc, ind_r);
	Real T = T0*pow(r/rc, ind_T);
	prim(IDN,k,j,ie+i) = rho;
        prim(IM1,k,j,ie+i) = 0.;
        prim(IM2,k,j,ie+i) = 0.;
        prim(IM3,k,j,ie+i) = 0.;
	if (NON_BAROTROPIC_EOS) {
          Real press = rho*T;
          prim(IEN,k,j,ie+i) = press/(gamma_gas - 1.0);
       }
      }
    }
  }
}
