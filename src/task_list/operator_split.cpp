//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file operator_split.cpp
//  \brief derived class for time integrator task list.  Can create task lists for one
//  of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "task_list.hpp"
#include "operator_split.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"

//----------------------------------------------------------------------------------------
//  OperatorSplitTaskList constructor

OperatorSplitTaskList::OperatorSplitTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{
  // Set the diffusion flags for operator-split diffusion
  // Note: the do_field_diffusion flag has already taken into account [B field enabled]
  DiffusionDriver *pdiff = pm->pdiff;
  do_field_diffusion = (pdiff->phys_def[OAD]  && pdiff->operator_split_def[OAD])  ||
                       (pdiff->phys_def[HALL] && pdiff->operator_split_def[HALL]);

  do_hydro_diffusion = (pdiff->phys_def[ISO_VIS]  && pdiff->operator_split_def[ISO_VIS])  ||
                       (pdiff->phys_def[ISO_COND] && pdiff->operator_split_def[ISO_COND]) ||
                       (pdiff->phys_def[COOL]     && pdiff->operator_split_def[COOL]);

  // Now assemble list of tasks for each step of time integrator
  {using namespace OperatorSplitTaskNames;

    AddTimeIntegratorTask(START_OSRECV,NONE);

    AddTimeIntegratorTask(CLEAN_FLX,NONE);

    // Main operator split operations
    // calculate hydro/field fluxes due to diffusive processess
    // note that hydro/field diffusion flags will be called before executing
    AddTimeIntegratorTask(DIFFUSE_FLD,START_OSRECV|CLEAN_FLX);

    AddTimeIntegratorTask(DIFFUSE_HYD,START_OSRECV|CLEAN_FLX);

    // Update field variables and set field boundary conditions
    if (do_field_diffusion) { // MHD
      AddTimeIntegratorTask(SEND_FLDFLX,DIFFUSE_FLD);
      AddTimeIntegratorTask(RECV_FLDFLX,DIFFUSE_FLD);
      AddTimeIntegratorTask(INT_FLD, RECV_FLDFLX);
      AddTimeIntegratorTask(SEND_FLD,INT_FLD);
      AddTimeIntegratorTask(RECV_FLD,START_OSRECV);
    }

    // Update hydro variables and set hydro boundary conditions
    if (do_hydro_diffusion || NON_BAROTROPIC_EOS){
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(SEND_HYDFLX,DIFFUSE_HYD|DIFFUSE_FLD);
        AddTimeIntegratorTask(RECV_HYDFLX,DIFFUSE_HYD|DIFFUSE_FLD);
        AddTimeIntegratorTask(INT_HYD, RECV_HYDFLX);
      }
      else
        AddTimeIntegratorTask(INT_HYD, DIFFUSE_HYD|DIFFUSE_FLD);

      AddTimeIntegratorTask(SEND_HYD,INT_HYD);
      AddTimeIntegratorTask(RECV_HYD,START_OSRECV);
    }

    if (do_field_diffusion) {
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG, (SEND_HYD|RECV_HYD|SEND_FLD|RECV_FLD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,(INT_HYD|RECV_HYD|INT_FLD|RECV_FLD));
      }
    } else {  // HYDRO
      if(pm->multilevel==true) { // SMR or AMR
        AddTimeIntegratorTask(PROLONG,(SEND_HYD|RECV_HYD));
        AddTimeIntegratorTask(CON2PRIM,PROLONG);
      } else {
        AddTimeIntegratorTask(CON2PRIM,INT_HYD|RECV_HYD);
      }
    }

    // everything else
    AddTimeIntegratorTask(PHY_BVAL,CON2PRIM);
    AddTimeIntegratorTask(CLEAR_OSRECV,PHY_BVAL);

  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void OperatorSplitTaskList::AddTimeIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace OperatorSplitTaskNames;
  switch((id)) {
    case (START_OSRECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::StartOSReceive);
      break;
    case (CLEAR_OSRECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::ClearOSReceive);
      break;

    case (SEND_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::FluxCorrectSend);
      break;
    case (SEND_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::EMFCorrectSend);
      break;

    case (RECV_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::FluxCorrectReceive);
      break;
    case (RECV_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::EMFCorrectReceive);
      break;

    case (INT_HYD):
      task_list_[ntasks].TaskFunc=
      static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
      (&OperatorSplitTaskList::HydroIntegrate);
      break;
    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::FieldIntegrate);
      break;

    case (SEND_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::HydroSend);
      break;
    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::FieldSend);
      break;

    case (RECV_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::HydroReceive);
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::FieldReceive);
      break;

    case (PROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::Prolongation);
      break;
    case (CON2PRIM):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::Primitives);
      break;
    case (PHY_BVAL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::PhysicalBoundary);
      break;

    case (CLEAN_FLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::CleanFluxes);
      break;
    case (DIFFUSE_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::HydroDiffusion);
      break;
    case (DIFFUSE_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&OperatorSplitTaskList::FieldDiffusion);
      break;


    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus OperatorSplitTaskList::StartOSReceive(MeshBlock *pmb, int step)
{
  Real time = pmb->pmy_mesh->time;
  pmb->pbval->StartReceivingAll(time);
  return TASK_SUCCESS;
}

enum TaskStatus OperatorSplitTaskList::ClearOSReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction step with AMR

enum TaskStatus OperatorSplitTaskList::FluxCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFluxCorrection(FLUX_HYDRO);
  return TASK_SUCCESS;
}

enum TaskStatus OperatorSplitTaskList::EMFCorrectSend(MeshBlock *pmb, int step)
{
  //pmb->pbval->SendEMFCorrection(0);
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

enum TaskStatus OperatorSplitTaskList::FluxCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveFluxCorrection(FLUX_HYDRO) == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus OperatorSplitTaskList::EMFCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus OperatorSplitTaskList::HydroIntegrate(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;
  Mesh *pm=pmb->pmy_mesh;
  Real weight = pm->pdiff->dt_sub[step-1]/pm->dt;

  ph->AddFluxDivergenceToAverage(ph->w,pf->bcc,weight,ph->u);
  return TASK_NEXT;
}

enum TaskStatus OperatorSplitTaskList::FieldIntegrate(MeshBlock *pmb, int step)
{
  Field *pf=pmb->pfield;
  Mesh *pm=pmb->pmy_mesh;
  Real weight = pm->pdiff->dt_sub[step-1]/pm->dt;

  pmb->pfield->CT(weight, pf->b);
  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to clean existing hydro fluxes and emfs
enum TaskStatus OperatorSplitTaskList::CleanFluxes(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  if (do_hydro_diffusion || NON_BAROTROPIC_EOS)
    ph->phdif->ClearHydroFlux(ph->flux);

  if (MAGNETIC_FIELDS_ENABLED)
    if (do_field_diffusion)
      pf->pfdif->ClearEMF(pf->e);

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to calculate diffusion fluxes
enum TaskStatus OperatorSplitTaskList::HydroDiffusion(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;

  if (do_hydro_diffusion)
    ph->phdif->CalcHydroDiffusionFlux(ph->w,ph->u,ph->flux);

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to calculate diffusion EMF
enum TaskStatus OperatorSplitTaskList::FieldDiffusion(MeshBlock *pmb, int step)
{
  Field *pf=pmb->pfield;

  if (do_field_diffusion)
    pf->pfdif->CalcFieldDiffusionEMF(pf->b,pf->bcc,pf->e);

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

enum TaskStatus OperatorSplitTaskList::HydroSend(MeshBlock *pmb, int step)
{
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);

  return TASK_SUCCESS;
}

enum TaskStatus OperatorSplitTaskList::FieldSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);

  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus OperatorSplitTaskList::HydroReceive(MeshBlock *pmb, int step)
{
  bool ret;
  ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);

  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus OperatorSplitTaskList::FieldReceive(MeshBlock *pmb, int step)
{
  bool ret;
  ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b);

  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

//--------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus OperatorSplitTaskList::Prolongation(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  Mesh *pm=pmb->pmy_mesh;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  dt = pm->pdiff->dt_sub[step];

  pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                              pm->time+dt, dt);

  return TASK_SUCCESS;
}

enum TaskStatus OperatorSplitTaskList::Primitives(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if(pmb->pbval->nblevel[1][1][0]!=-1) is-=NGHOST;
  if(pmb->pbval->nblevel[1][1][2]!=-1) ie+=NGHOST;
  if(pmb->pbval->nblevel[1][0][1]!=-1) js-=NGHOST;
  if(pmb->pbval->nblevel[1][2][1]!=-1) je+=NGHOST;
  if(pmb->pbval->nblevel[0][1][1]!=-1) ks-=NGHOST;
  if(pmb->pbval->nblevel[2][1][1]!=-1) ke+=NGHOST;

  pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                    phydro->w, pfield->bcc, pmb->pcoord,
                                    is, ie, js, je, ks, ke);

  return TASK_SUCCESS;
}

enum TaskStatus OperatorSplitTaskList::PhysicalBoundary(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  Mesh *pm=pmb->pmy_mesh;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  dt = pm->pdiff->dt_sub[step];
  pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                 pm->time+dt, dt);
  return TASK_SUCCESS;
}
