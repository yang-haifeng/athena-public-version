#ifndef OPERATOR_SPLIT_TASK_LIST_HPP
#define OPERATOR_SPLIT_TASK_LIST_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file operatorsplit_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

#include <stdint.h>

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//----------------------------------------------------------------------------------------
//! \class OperatorSplitTaskList
//  \brief data and function definitions for OperatorSplitTaskList derived class

class OperatorSplitTaskList : public TaskList {
public:
    OperatorSplitTaskList(ParameterInput *pin, Mesh *pm);
    ~OperatorSplitTaskList() {};

    // data
    Real diff_dt_hydro, diff_dt_field;
    bool super_timestepping;
    bool do_field_diffusion, do_hydro_diffusion;

    // functions
    void AddTimeIntegratorTask(uint64_t id, uint64_t dep);

    void EstimateNstage(Mesh *pm);

    enum TaskStatus StartOSReceive(MeshBlock *pmb, int step);
    enum TaskStatus ClearOSReceive(MeshBlock *pmb, int step);

    enum TaskStatus FluxCorrectSend(MeshBlock *pmb, int step);
    enum TaskStatus EMFCorrectSend(MeshBlock *pmb, int step);

    enum TaskStatus FluxCorrectReceive(MeshBlock *pmb, int step);
    enum TaskStatus EMFCorrectReceive(MeshBlock *pmb, int step);

    enum TaskStatus HydroIntegrate(MeshBlock *pmb, int step);
    enum TaskStatus FieldIntegrate(MeshBlock *pmb, int step);

    enum TaskStatus CleanFluxes(MeshBlock *pmb, int step);
    enum TaskStatus HydroDiffusion(MeshBlock *pmb, int step);
    enum TaskStatus FieldDiffusion(MeshBlock *pmb, int step);

    enum TaskStatus HydroSend(MeshBlock *pmb, int step);
    enum TaskStatus FieldSend(MeshBlock *pmb, int step);

    enum TaskStatus HydroReceive(MeshBlock *pmb, int step);
    enum TaskStatus FieldReceive(MeshBlock *pmb, int step);

    enum TaskStatus Prolongation(MeshBlock *pmb, int step);
    enum TaskStatus Primitives(MeshBlock *pmb, int step);
    enum TaskStatus PhysicalBoundary(MeshBlock *pmb, int step);
};

//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.

namespace OperatorSplitTaskNames {
    const uint64_t NONE=0;
    const uint64_t START_OSRECV=1LL<<0;
    const uint64_t CLEAR_OSRECV=1LL<<1;

    const uint64_t SEND_HYDFLX=1LL<<2;
    const uint64_t SEND_FLDFLX=1LL<<3;

    const uint64_t RECV_HYDFLX=1LL<<4;
    const uint64_t RECV_FLDFLX=1LL<<5;

    const uint64_t INT_HYD=1LL<<6;
    const uint64_t INT_FLD=1LL<<7;

    const uint64_t SEND_HYD=1LL<<8;
    const uint64_t SEND_FLD=1LL<<9;

    const uint64_t RECV_HYD=1LL<<10;
    const uint64_t RECV_FLD=1LL<<11;

    const uint64_t PROLONG =1LL<<12;
    const uint64_t CON2PRIM=1LL<<13;
    const uint64_t PHY_BVAL=1LL<<14;

    const uint64_t CLEAN_FLX=1LL<<15;
    const uint64_t DIFFUSE_HYD=1LL<<16;
    const uint64_t DIFFUSE_FLD=1LL<<17;
};
#endif // OPERATOR_SPLIT_TASK_LIST_HPP
