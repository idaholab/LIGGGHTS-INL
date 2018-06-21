/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_set_vel.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define FIX_SET_VEL_V 1
#define FIX_SET_VEL_W 2

/* ---------------------------------------------------------------------- */

FixSetVel::FixSetVel(LAMMPS *lmp, int narg, char **arg) :
    FixBaseLiggghts(lmp, narg, arg),
    dim_(-1),
    type_(0),
    vel_(0.),
    constant_(true),
    varid_(0),
    ivar_(-1),
    itype_(-1),
    fix_ms_(NULL),
    multisphere_(NULL)
{
    int iarg = 3;
    do_support_multisphere();

    // parse args
    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg],"dim") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'dim'");
            if(strcmp(arg[iarg+1],"x") == 0)
              dim_ = 0;
            else if(strcmp(arg[iarg+1],"y") == 0)
              dim_ = 1;
            else if(strcmp(arg[iarg+1],"z") == 0)
              dim_ = 2;
            else
              error->fix_error(FLERR,this,"expecting 'x', 'y' or 'z' after 'dim'");
            iarg += 2;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"vel") == 0)
        {
              type_ = FIX_SET_VEL_V;
              if(narg < iarg+2)
                  error->fix_error(FLERR,this,"not enough arguments for 'vel'");

              if (strstr(arg[iarg+1],"v_") == arg[iarg+1])
              {
                  int n = strlen(&arg[iarg+1][2]) + 1;
                  varid_ = new char[n];
                  strcpy(varid_,&arg[iarg+1][2]);
                  constant_ = false;
              }
              else
              {
                  vel_ = atof(arg[iarg+1]);
                  constant_ = true;
              }
              iarg += 2;
              hasargs = true;
        }
        else if(strcmp(arg[iarg],"omega") == 0)
        {
              type_ = FIX_SET_VEL_W;
              if(narg < iarg+2)
                  error->fix_error(FLERR,this,"not enough arguments for 'vel'");

              if (strstr(arg[iarg+1],"v_") == arg[iarg+1])
              {
                  int n = strlen(&arg[iarg+1][2]) + 1;
                  varid_ = new char[n];
                  strcpy(varid_,&arg[iarg+1][2]);
                  constant_ = false;
              }
              else
              {
                  vel_ = atof(arg[iarg+1]);
                  constant_ = true;
              }
              iarg += 2;
              hasargs = true;
        }
        else if(strcmp(arg[iarg],"region") == 0)
        {
              if(narg < iarg+2)
                  error->fix_error(FLERR,this,"not enough arguments for 'region'");
              process_region(arg[iarg+1]);
              iarg += 2;
              hasargs = true;
        }
        else if(strcmp(arg[iarg],"type") == 0)
        {
              if(narg < iarg+2)
                  error->fix_error(FLERR,this,"not enough arguments for 'type'");
              itype_ = atoi(arg[iarg+1]);
              iarg += 2;
              hasargs = true;
        }
        else error->fix_error(FLERR,this,"unknown keyword");
    }

    //time_integrate = 1;
    scalar_flag = 1;
    global_freq = 1;
    extscalar = 1;
    force_flag = 0;
    // error checks

    if (!type_)
        error->fix_error(FLERR, this, "Could not find keyword vel or omega");
    if (dim_ == -1)
        error->fix_error(FLERR, this, "definition of 'dim' required");
}

/* ---------------------------------------------------------------------- */

FixSetVel::~FixSetVel()
{
    if(varid_)
        delete [] varid_;
}

/* ---------------------------------------------------------------------- */

int FixSetVel::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetVel::init()
{
    FixBaseLiggghts::init();

    // check variables

    if (varid_)
    {
        ivar_ = input->variable->find(varid_);
        if (ivar_ < 0)
            error->fix_error(FLERR,this,"Variable name does not exist");
        if (!input->variable->equalstyle(ivar_))
            error->fix_error(FLERR,this,"Variable must be of style 'equal'");
    }

    fix_ms_ = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
    if(fix_ms_)
    {
        multisphere_ = &fix_ms_->data();
        if (iregion_ != -1)
            error->fix_error(FLERR, this, "Fix set/vel does not support region for multispheres, only keyword type is supported");
    }

    int myindex = modify->my_index(this);
    int intindex = modify->index_first_fix_with_function(INITIAL_INTEGRATE,true);

    if (intindex < myindex && intindex >= 0)
        error->fix_error(FLERR,modify->fix[intindex],"has to be located after fix set/vel");

    if (fix_ms_)
    {
        int intmsindex = modify->index_first_fix_of_style("multisphere");
        if (intmsindex < myindex)
            error->fix_error(FLERR,modify->fix[intmsindex],"has to be located after fix set/vel");
    }
}

/* ---------------------------------------------------------------------- */

void FixSetVel::setup(int vflag)
{
    FixBaseLiggghts::setup(vflag);

    if (strstr(update->integrate_style,"verlet"))
        initial_integrate(vflag);
    else
        error->fix_error(FLERR,this,"respa not implemented");
}

/* ---------------------------------------------------------------------- */

void FixSetVel::min_setup(int vflag)
{
    initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetVel::initial_integrate(int vflag)
{
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    int * const mask = atom->mask;
    int * const type = atom->type;
    int nlocal = atom->nlocal;

    // set force to zero and set velocity

    if (!constant_ || (iregion_ >= 0 && domain->regions[iregion_]->has_varshape()))
    {
        modify->clearstep_compute();
        if (!constant_)
            vel_ = input->variable->compute_equal(ivar_);
        if (iregion_ >= 0 && domain->regions[iregion_]->has_varshape())
            domain->regions[iregion_]->update_region();
        modify->addstep_compute(update->ntimestep + 1);
        
    }

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (iregion_ >= 0 &&
                !domain->regions[iregion_]->match(x[i][0],x[i][1],x[i][2]))
                continue;
            if (itype_ > 0 && type[i] != itype_)
                continue;
            switch(type_)
            {
            case FIX_SET_VEL_V:
                f[i][dim_] = 0.;
                v[i][dim_] = vel_;
                break;
            case FIX_SET_VEL_W:
                torque[i][dim_] = 0.;
                omega[i][dim_] = vel_;
                break;
            }
        }
    }

    if(fix_ms_)
    {
        const int nbody_local = multisphere_->n_body();

        for(int ibody_local = 0; ibody_local < nbody_local ; ibody_local++)
        {
            if (itype_ > 0 && multisphere_->atomtype(ibody_local) != itype_)
                continue;
            switch(type_)
            {
            case FIX_SET_VEL_V:
                double v[3];
                multisphere_->vcm(v, ibody_local);
                v[dim_] = vel_;
                multisphere_->set_v_body(ibody_local, v);
                break;
            case FIX_SET_VEL_W:
                double w[3];
                multisphere_->omega(w, ibody_local);
                w[dim_] = vel_;
                multisphere_->set_angmom_via_omega_body(ibody_local, w);
                break;
            }
        }
    }
}
/* ---------------------------------------------------------------------- */

void FixSetVel::final_integrate()
{
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    int * const mask = atom->mask;
    int * const type = atom->type;
    int nlocal = atom->nlocal;

    // set force to zero and set velocity
    // does not really impact on position, just on velocity

    foriginal = 0.0;
    force_flag = 0;

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (iregion_ >= 0 &&
                !domain->regions[iregion_]->match(x[i][0],x[i][1],x[i][2]))
                continue;
            if (itype_ > 0 && type[i] != itype_)
                continue;
            switch(type_)
            {
            case FIX_SET_VEL_V:
                foriginal += f[i][dim_];
                f[i][dim_] = 0.;
                v[i][dim_] = vel_;
                break;
            case FIX_SET_VEL_W:
                foriginal += torque[i][dim_];
                torque[i][dim_] = 0.;
                omega[i][dim_] = vel_;
                break;
            }
        }
    }

    if(fix_ms_)
    {
        const int nbody_local = multisphere_->n_body();

        for(int ibody_local = 0; ibody_local < nbody_local ; ibody_local++)
        {
            if (itype_ > 0 && multisphere_->atomtype(ibody_local) != itype_)
                continue;
            switch(type_)
            {
            case FIX_SET_VEL_V:
                double v[3];
                multisphere_->vcm(v, ibody_local);
                v[dim_] = vel_;
                multisphere_->set_v_body(ibody_local, v);
                break;
            case FIX_SET_VEL_W:
                double w[3];
                multisphere_->omega(w, ibody_local);
                w[dim_] = vel_;
                multisphere_->set_angmom_via_omega_body(ibody_local, w);
                break;
            }
        }
    }
}

double FixSetVel::compute_scalar()
{
    // only sum across procs one time
    if (force_flag == 0)
    {
        MPI_Sum_Scalar(foriginal,world);
        force_flag = 1;
    }
    return foriginal;
}
