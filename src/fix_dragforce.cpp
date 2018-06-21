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

#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "fix_dragforce.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "region.h"
#include "vector_liggghts.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

double FixDragforce::small_ =  0.00000001;

/* ---------------------------------------------------------------------- */

FixDragforce::FixDragforce(LAMMPS *lmp, int narg, char **arg) :
  FixBaseLiggghts(lmp, narg, arg),
  iarg_(3),
  draglaw_(DRAGLAW_NONE),
  Cd_(0.),
  viscosity_(0.),
  density_(0.)
{
  do_support_respa();

  vectorZeroize3D(U_fluid_);

  if (narg < 6) error->fix_error(FLERR,this,"not enough arguments");

  // args
  bool hasargs = true;
  while (iarg_ < narg && hasargs) {

    hasargs = false;

    if (strcmp(arg[iarg_],"Schiller_Naumann") == 0) {
      draglaw_ = DRAGLAW_SCHILLER_NAUMANN;
      if (iarg_+5 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'Schiller_Naumann'");
      if (strcmp(arg[iarg_+1],"viscosity"))
        error->fix_error(FLERR,this,"expecting 'viscosity' after 'Schiller_Naumann'");
      viscosity_ = atof(arg[iarg_+2]);
      if (viscosity_ <= 0.)
        error->fix_error(FLERR,this,"'viscosity' > 0 required");
      if (strcmp(arg[iarg_+3],"density"))
        error->fix_error(FLERR,this,"expecting 'density' after 'viscosity'");
      density_ = atof(arg[iarg_+4]);
      if (density_ <= 0.)
        error->fix_error(FLERR,this,"'density' > 0 required");
      iarg_ += 5;
      hasargs = true;
    } else if (strcmp(arg[iarg_],"const_Cd") == 0) {
      draglaw_ = DRAGLAW_CD_CONST;
      if (iarg_+7 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'const_Cd'");
      if (strcmp(arg[iarg_+1],"Cd"))
        error->fix_error(FLERR,this,"expecting 'Cd' after 'const_Cd'");
      Cd_ = atof(arg[iarg_+2]);
      if (Cd_ <= 0.)
        error->fix_error(FLERR,this,"'Cd' > 0 required");
      if (strcmp(arg[iarg_+3],"viscosity"))
        error->fix_error(FLERR,this,"expecting 'viscosity' after 'Cd'");
      viscosity_ = atof(arg[iarg_+4]);
      if (viscosity_ <= 0.)
        error->fix_error(FLERR,this,"'viscosity' > 0 required");
      if (strcmp(arg[iarg_+5],"density"))
        error->fix_error(FLERR,this,"expecting 'density' after 'viscosity'");
      density_ = atof(arg[iarg_+6]);
      if (density_ <= 0.)
        error->fix_error(FLERR,this,"'density' > 0 required");
      iarg_ += 7;
      hasargs = true;
    } else if (strcmp(arg[iarg_],"U_fluid") == 0) {
      if (iarg_+4 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'U_fluid'");
      U_fluid_[0] = atof(arg[iarg_+1]);
      U_fluid_[1] = atof(arg[iarg_+2]);
      U_fluid_[2] = atof(arg[iarg_+3]);
      iarg_ += 4;
      hasargs = true;
    } else if (strcmp(arg[iarg_],"region") == 0) {
      process_region(arg[iarg_+1]);
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(style,"dragforce") == 0) {
      char *errmsg = new char[strlen(arg[iarg_])+50];
      sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
      error->fix_error(FLERR,this,errmsg);
      delete []errmsg;
    }
  }

  if(DRAGLAW_NONE == draglaw_)
    error->fix_error(FLERR,this,"have to specify either 'Schiller_Naumann' or 'const_Cd'");

  do_support_multisphere();
}

/* ---------------------------------------------------------------------- */

FixDragforce::~FixDragforce()
{
}

/* ---------------------------------------------------------------------- */

int FixDragforce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDragforce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDragforce::post_force(int vflag)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double uRel[3],uMag, Cd=0.0, A, drag;
  double invvisc = 1./viscosity_;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      if (region_ && !region_->match(x[i][0],x[i][1],x[i][2]))
          continue;

      vectorSubtract3D(v[i],U_fluid_,uRel);
      uMag = vectorMag3D(uRel);
      if(uMag > small_)
      {
          double cg = force->cg(type[i]);
          double invcg = 1./cg;
          if(DRAGLAW_SCHILLER_NAUMANN == draglaw_)
          {
              double Re = uMag*(2.*radius[i]*invcg)*invvisc;
              Cd = std::max(24./Re*(1.+0.15*pow(Re,0.687)),0.44);
          }
          else if(DRAGLAW_CD_CONST == draglaw_)
          {
              Cd = Cd_;
          }

          A = (radius[i]*invcg)*(radius[i]*invcg)*M_PI;
          drag = density_*uMag*0.5*A*Cd*cg*cg*cg;

          f[i][0] -= drag*uRel[0];
          f[i][1] -= drag*uRel[1];
          f[i][2] -= drag*uRel[2];
          
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixDragforce::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa_-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDragforce::min_post_force(int vflag)
{
  post_force(vflag);
}
