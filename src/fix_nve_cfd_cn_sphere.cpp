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

    Contributing author and copyright for this file:
    Josef Kerbl (DCS Computing GmbH, Linz)

    Copyright 2018-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_cfd_cn_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h"
#include "fix_cfd_coupling_force.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVECfdCnSphere::FixNVECfdCnSphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg),
//    FixCfdCouplingForce(lmp,narg,arg),
    CNalpha_(0.5),
    CAddRhoFluid_(0.),
    fix_Ksl_(0),
    fix_uf_(0)
//    use_implicitAnisotropic_(false),
//    fix_KslExtra_(0)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/cfd_cn/sphere command");

  // process extra keywords

  extra = NONE;

  int iarg = 3;
  while (iarg < narg) {
      if (strcmp(arg[iarg],"update") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/cfd_cn/sphere command");
          iarg++;
          if (strcmp(arg[iarg],"dipole") == 0)
          {
              extra = DIPOLE;
              iarg++;
          }
      }
      else if(strcmp(arg[iarg],"CrankNicolson") == 0)
      {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CrankNicolson'");
            iarg++;
            CNalpha_ = atof(arg[iarg]);
            if(CNalpha_<0 || CNalpha_>1)
                error->fix_error(FLERR,this,"incorrect choice for 'CrankNicolson': setting CNalpha_<0 or CNalpha_>1 is not appropriate");

            if (comm->me == 0 && screen) fprintf(screen,"nve_cfd_cn_sphere will use Crank-Nicolson scheme with %f\n", CNalpha_);
            iarg++;
        }

//      else if(strcmp(arg[iarg],"implicitAnisotropic") == 0)
//      {
//          if(narg < iarg+2)
//              error->fix_error(FLERR,this,"not enough arguments for 'implicitAnisotropic'");
//          iarg++;
//          if(strcmp(arg[iarg],"yes") == 0)
//              use_implicitAnisotropic_ = true;
//          else if(strcmp(arg[iarg],"no") == 0)
//              use_implicitAnisotropic_ = false;
//          else
//              error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'implicitAnisotropic'");
//          iarg++;
//      }
      else error->all(FLERR,"Illegal fix nve/cfd_cn/sphere command");
  }

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nve/cfd_cn/sphere requires atom style sphere");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix nve/cfd_cn/sphere requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

void FixNVECfdCnSphere::init()
{
  FixNVE::init();

  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (domain->dimension == 2)
      error->one(FLERR,"Fix nve/cfd_cn/sphere does not allow 2D simulations");

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/cfd_cn/sphere requires extended particles");

//  class FixCfdCouplingForce *fix_coupling_force_ = static_cast<FixCfdCouplingForce*>(modify->find_fix_style("couple/cfd/force",0));
  class FixCfdCouplingForce *fix_coupling_force_ = (FixCfdCouplingForce*)modify->find_fix_style("couple/cfd/force",0);
  if (fix_coupling_force_)
  {
      CAddRhoFluid_ = fix_coupling_force_->getCAddRhoFluid();
      if (comm->me == 0 && screen) fprintf(screen,"CAddRhoFluid is %e\n",CAddRhoFluid_);
  }

  fix_Ksl_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("Ksl","property/atom","scalar",0,0,style,false));
  fix_uf_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("uf","property/atom","vector",0,0,style,false));
  fix_dragforce_implicit_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dragforce_implicit","property/atom","vector",0,0,style,false));

  if (!fix_Ksl_)
  {
        const char* fixarg[9];
        fixarg[0]="Ksl";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="Ksl";
        fixarg[4]="scalar";
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_Ksl_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  if (!fix_uf_)
  {
        const char* fixarg[11];
        fixarg[0]="uf";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="uf";
        fixarg[4]="vector";
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_uf_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
  }

  if (!fix_dragforce_implicit_)
  {
        const char* fixarg[11];
        fixarg[0]="dragforce_implicit";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dragforce_implicit";
        fixarg[4]="vector";
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_dragforce_implicit_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

void FixNVECfdCnSphere::initial_integrate(int vflag)
{
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *density  = atom->density;
    double *Ksl = fix_Ksl_->vector_atom;
    double **uf = fix_uf_->array_atom;
    //double **dragforce_implicit = fix_dragforce_implicit_->array_atom;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    double frc[3],KslCurr[3];

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            //Implicit 1/2 step of implicit velocity update
            KslCurr[0] = Ksl[i]; KslCurr[1] = Ksl[i]; KslCurr[2] = Ksl[i];
            implicitVelocityUpdate( dtf,rmass[i], density[i], v[i], f[i], KslCurr, uf[i], frc);

            x[i][0] += dtv  * v[i][0];
            x[i][1] += dtv  * v[i][1];
            x[i][2] += dtv  * v[i][2];
        }
    }

    rotationUpdate();
    dipoleUpdate();
}

/* ---------------------------------------------------------------------- */

void FixNVECfdCnSphere::final_integrate()
{
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *density  = atom->density;
    double *Ksl = fix_Ksl_->vector_atom;
    double **uf = fix_uf_->array_atom;
    double **dragforce_implicit = fix_dragforce_implicit_->array_atom;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    double frc[3],KslCurr[3];
    
    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            //Implicit 1/2 step of implicit velocity update, frc is also half force!
            KslCurr[0] = Ksl[i]; KslCurr[1] = Ksl[i]; KslCurr[2] = Ksl[i];
            implicitVelocityUpdate( dtf,rmass[i], density[i], v[i], f[i], KslCurr, uf[i], frc);
            for (int j = 0; j < 3; j++)
            {
                dragforce_implicit[i][j] = frc[j];
            }
        }
    }

    rotationUpdate();
}

/* --------------------------------------------------------------------- */
void FixNVECfdCnSphere::implicitVelocityUpdate
            (
                double dtf, double rmass, double particleDensity,
                double *v, double *f, double *Ksl, double *uf,
                double *frc
            )
{
      double vN32[3], deltaU, dtfm, KslMDeltaT;

      dtfm = dtf / ( rmass* (1+CAddRhoFluid_/particleDensity) );

      for(int dirI=0;dirI<3;dirI++)
      {
            KslMDeltaT = Ksl[dirI] * dtfm;

            //calculate new velocity
            vN32[dirI] = (  v[dirI] + f[dirI] * dtfm + KslMDeltaT * ( uf[dirI] - (1.0-CNalpha_)*v[dirI] )   )
                         /
                         (1.0+KslMDeltaT*CNalpha_);

            //calculate velocity difference and force
            deltaU    =  uf[dirI] - ( (1.0-CNalpha_)*v[dirI] + CNalpha_ *vN32[dirI]  );
            frc[dirI]  = Ksl[dirI] * deltaU; 

            //update the particle velocity
            v[dirI] = vN32[dirI];  //update velocity for a half step!
      }
}

/* ---------------------------------------------------------------------- */
void FixNVECfdCnSphere::dipoleUpdate()
{
  if(extra != DIPOLE)
    return;

  double **omega   = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double g[3];
  double msq,scale;
  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length
  double **mu = atom->mu;
  for (int i = 0; i < nlocal; i++)
  {
      if (mask[i] & groupbit)
      {
          if (mu[i][3] > 0.0)
          {
                g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
                g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
                g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
                msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
                scale = mu[i][3]/sqrt(msq);
                mu[i][0] = g[0]*scale;
                mu[i][1] = g[1]*scale;
                mu[i][2] = g[2]*scale;
          }
      }
  }
}

/* ---------------------------------------------------------------------- */
void FixNVECfdCnSphere::rotationUpdate()
{
  // This is okay only for SPHERICAL particles
  // update 1/2 step for omega
  double **omega   = atom->omega;
  double **torque  = atom->torque;
  double *rmass    = atom->rmass;
  double *radius   = atom->radius;
  double *density   = atom->density;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dtirotate;
  double dtfrotate  = dtf / INERTIA;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
    {
      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]*(1+CAddRhoFluid_/density[i]));
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }
}
