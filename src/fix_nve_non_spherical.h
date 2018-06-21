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
    Alexander Podlozhnyuk, DCS Computing GmbH, Linz
    Christoph Kloss, DCS Computing GmbH, Linz

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/nonspherical,FixNVENonSpherical)
FixStyle(nve/convexhull,FixNVENonSpherical)
FixStyle(nve/superquadric,FixNVENonSpherical)

#else

#ifndef LMP_FIX_NVE_NONSPHE_BASE_H
#define LMP_FIX_NVE_NONSPHE_BASE_H

#include "fix_nve.h"
#include "contact_force_corrector.h"
#include "fix_multisphere.h"

namespace LAMMPS_NS {

class FixNVENonSpherical : public FixNVE {
 public:
  FixNVENonSpherical(class LAMMPS *, int, char **);
  virtual ~FixNVENonSpherical();
  virtual void init();
  virtual void initial_integrate(int);
  void initial_integrate_clump(int);
  void initial_integrate_single(int);
  virtual void final_integrate();
  void final_integrate_clump();
  void final_integrate_single();
  void dynamic_euler(double *wbody, double *tbody, double *inertia, double *result);
  void integrate_dynamic_euler(double dt, double *wbody, double *tbody, double *inertia);
  void integrate_quaternion(double *quat, double *wbody);
  int setmask();
  void setup(int);
  void pre_final_integrate();
  void convexHullInit();
  bool concave_;
  void reset_dt();
  void richardson_initial_integrate(double dt, double *angmom, double *omega, double *torque, double *quat, double *inertia,
                                    bool *tflag, double *ex_space, double *ey_space, double *ez_space);
  void richardson_final_integrate(double dt, double *angmom, double *omega, double *torque,  double *quat, double *inertia,
                                    bool *tflag, double *ex_space, double *ey_space, double *ez_space);
  
  void symplectic_initial_integrate(double dt, double *angmom, double *quat, double *omega, double *torque, double *inertia);
  void symplectic_final_integrate(double dtf2, double *angmom, double *quat, double *omega, double *torque, double *inertia);

  void predictor_corrector_initial_integrate(double dtf, double *angmom, double *quat, double *omega, double *torque, double *inertia);
  void predictor_corrector_final_integrate(double dtf, double *angmom, double *quat, double *omega, double *torque, double *inertia);
  void woodem_initial_integrate(double dtf, double *angmom, double *quat, double *omega, double *torque, double *inertia);

  class Multisphere* multisphere_;
  class FixMultisphere *fix_multisphere_;

protected:
  double CAddRhoFluid_; //CAddRhoFluid
  double dtq;

 private:
  int integration_scheme;
  std::list<ContactForceCorrector*> correctors_;
  void modify_body_forces_torques();
  int particle_type;
  
  bool integrate_concave;
  bool integrate_superquadric;
  bool integrate_multisphere;
  bool integrate_convexhull;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/superquadric requires atom style superquadric

Self-explanatory.

*/

