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

#ifdef FIX_CLASS

FixStyle(nve/cfd_cn/sphere,FixNVECfdCnSphere)

#else

#ifndef LMP_FIX_NVE_CFD_CN_SPHERE_H
#define LMP_FIX_NVE_CFD_CN_SPHERE_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVECfdCnSphere : public FixNVE {
   public:
      FixNVECfdCnSphere(class LAMMPS *, int, char **);
      virtual ~FixNVECfdCnSphere() {}
      virtual void init();
      virtual void initial_integrate(int);
      virtual void final_integrate();

   protected:
      void rotationUpdate();

      void dipoleUpdate();

      inline void implicitVelocityUpdate( double dtf, double rmass, double particleDensity,
                    double *v, double *f, double* Ksl, double *uf, double *frc );

      int extra;

      double CNalpha_;  //Crank-Nicolson blending factor
      double CAddRhoFluid_; //CAddRhoFluid
      class FixPropertyAtom* fix_Ksl_;
      class FixPropertyAtom* fix_uf_;
      class FixPropertyAtom* fix_dragforce_implicit_;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/sphere requires atom style sphere

Self-explanatory.

E: Fix nve/sphere requires atom attribute mu

An atom style with this attribute is needed.

E: Fix nve/sphere requires extended particles

This fix can only be used for particles of a finite size.

*/
