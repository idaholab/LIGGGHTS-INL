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

#ifdef FIX_CLASS

FixStyle(couple/cfd/force,FixCfdCouplingForce)
FixStyle(couple/cfd/force/implicit,FixCfdCouplingForce)

#else

#ifndef LMP_FIX_CFD_COUPLING_FORCE_H
#define LMP_FIX_CFD_COUPLING_FORCE_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

enum cfdCoupleType {SCALAR,VECTOR,VECTOR2D,QUATERNION,SCALARMULTISPHERE,VECTORMULTISPHERE,SCALARGLOB,VECTORGLOB,MATRIXGLOB};

// starts with property atoms, later <= quaternion is used to check this, so don't mess around

struct coupleParams
{
//coupleParams shall have fix property pointer, bool push, bool pull, type cfdCoupleType
  class Fix* fix_property_;
  bool push;
  bool pull;
  cfdCoupleType type;
};

class FixCfdCouplingForce : public Fix  {
 public:
  FixCfdCouplingForce(class LAMMPS *, int, char **);
  ~FixCfdCouplingForce();
  void post_create();
  void pre_delete(bool unfixflag);
  void registerProp(std::string propName, bool push, bool pull, cfdCoupleType type);

  double getCAddRhoFluid() { return CAddRhoFluid_; };

  int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void post_force(int);
  double compute_vector(int n);

 protected:

  int iarg;

// coupleList shall have fix property pointer, bool push, bool pull, bool vector
  std::map<std::string, coupleParams> coupleList_;

  double CAddRhoFluid_; //Cadd*rhofluid

  class FixPropertyAtom* fix_dragforce_; //explicit part!
  class FixPropertyAtom* fix_dragforce_implicit_;
  class FixPropertyAtom* fix_dragforce_total_; //explicit+implicit
  double dragforce_total[3];
  class FixPropertyAtom* fix_hdtorque_;
  class FixPropertyAtom* fix_hdtorque_implicit_;
  double hdtorque_total[3];

  class FixCfdCoupling* fix_coupling_;

  bool use_torque_;
  bool dragforce_implicit_;
};

}

#endif
#endif
