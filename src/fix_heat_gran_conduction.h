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

FixStyle(heat/gran/conduction,FixHeatGranCond)
FixStyle(heat/gran,FixHeatGranCond)

#else

#ifndef LMP_FIX_HEATGRAN_CONDUCTION_H
#define LMP_FIX_HEATGRAN_CONDUCTION_H

#include "fix_heat_gran.h"
#include "physicsheatconduction.h"

namespace LAMMPS_NS {

  class FixHeatGranCond : public FixHeatGran {
  public:
    FixHeatGranCond(class LAMMPS *, int, char **);
    ~FixHeatGranCond();
    virtual void post_create();

    int setmask();
    void init();
    virtual void pre_force(int vflag);
    virtual void post_force(int vflag);

  protected:
    virtual void updatePtrs();
    template <int> void post_force_eval(int);

    int iarg_;

    bool store_contact_data_;
    class FixPropertyAtom* fix_conduction_contact_area_;
    class FixPropertyAtom* fix_n_conduction_contacts_;
    class FixPropertyAtom* fix_wall_heattransfer_coeff_;
    class FixPropertyAtom* fix_wall_temperature_;
    double *conduction_contact_area_;
    double *n_conduction_contacts_;
    double *wall_heattransfer_coeff_;
    double *wall_temp_;

    double temp_max_;

  private:
    PhysicsHeatConduction heatPhy;
  };

}

#endif
#endif

