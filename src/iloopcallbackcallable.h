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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Andreas Aigner (DCS Computing GmbH, Linz)

    Copyright 2018-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef ILOOPCALLBACKCALLABLE_H
#define ILOOPCALLBACKCALLABLE_H

#include "contact_interface.h"
#include "lmptype.h"

namespace LCM = LIGGGHTS::ContactModels;

namespace LIGGGHTS{
class ILoopCallbackCallable
{
public:
    /* tightLoopCallback interface */
    virtual void deleteReference() {}
    virtual void beginPass() {}
    virtual void compute_post_force(const LCM::SurfacesIntersectData & sidata, const LCM::ForceData & i_forces, const LCM::ForceData & j_forces) {}
    virtual void add_heat_pp(const int i, const int j, const double heatflux) {}
    virtual void add_heat_pw(const LCM::SurfacesIntersectData & sidata, const double heatflux) {}
    virtual void endPass() {}
    virtual void registerNextCall(LAMMPS_NS::bigint /*step*/) {}
    virtual bool callRequired(LAMMPS_NS::bigint /*cstep*/) { return false; }
};
}

#endif // ILOOPCALLBACKCALLABLE_H
