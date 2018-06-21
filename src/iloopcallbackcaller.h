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

#ifndef ILOOPCALLBACKCALLER_H
#define ILOOPCALLBACKCALLER_H

#include "contact_interface.h"
#include <vector>

namespace LCM = LIGGGHTS::ContactModels;

// forward declaration in different namespace
namespace LAMMPS_NS {
class LAMMPS;
}

namespace LIGGGHTS {
// forward declaration
class ILoopCallbackCallable;

class ILoopCallbackCaller
{
public:
    /* tightLoopCallback interface */
    ILoopCallbackCaller(LAMMPS_NS::LAMMPS *lmp);

    virtual void register_post_force_callback(ILoopCallbackCallable *caller);
    void unregister_post_force_callback(ILoopCallbackCallable *caller);

    void referenceDeleted();
    void updateRequestedCallbacks();
    void callBeginPass() const;
    void call_post_force_callback(const LCM::SurfacesIntersectData &sidata, const LCM::ForceData &i_forces, const LCM::ForceData &j_forces) const;
    
    void call_add_heat_pp_callback(const int i, const int j, const double heatflux);
    void call_add_heat_pw_callback(const LCM::SurfacesIntersectData &sidata, const double heatflux);
    void callEndPass() const;
private:
    std::vector<ILoopCallbackCallable *> registeredCallbacks_; // contains all registered callbacks
    std::vector<ILoopCallbackCallable *> requestedCallbacks_;  // contains callbacks active for the current timestep

    LAMMPS_NS::LAMMPS *lmp_;
};
}

#endif // ILOOPCALLBACKCALLER_H
