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

#include "iloopcallbackcaller.h"

#include "error.h"
#include "update.h"
#include "lmptype.h"

#include <algorithm>

using namespace LIGGGHTS;
using LIGGGHTS::ContactModels::SurfacesIntersectData;
using LIGGGHTS::ContactModels::ForceData;

ILoopCallbackCaller::ILoopCallbackCaller(LAMMPS_NS::LAMMPS *lmp) :
    lmp_(lmp)
{

}

void ILoopCallbackCaller::register_post_force_callback(ILoopCallbackCallable *caller)
{
    if (std::find(registeredCallbacks_.begin(), registeredCallbacks_.end(), caller) == registeredCallbacks_.end() )
        registeredCallbacks_.push_back(caller);
    else
        lmp_->error->all(FLERR,"douple registration of postForceCallback\n");
}

void ILoopCallbackCaller::unregister_post_force_callback(ILoopCallbackCallable *caller)
{
    std::vector<ILoopCallbackCallable*>::iterator it = std::find(registeredCallbacks_.begin(), registeredCallbacks_.end(), caller);
    if (it != registeredCallbacks_.end() )
        registeredCallbacks_.erase(it);
    else
        lmp_->error->all(FLERR,"try to unregister not-existing postForceCallback\n");
}

void ILoopCallbackCaller::referenceDeleted()
{
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = requestedCallbacks_.begin(); it != requestedCallbacks_.end(); ++it)
    {
        (*it)->deleteReference();
    }
}

void ILoopCallbackCaller::updateRequestedCallbacks()
{
    requestedCallbacks_.clear();
    const LAMMPS_NS::bigint currentStep = lmp_->update->ntimestep;
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = registeredCallbacks_.begin(); it != registeredCallbacks_.end(); ++it)
    {
        if ((*it)->callRequired(currentStep))
            requestedCallbacks_.push_back((*it));
    }
}

void ILoopCallbackCaller::callBeginPass() const
{
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = requestedCallbacks_.begin(); it != requestedCallbacks_.end(); ++it)
    {
        (*it)->beginPass();
    }
}

void ILoopCallbackCaller::call_post_force_callback(const SurfacesIntersectData & sidata, const ForceData & i_forces, const ForceData & j_forces) const
{
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = requestedCallbacks_.begin(); it != requestedCallbacks_.end(); ++it)
    {
        (*it)->compute_post_force(sidata,i_forces,j_forces);
    }
}

void ILoopCallbackCaller::call_add_heat_pp_callback(const int i, const int j, const double heatflux)
{
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = requestedCallbacks_.begin(); it != requestedCallbacks_.end(); ++it)
    {
        (*it)->add_heat_pp(i,j,heatflux);
    }
}

void ILoopCallbackCaller::call_add_heat_pw_callback(const ContactModels::SurfacesIntersectData &sidata, const double heatflux)
{
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = requestedCallbacks_.begin(); it != requestedCallbacks_.end(); ++it)
    {
        (*it)->add_heat_pw(sidata,heatflux);
    }
}

void ILoopCallbackCaller::callEndPass() const
{
    for (std::vector<ILoopCallbackCallable*>::const_iterator it = requestedCallbacks_.begin(); it != requestedCallbacks_.end(); ++it)
    {
        (*it)->endPass();
    }
}
