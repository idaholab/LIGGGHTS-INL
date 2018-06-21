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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_SJKR,sjkr,1)
#else

#ifndef COHESION_MODEL_SJKR_H_
#define COHESION_MODEL_SJKR_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_SJKR> : public CohesionModelBase {
  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * c) :
        CohesionModelBase(lmp, hsetup, c),
        cohEnergyDens(NULL)
    {
        
    }

    void registerSettings(Settings& settings) 
    {
        settings.registerOnOff("tangential_reduce",tangentialReduce_,false);
    }

    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {
        registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
        registry.connect("cohEnergyDens", cohEnergyDens,"cohesion_model sjkr");
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      //r is the distance between the sphere's centers
      const double r = sidata.r;
      const double ri = sidata.radi;
      const double rj = sidata.radj;

      double Acont;
      if(sidata.is_wall)
        Acont = (ri*ri-r*r)*M_PI*sidata.area_ratio; //contact area sphere-wall
      else
        Acont = - M_PI/4 * ( (r-ri-rj)*(r+ri-rj)*(r-ri+rj)*(r+ri+rj) )/(r*r); //contact area of the two spheres
      const double Fn_coh = -cohEnergyDens[sidata.itype][sidata.jtype]*Acont;
      if(tangentialReduce_) sidata.Fn += Fn_coh; 

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
      double torque_i[3] = {0., 0., 0.};
      double Fn_i[3] = { Fn_coh * sidata.en[0], Fn_coh * sidata.en[1], Fn_coh * sidata.en[2]};
      if (sidata.is_non_spherical)
      {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
      }
      #endif

      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn_coh * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];

        #ifdef NONSPHERICAL_ACTIVE_FLAG
        if(sidata.is_non_spherical)
        {
            //for non-spherical particles normal force can produce torque!
            i_forces.delta_torque[0] += sidata.area_ratio*torque_i[0];
            i_forces.delta_torque[1] += sidata.area_ratio*torque_i[1];
            i_forces.delta_torque[2] += sidata.area_ratio*torque_i[2];
        }
        #endif
      } else {
        const double fx = Fn_coh * sidata.en[0];
        const double fy = Fn_coh * sidata.en[1];
        const double fz = Fn_coh * sidata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;

        #ifdef NONSPHERICAL_ACTIVE_FLAG
        if(sidata.is_non_spherical)
        {
            //for non-spherical particles normal force can produce torque!
            double xcj[3], torque_j[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
            vectorCross3D(Fn_i, xcj, torque_j);

            i_forces.delta_torque[0] += torque_i[0];
            i_forces.delta_torque[1] += torque_i[1];
            i_forces.delta_torque[2] += torque_i[2];

            j_forces.delta_torque[0] += torque_j[0];
            j_forces.delta_torque[1] += torque_j[1];
            j_forces.delta_torque[2] += torque_j[2];
        }
        #endif
      }
    }

    void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesClose(SurfacesCloseData& scdata, ForceData&, ForceData&)
    {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;

    }

  private:
    double ** cohEnergyDens;
    bool tangentialReduce_;
  };
}
}

#endif // COHESION_MODEL_SJKR_H_
#endif
