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

    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SURFACE_MODEL
SURFACE_MODEL(SURFACE_SUPERQUADRIC_ORTHOGONAL,superquadric/orthogonal,6)
#else
#ifndef SURFACE_MODEL_SUPERQUADRIC_ORTHOGONAL_H_
#define SURFACE_MODEL_SUPERQUADRIC_ORTHOGONAL_H_

#include "surface_model_base.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG

#include "contact_models.h"
#include <cmath>
#include <algorithm>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "math_extra_liggghts_superquadric.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class SurfaceModel<SURFACE_SUPERQUADRIC_ORTHOGONAL> : public SurfaceModelBase
  {
    int inequality_start_offset;
    int particles_were_in_contact_offset;
    int contact_point_offset;
    int alpha1_offset;
    int alpha2_offset;
    int contact_point_i_offset;
    int contact_point_j_offset;
    int reff_offset;
    Superquadric particle_i;
    Superquadric particle_j;
    enum {SURFACES_FAR, SURFACES_CLOSE, SURFACES_INTERSECT};
  public:
    SurfaceModel(LAMMPS * lmp, IContactHistorySetup* hsetup,class ContactModelBase *cmb) :
        SurfaceModelBase(lmp, hsetup, cmb)
    {
      if(!atom->superquadric_flag)
        error->one(FLERR,"Applying surface model superquadric to a non-superquadric particle!");
      inequality_start_offset = hsetup->add_history_value("inequality_obb","0");
      particles_were_in_contact_offset = hsetup->add_history_value("particles_were_in_contact","0");
      contact_point_offset = hsetup->add_history_value("cpx", "0");
      hsetup->add_history_value("cpy", "0");
      hsetup->add_history_value("cpz", "0");
      alpha1_offset = hsetup->add_history_value("a1", "0");
      alpha2_offset = hsetup->add_history_value("a2", "0");
      contact_point_i_offset = hsetup->add_history_value("с1x", "0");
      hsetup->add_history_value("с1y", "0");
      hsetup->add_history_value("с1z", "0");
      contact_point_j_offset = hsetup->add_history_value("с2x", "0");
      hsetup->add_history_value("с2y", "0");
      hsetup->add_history_value("с2z", "0");
    }

    void registerSettings(Settings& settings) {
      settings.registerDoubleSetting("curvatureLimitFactor",curvatureLimitFactor, 0.0);
      settings.registerYesNo("meanCurvature", meanCurvature, false);
      settings.registerYesNo("gaussianCurvature", gaussianCurvature, false);
      if(curvatureLimitFactor < 0.0)
        error->one(FLERR,"Curvature limiter cannot be negative!");
      if(!meanCurvature && !gaussianCurvature)
        curvatureLimitFactor = 0.0;
      if(meanCurvature && gaussianCurvature)
        error->one(FLERR,"meanCurvature and gaussianCurvature cannot be simultaneously yes !");
    }

    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}
    void connectToProperties(PropertyRegistry&) {}

    bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
    {
      sidata.is_non_spherical = true;
      bool particles_in_contact = false;
      double *const prev_step_point = &sidata.contact_history[contact_point_offset]; //contact points
      double *const inequality_start = &sidata.contact_history[inequality_start_offset];
      double *const particles_were_in_contact = &sidata.contact_history[particles_were_in_contact_offset];

      const int iPart = sidata.i;
      const int jPart = sidata.j;

#ifdef LIGGGHTS_DEBUG
      if(std::isnan(vectorMag3D(atom->x[iPart])))
        error->one(FLERR,"atom->x[iPart] is NaN!");
      if(std::isnan(vectorMag4D(atom->quaternion[iPart])))
        error->one(FLERR,"atom->quaternion[iPart] is NaN!");
      if(std::isnan(vectorMag3D(atom->x[jPart])))
        error->one(FLERR,"atom->x[jPart] is NaN!");
      if(std::isnan(vectorMag4D(atom->quaternion[jPart])))
        error->one(FLERR,"atom->quaternion[jPart] is NaN!");
#endif
      particle_i.set(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
      particle_j.set(atom->x[jPart], atom->quaternion[jPart], atom->shape[jPart], atom->blockiness[jPart]);

      unsigned int int_inequality_start = MathExtraLiggghtsNonspherical::round_int(*inequality_start);

      bool obb_intersect = false;
      if(*particles_were_in_contact == SURFACES_INTERSECT)
        obb_intersect = true; //particles had overlap on the previous time step, skipping OBB intersection check
      else
        obb_intersect = MathExtraLiggghtsNonspherical::obb_intersect(&particle_i, &particle_j, int_inequality_start);
      if(obb_intersect) {//OBB intersect particles in possible contact

        double fi, fj;
        const double ri = cbrt(particle_i.shape[0]*particle_i.shape[1]*particle_i.shape[2]);
        const double rj = cbrt(particle_j.shape[0]*particle_j.shape[1]*particle_j.shape[2]);
        double ratio = ri / (ri + rj);

        double contact_point[3];
        if(*particles_were_in_contact == SURFACES_FAR)
          calc_contact_point_if_no_previous_point_avaialable(sidata, &particle_i, &particle_j, contact_point, fi, fj,this->error);
        else
          calc_contact_point_using_prev_step(sidata, &particle_i, &particle_j, ratio, update->dt, prev_step_point, contact_point, fi, fj, this->error);
        vectorCopy3D(contact_point, prev_step_point); //save contact point in contact history for the next DEM time step

#ifdef LIGGGHTS_DEBUG
        if(std::isnan(vectorMag3D(sidata.contact_point)))
          error->one(FLERR,"sidata.contact_point is NaN!");
#endif
        particles_in_contact = std::max(fi, fj) < 0.0;

        if(particles_in_contact) {
          double contact_point_i[3], contact_point_j[3], en[3];
          vectorSubtract3D(particle_j.gradient, particle_i.gradient, en);
          vectorNormalize3D(en);

          double *const alpha_i = &sidata.contact_history[alpha1_offset];
          double *const alpha_j = &sidata.contact_history[alpha2_offset];

          double delta[3];
          extended_overlap_algorithm(&particle_i, &particle_j, en, alpha_i, alpha_j,
              contact_point, contact_point_i, contact_point_j, delta);

          double *contact_point_i_local = &sidata.contact_history[contact_point_i_offset];
          double *contact_point_j_local = &sidata.contact_history[contact_point_j_offset];

          bool particles_were_in_contact_flag = (*particles_were_in_contact == SURFACES_INTERSECT);

          sidata.deltan = common_normal(sidata, &particle_i, &particle_j, particles_were_in_contact_flag,
                 contact_point_i_local, contact_point_j_local, contact_point_i, contact_point_j);

          sidata.reff = sidata.radi*sidata.radj / (sidata.radi + sidata.radj);
          if(curvatureLimitFactor > 0.0) {
              int curvature_radius_method = meanCurvature ? 0 : 1;
              double koefi = particle_i.calc_curvature_coefficient(curvature_radius_method, contact_point_i); //mean curvature coefficient
              double koefj = particle_j.calc_curvature_coefficient(curvature_radius_method, contact_point_j); //mean curvature coefficient
#ifdef LIGGGHTS_DEBUG
              if(std::isnan(koefi))
                error->one(FLERR,"sidata.koefi is NaN!");
              if(std::isnan(koefj))
                error->one(FLERR,"sidata.koefj is NaN!");
              if(std::isnan(sidata.reff))
                error->one(FLERR,"sidata.koefj is NaN!");
#endif
              sidata.reff = get_effective_radius(sidata, atom->blockiness[iPart], atom->blockiness[jPart], koefi, koefj, curvatureLimitFactor, this->error);
          }
#ifdef LIGGGHTS_DEBUG
          if(std::isnan(vectorMag3D(contact_point_i)))
            error->one(FLERR,"sidata.contact_point_i is NaN!");
          if(std::isnan(vectorMag3D(contact_point_j)))
            error->one(FLERR,"sidata.contact_point_j is NaN!");

          if(std::isnan(*alpha_i))
            error->one(FLERR,"alpha_i is NaN!");
          if(std::isnan(*alpha_j))
            error->one(FLERR,"alpha_j is NaN!");
          if(std::isnan(sidata.deltan))
            error->one(FLERR,"sidata.deltan is NaN!");
#endif

          *particles_were_in_contact = SURFACES_INTERSECT;
        } else
          *particles_were_in_contact = SURFACES_CLOSE;
       } else
         *particles_were_in_contact = SURFACES_FAR;
      *inequality_start = static_cast<double>(int_inequality_start);
      return particles_in_contact;
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData&, ForceData&)
    {
      if(sidata.is_wall) {
        sidata.reff = sidata.radi;
        if(curvatureLimitFactor > 0.0) {
          const int iPart = sidata.i;
          particle_i.set(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
          int curvature_radius_method = meanCurvature ? 0 : 1;
          double koefi = particle_i.calc_curvature_coefficient(curvature_radius_method, sidata.contact_point); //mean curvature coefficient
          sidata.reff = get_effective_radius_wall(sidata, atom->blockiness[iPart], koefi, curvatureLimitFactor, this->error);
        }
      }
      MathExtraLiggghtsNonspherical::surfacesIntersectNonSpherical(sidata, atom->x);
    }

    void endSurfacesIntersect(SurfacesIntersectData&,TriMesh *, double * const) {}
    void surfacesClose(SurfacesCloseData&, ForceData&, ForceData&){}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void tally_pp(double,int,int,int) {}
    void tally_pw(double,int,int,int) {}

  protected:
     double curvatureLimitFactor;
     bool meanCurvature;
     bool gaussianCurvature;
  };
}
}

#else // SUPERQUADRIC_ACTIVE_FLAG

// add dummy class
SURFACE_MODEL_DUMMY(SURFACE_SUPERQUADRIC_ORTHOGONAL, "Surface model used without support for superquadric particles")

#endif // SUPERQUADRIC_ACTIVE_FLAG

#endif // SURFACE_MODEL_SUPERQUADRIC_H_
#endif // SURFACE_MODEL
