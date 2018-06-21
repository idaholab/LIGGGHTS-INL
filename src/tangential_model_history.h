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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef TANGENTIAL_MODEL
TANGENTIAL_MODEL(TANGENTIAL_HISTORY,history,2)
#else
#ifndef TANGENTIAL_MODEL_HISTORY_H_
#define TANGENTIAL_MODEL_HISTORY_H_
#include "contact_models.h"
#include "tangential_model_base.h"
#include <cmath>
#include "update.h"
#include "global_properties.h"
#include "atom.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class TangentialModel<TANGENTIAL_HISTORY> : public TangentialModelBase
  {
  public:
    TangentialModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
        TangentialModelBase(lmp, hsetup, c),
        coeffFrict(NULL),
        heating(false),
        heating_track(false),
        cmb(c),
        iradius_offset_(-1),
        elastic_potential_offset_(-1),
        elasticpotflag_(false),
        dissipation_history_offset_(-1),
        dissipatedflag_(false),
        fix_dissipated_(NULL),
        torsionflag_(false),
        torsion_offset_(-1)
    {
        history_offset = hsetup->add_history_value("shearx", "1");
        hsetup->add_history_value("sheary", "1");
        hsetup->add_history_value("shearz", "1");

    }

    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
        iradius_offset_ = cmb->get_history_offset("intersection_radius");
        if (elasticpotflag_)
        {
            elastic_potential_offset_ = cmb->get_history_offset("elastic_potential_normal");
            if (elastic_potential_offset_ == -1)
                error->all(FLERR, "Require normal model with elastic potential computation");
        }
        if (dissipatedflag_)
        {
            if (cmb->is_wall())
            {
                fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_wall_tangential", "property/atom", "vector", 0, 0, "dissipated energy"));
                dissipation_history_offset_ = cmb->get_history_offset("dissipation_force");
                if (!dissipation_history_offset_)
                    error->one(FLERR, "Internal error: Could not find dissipation history offset");
            }
            else
                fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_tangential", "property/atom", "vector", 0, 0, "dissipated energy"));
            if (!fix_dissipated_)
                error->one(FLERR, "Surface model has not registered dissipated_energy fix");
        }
        if (torsionflag_)
        {
            torsion_offset_ = hsetup->add_history_value("shear_omega_n", "1");
            if (elasticpotflag_ || dissipatedflag_)
                error->all(FLERR, "Elastic and dissipated energy calculation not implemented for torsion");
        }
    }

    void registerSettings(Settings& settings)
    {
        settings.registerOnOff("heating_tangential_history",heating,false);
        settings.registerOnOff("heating_tracking",heating_track,false);
        settings.registerOnOff("computeElasticPotential", elasticpotflag_, false);
        settings.registerOnOff("computeDissipatedEnergy", dissipatedflag_, false);
        settings.registerOnOff("torsion", torsionflag_, false);
        //TODO error->one(FLERR,"TODO here also check if right surface model used");
    }

    void connectToProperties(PropertyRegistry & registry)
    {
        registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
        registry.connect("coeffFrict", coeffFrict,"tangential_model history");
        if ((elasticpotflag_ || dissipatedflag_) && cmb->is_wall())
            
            error->warning(FLERR, "Energy conservation might not be computed correctly if mesh/surface is used without the weighted_wall_formulation keyword turned off (default).");
    }

    void surfacesIntersect(const SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
        // normal forces = Hookian contact + normal velocity damping
        const double enx = sidata.en[0];
        const double eny = sidata.en[1];
        const double enz = sidata.en[2];

        if(sidata.contact_flags)
            *sidata.contact_flags |= CONTACT_TANGENTIAL_MODEL;

        // shear history effects
        double * const shear = &sidata.contact_history[history_offset];
        double shear_old[3];
        const bool update_history = sidata.computeflag && sidata.shearupdate;
        if (update_history && elasticpotflag_)
            vectorCopy3D(shear, shear_old);

        if (update_history)
        {
            const double dt = update->dt;
            shear[0] += sidata.vtr1 * dt;
            shear[1] += sidata.vtr2 * dt;
            shear[2] += sidata.vtr3 * dt;

            // rotate shear displacements

            double rsht = shear[0]*enx + shear[1]*eny + shear[2]*enz;
            shear[0] -= rsht * enx;
            shear[1] -= rsht * eny;
            shear[2] -= rsht * enz;

            // add rotation angle
            if (torsionflag_)
                sidata.contact_history[torsion_offset_] += dt*(enx*sidata.wr1+eny*sidata.wr2+enz*sidata.wr3);
        }

        const double shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);
        const double kt = sidata.kt;
        const double xmu = coeffFrict[sidata.itype][sidata.jtype];

        // tangential forces = shear + tangential velocity damping
        double Ft_ela1 = -(kt * shear[0]);
        double Ft_ela2 = -(kt * shear[1]);
        double Ft_ela3 = -(kt * shear[2]);
        double Ft1 = Ft_ela1;
        double Ft2 = Ft_ela2;
        double Ft3 = Ft_ela3;

        // rescale frictional displacements and forces if needed
        const double Ft_shear = kt * shrmag; // sqrt(Ft1 * Ft1 + Ft2 * Ft2 + Ft3 * Ft3);
        const double Ft_friction = xmu * fabs(sidata.Fn);

        // energy loss from sliding or damping
        if (Ft_shear > Ft_friction)
        {
            if (shrmag != 0.0)
            {
                const double ratio = Ft_friction / Ft_shear;
                
                if(heating)
                {
                    const double P_diss_local = (Ft_shear - Ft_friction)*(Ft_shear + Ft_friction) / (update->dt*kt); 
                    sidata.P_diss += P_diss_local;
                    if(heating_track && sidata.is_wall)
                        cmb->tally_pw(P_diss_local, sidata.i, sidata.jtype, 2);
                    if(heating_track && !sidata.is_wall)
                        cmb->tally_pp(P_diss_local, sidata.i, sidata.j, 2);
                }
                Ft1 *= ratio;
                Ft2 *= ratio;
                Ft3 *= ratio;
                
                if (update_history)
                {
                    shear[0] = -Ft1/kt;
                    shear[1] = -Ft2/kt;
                    shear[2] = -Ft3/kt;
                    if (elasticpotflag_ || dissipatedflag_)
                    {
                        
                        const double weight = 1.0 - vectorMag3D(shear_old)/shrmag;
                        Ft_ela1 = Ft1 * weight;
                        Ft_ela2 = Ft2 * weight;
                        Ft_ela3 = Ft3 * weight;
                    }
                }
            }
            else
                Ft1 = Ft2 = Ft3 = 0.0;
        }
        else
        {
            const double gammat = sidata.gammat;
            Ft1 -= (gammat*sidata.vtr1);
            Ft2 -= (gammat*sidata.vtr2);
            Ft3 -= (gammat*sidata.vtr3);
            if(heating)
            {
                const double P_diss_local = gammat*(sidata.vtr1*sidata.vtr1+sidata.vtr2*sidata.vtr2+sidata.vtr3*sidata.vtr3); 
                sidata.P_diss += P_diss_local;
                if(heating_track && sidata.is_wall)
                    cmb->tally_pw(P_diss_local, sidata.i, sidata.jtype, 1);
                if(heating_track && !sidata.is_wall)
                    cmb->tally_pp(P_diss_local, sidata.i, sidata.j, 1);
            }
        }

        double Tn_shear = 0.;
        if (torsionflag_)
        {
            // moments
            // radius of intersection area
            double cap_radius = 0.;
            if (atom->shapetype_flag && iradius_offset_ >= 0)
                cap_radius = sidata.contact_history[iradius_offset_];
            else
            {
                // for spheres it is the radius of the spherical cap
                if (sidata.is_wall)
                {
                    const double cap_radi = sqrt(sidata.radi*sidata.radi - (sidata.radi-sidata.deltan)*(sidata.radi-sidata.deltan));
                    cap_radius = cap_radi;
                }
                else
                {
                    const double half_delta = sidata.deltan*0.5;
                    const double cap_radi = sqrt(sidata.radi*sidata.radi - (sidata.radi-half_delta)*(sidata.radi-half_delta));
                    const double cap_radj = sqrt(sidata.radj*sidata.radj - (sidata.radj-half_delta)*(sidata.radj-half_delta));
                    cap_radius = (cap_radi+cap_radj)*0.5;
                }
            }
            Tn_shear = -kt*cap_radius*cap_radius*shear[3];
            const double Tn_friction = Ft_friction*cap_radius;

            // sliding friction
            const double abs_Tn_shear = fabs(Tn_shear);
            if (abs_Tn_shear > Tn_friction)
            {
                if (abs_Tn_shear > 1e-16*Tn_friction)
                {
                    const double factor = Tn_friction/abs_Tn_shear;
                    Tn_shear *= factor;
                    sidata.contact_history[torsion_offset_] = Tn_shear/(-kt*cap_radius*cap_radius);
                }
                else
                    Tn_shear = 0.;
            }
            // static friction (with damping)
            else
                Tn_shear -= sidata.gammat*cap_radius*cap_radius*(enx*sidata.wr1+eny*sidata.wr2+enz*sidata.wr3);
        }

        // forces & torques
        const double tor1 = eny * Ft3 - enz * Ft2;
        const double tor2 = enz * Ft1 - enx * Ft3;
        const double tor3 = enx * Ft2 - eny * Ft1;

        double torque_i[3];
        double torque_j[3] = {0.,0.,0.};
        #ifdef NONSPHERICAL_ACTIVE_FLAG
        if(sidata.is_non_spherical)
        {
            double xci[3];
            double Ft_i[3] = { Ft1,  Ft2,  Ft3 };
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Ft_i, torque_i);
            if (torsionflag_)
            {
                torque_i[0] += Tn_shear*enx;
                torque_i[1] += Tn_shear*eny;
                torque_i[2] += Tn_shear*enz;
            }
            if (!sidata.is_wall)
            {
                double xcj[3];
                vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
                double Ft_j[3] = { -Ft1,  -Ft2,  -Ft3 };
                vectorCross3D(xcj, Ft_j, torque_j);
                if (torsionflag_)
                {
                    torque_j[0] -= Tn_shear*enx;
                    torque_j[1] -= Tn_shear*eny;
                    torque_j[2] -= Tn_shear*enz;
                }
            }
        }
        else
        #endif
        {
            torque_i[0] = -sidata.cri * tor1 + Tn_shear*enx;
            torque_i[1] = -sidata.cri * tor2 + Tn_shear*eny;
            torque_i[2] = -sidata.cri * tor3 + Tn_shear*enz;
            if (!sidata.is_wall)
            {
                torque_j[0] = -sidata.crj * tor1 - Tn_shear*enx;
                torque_j[1] = -sidata.crj * tor2 - Tn_shear*eny;
                torque_j[2] = -sidata.crj * tor3 - Tn_shear*enz;
            }
        }

        if (update_history && (elasticpotflag_ || dissipatedflag_))
        {
            double torque_ela_i[3], torque_ela_j[3];
            torque_ela_i[0] = eny * Ft_ela3 - enz * Ft_ela2;
            torque_ela_i[1] = enz * Ft_ela1 - enx * Ft_ela3;
            torque_ela_i[2] = enx * Ft_ela2 - eny * Ft_ela1;
            vectorScalarMult3D(torque_ela_i, -sidata.crj, torque_ela_j);
            vectorScalarMult3D(torque_ela_i, -sidata.cri);
            if (elasticpotflag_)
            {
                double * const elastic_pot = &sidata.contact_history[elastic_potential_offset_];
                if (sidata.is_wall)
                {
                    double delta[3] = {0.,0.,0.};
                    vectorCopy3D(sidata.v_j,delta);
                    vectorScalarMult3D(delta, update->dt);
                    
                    elastic_pot[10] -= (delta[0]*Ft_ela1 +
                                        delta[1]*Ft_ela2 +
                                        delta[2]*Ft_ela3)*0.5;
                }
                elastic_pot[1] -= Ft_ela1;
                elastic_pot[2] -= Ft_ela2;
                elastic_pot[3] -= Ft_ela3;
                elastic_pot[4] -= torque_ela_i[0];
                elastic_pot[5] -= torque_ela_i[1];
                elastic_pot[6] -= torque_ela_i[2];
                elastic_pot[7] -= torque_ela_j[0];
                elastic_pot[8] -= torque_ela_j[1];
                elastic_pot[9] -= torque_ela_j[2];
            }
            if (dissipatedflag_)
            {
                double * const * const dissipated = fix_dissipated_->array_atom;

                double * const dissipated_i = dissipated[sidata.i];
                double * const dissipated_j = dissipated[sidata.j];
                dissipated_i[1] += -(Ft1 - Ft_ela1);
                dissipated_i[2] += -(Ft2 - Ft_ela2);
                dissipated_i[3] += -(Ft3 - Ft_ela3);
                const double dTorqueDamp_i[3] = {torque_i[0]-torque_ela_i[0], torque_i[1]-torque_ela_i[1], torque_i[2]-torque_ela_i[2]};
                dissipated_i[4] += -dTorqueDamp_i[0];
                dissipated_i[5] += -dTorqueDamp_i[1];
                dissipated_i[6] += -dTorqueDamp_i[2];
                if (sidata.j < atom->nlocal && !sidata.is_wall)
                {
                    dissipated_j[1] -= -(Ft1 - Ft_ela1);
                    dissipated_j[2] -= -(Ft2 - Ft_ela2);
                    dissipated_j[3] -= -(Ft3 - Ft_ela3);
                    const double dTorqueDamp_j[3] = {torque_j[0]-torque_ela_j[0], torque_j[1]-torque_ela_j[1], torque_j[2]-torque_ela_j[2]};
                    dissipated_j[4] += -dTorqueDamp_j[0];
                    dissipated_j[5] += -dTorqueDamp_j[1];
                    dissipated_j[6] += -dTorqueDamp_j[2];
                }
                else if (sidata.is_wall)
                {
                    double * const diss_force = &sidata.contact_history[dissipation_history_offset_];
                    diss_force[0] += Ft1 - Ft_ela1;
                    diss_force[1] += Ft2 - Ft_ela2;
                    diss_force[2] += Ft3 - Ft_ela3;
                }
            }
        }

        // return resulting forces
        if(sidata.is_wall)
        {
            const double area_ratio = sidata.area_ratio;
            i_forces.delta_F[0] += Ft1 * area_ratio;
            i_forces.delta_F[1] += Ft2 * area_ratio;
            i_forces.delta_F[2] += Ft3 * area_ratio;
            i_forces.delta_torque[0] += torque_i[0] * area_ratio;
            i_forces.delta_torque[1] += torque_i[1] * area_ratio;
            i_forces.delta_torque[2] += torque_i[2] * area_ratio;
        }
        else
        {
            i_forces.delta_F[0] += Ft1;
            i_forces.delta_F[1] += Ft2;
            i_forces.delta_F[2] += Ft3;
            j_forces.delta_F[0] += -Ft1;
            j_forces.delta_F[1] += -Ft2;
            j_forces.delta_F[2] += -Ft3;
            i_forces.delta_torque[0] += torque_i[0];
            i_forces.delta_torque[1] += torque_i[1];
            i_forces.delta_torque[2] += torque_i[2];

            j_forces.delta_torque[0] += torque_j[0];
            j_forces.delta_torque[1] += torque_j[1];
            j_forces.delta_torque[2] += torque_j[2];
        }
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
        // unset non-touching neighbors
        // TODO even if shearupdate == false?
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_TANGENTIAL_MODEL;
        if(!scdata.contact_history)
          return; //DO NOT access contact_history if not available
        double * const shear = &scdata.contact_history[history_offset];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

   protected:
    double ** coeffFrict;
    int history_offset;
    bool heating;
    bool heating_track;
    class ContactModelBase *cmb;
    int iradius_offset_;
    int elastic_potential_offset_;
    bool elasticpotflag_;
    int dissipation_history_offset_;
    bool dissipatedflag_;
    FixPropertyAtom *fix_dissipated_;
    bool torsionflag_;
    int torsion_offset_;
  };
}
}
#endif // TANGENTIAL_MODEL_HISTORY_H_
#endif
