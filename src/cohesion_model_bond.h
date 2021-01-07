/* ----------------------------------------------------------------------
   Copyright 2021, Battelle Energy Alliance, LLC  All Rights Reserved
-------------------------------------------------------------------------

   Contributing author and copyright for this file:

   Yidong Xia (Idaho National Laboratory) - fixed minor bugs
------------------------------------------------------------------------- */
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

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_BOND,bond,14)
#else

#ifndef COHESION_MODEL_BOND_H_
#define COHESION_MODEL_BOND_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include "compute_bond_counter.h"
#include <cmath>
#include <algorithm>
#include "error.h"
#include "math_extra_liggghts.h"
#include "fix_property_atom.h"
#include "neighbor.h"

static const double SMALL_COHESION_MODEL_BOND = 1.e-6;

// required model properties
namespace MODEL_PARAMS
{
  static const char * DISSIPATION_NORM_FORCE_BOND = "dissipationNormalForceBond";
  static const char * DISSIPATION_TANG_FORCE_BOND = "dissipationTangentialForceBond";
  static const char * DISSIPATION_NORM_TORQUE_BOND = "dissipationNormalTorqueBond";
  static const char * DISSIPATION_TANG_TORQUE_BOND = "dissipationTangentialTorqueBond";

  inline static MatrixProperty* createRadiusMultiplierBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "radiusMultiplierBond", caller, sanity_checks, 0.0);
  }

  /* ---------------------------------------------------------------------- */

  inline static MatrixProperty* createNormalBondStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "normalBondStiffnessPerUnitArea", caller, sanity_checks, 0.0);
  }

  /* ---------------------------------------------------------------------- */

  inline static MatrixProperty* createTangentialBondStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "tangentialBondStiffnessPerUnitArea", caller, sanity_checks, 0.0);
  }

  /* ---------------------------------------------------------------------- */

  inline static MatrixProperty* createDissipationNormalForceBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createDissipationMatrix(registry, DISSIPATION_NORM_FORCE_BOND, caller, sanity_checks);
  }

  inline static MatrixProperty* createDissipationTangentialForceBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createDissipationMatrix(registry, DISSIPATION_TANG_FORCE_BOND, caller, sanity_checks);
  }

  inline static MatrixProperty* createDissipationNormalTorqueBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createDissipationMatrix(registry, DISSIPATION_NORM_TORQUE_BOND, caller, sanity_checks);
  }

  inline static MatrixProperty* createDissipationTangentialTorqueBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createDissipationMatrix(registry, DISSIPATION_TANG_TORQUE_BOND, caller, sanity_checks);
  }

  /* ---------------------------------------------------------------------- */

  inline static MatrixProperty* createDampingNormalForceBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "dampingNormalForceBond", caller, sanity_checks, 0.0, 1.0);
  }

  inline static MatrixProperty* createDampingTangentialForceBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "dampingTangentialForceBond", caller, sanity_checks, 0.0, 1.0);
  }

  inline static MatrixProperty* createDampingNormalTorqueBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "dampingNormalTorqueBond", caller, sanity_checks, 0.0, 1.0);
  }

  inline static MatrixProperty* createDampingTangentialTorqueBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "dampingTangentialTorqueBond", caller, sanity_checks, 0.0, 1.0);
  }

  /* ---------------------------------------------------------------------- */
  // add additional properties without senity check
  /* ---------------------------------------------------------------------- */
  inline static MatrixProperty* createMaxDistanceBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "maxDistanceBond", caller, sanity_checks, 0.0);
  }

  inline static MatrixProperty* createMaxSigmaBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "maxSigmaBond", caller, sanity_checks, 0.0);
  }

  inline static MatrixProperty* createMaxTauBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "maxTauBond", caller, sanity_checks, 0.0);
  }

  inline static MatrixProperty* createRatioTensionCompression(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "ratioTensionCompression", caller, sanity_checks, 0.0);
  }

  inline static MatrixProperty* createFrictionAngle(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "frictionAngle", caller, sanity_checks, 0.0);
  }

  /* ---------------------------------------------------------------------- */

  inline static ScalarProperty* createTsCreateBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createScalarProperty(registry, "tsCreateBond", caller, sanity_checks, 0.0);
  }

  inline static MatrixProperty* createCreateDistanceBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "createDistanceBond", caller, sanity_checks, 0.0);
  }

} // ns MODEL_PARAMS

/* ---------------------------------------------------------------------- */

namespace LIGGGHTS {

namespace ContactModels {
  enum{
       BREAKSTYLE_SIMPLE,
       BREAKSTYLE_STRESS,
       BREAKSTYLE_STRESS_TEMP
      };

  template<>
  class CohesionModel<COHESION_BOND> : public CohesionModelBase {

  public:

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
      CohesionModelBase(lmp, hsetup, cmb),
      history_offset_(0),
      elastic_potential_offset_(-1),
      dissipation_history_offset_(-1),
      fix_dissipated_(NULL),
      radiusMultiplier_(NULL),
      normalBondStiffness_(NULL),
      tangentialBondStiffness_(NULL),
      fn_dissipation_(NULL),
      ft_dissipation_(NULL),
      tn_dissipation_(NULL),
      tt_dissipation_(NULL),
      maxDist_(NULL),
      maxSigma_(NULL),
      maxTau_(NULL),
      ratioTensionCompression_(NULL),
      tanFrictionAngle_(NULL),
      tsCreateBond_(0.0),
      createDistanceBond_(NULL),
      stressbreakflag_(false),
      heatbreakflag_(false),
      tensionflag_(true),
      compressionflag_(true),
      shearflag_(true),
      ntorqueflag_(true),
      ttorqueflag_(true),
      createbond_always_flag_(false),
      dampingflag_(true),
      dampingsmooth_(false),
      dissipationflag_(false),
      elasticpotflag_(false),
      dissipatedflag_(false),
      breakmode_(BREAKSTYLE_SIMPLE),
      useRatioTC_(false),
      druckerPrager_(false),
      fix_strenghtening_factor_(NULL),
      fix_bond_random_id_(NULL),
      compute_bond_counter_(NULL)
    {
      history_offset_ = hsetup->add_history_value("contflag", "0");
      hsetup->add_history_value("initial_dist", "0");
      hsetup->add_history_value("contactPosX", "0");
      hsetup->add_history_value("contactPosY", "0");
      hsetup->add_history_value("contactPosZ", "0");
      hsetup->add_history_value("ftx", "1");
      hsetup->add_history_value("fty", "1");
      hsetup->add_history_value("ftz", "1");
      hsetup->add_history_value("torquenx", "1");
      hsetup->add_history_value("torqueny", "1");
      hsetup->add_history_value("torquenz", "1");
      hsetup->add_history_value("torquetx", "1");
      hsetup->add_history_value("torquety", "1");
      hsetup->add_history_value("torquetz", "1");

      cmb->add_history_offset("bond_contactflag", history_offset_);
    }

    /* ---------------------------------------------------------------------- */

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("stressBreak", stressbreakflag_, false);
      settings.registerOnOff("temperatureBreak", heatbreakflag_, false);
      settings.registerOnOff("tensionStress", tensionflag_, true);
      settings.registerOnOff("compressionStress", compressionflag_, true);
      settings.registerOnOff("shearStress", shearflag_, true);
      settings.registerOnOff("normalTorqueStress", ntorqueflag_, true);
      settings.registerOnOff("shearTorqueStress", ttorqueflag_, true);

      settings.registerOnOff("createBondAlways", createbond_always_flag_, false); // if false, bonds are created if within range, else at a certain timestep and within range
      settings.registerOnOff("dampingBond", dampingflag_, true);
      settings.registerOnOff("dampingBondSmooth", dampingsmooth_, false);
      settings.registerOnOff("dissipationBond", dissipationflag_, false);
      settings.registerOnOff("computeElasticPotential", elasticpotflag_, false);
      settings.registerOnOff("computeDissipatedEnergy", dissipatedflag_, false);
      settings.registerOnOff("ratioTensionCompression", useRatioTC_, false);
      settings.registerOnOff("druckerPrager", druckerPrager_, false);
    }

    /* ---------------------------------------------------------------------- */

    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
      if (elasticpotflag_)
      {
        elastic_potential_offset_ = cmb->get_history_offset("elastic_potential_cohesion");
        if (elastic_potential_offset_ == -1)
        {
          elastic_potential_offset_ = hsetup->add_history_value("elastic_potential_cohesion", "0");
          hsetup->add_history_value("elastic_force_cohesion_0", "1");
          hsetup->add_history_value("elastic_force_cohesion_1", "1");
          hsetup->add_history_value("elastic_force_cohesion_2", "1");
          hsetup->add_history_value("elastic_torque_cohesion_i_0", "0");
          hsetup->add_history_value("elastic_torque_cohesion_i_1", "0");
          hsetup->add_history_value("elastic_torque_cohesion_i_2", "0");
          hsetup->add_history_value("elastic_torque_cohesion_j_0", "0");
          hsetup->add_history_value("elastic_torque_cohesion_j_1", "0");
          hsetup->add_history_value("elastic_torque_cohesion_j_2", "0");
          if (cmb->is_wall())
              hsetup->add_history_value("elastic_potential_wall", "0");
          cmb->add_history_offset("elastic_potential_cohesion", elastic_potential_offset_);
        }
      }
      if (dissipatedflag_)
      {
        if (cmb->is_wall())
        {
            fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_wall", "property/atom", "vector", 0, 0, "dissipated energy"));
            dissipation_history_offset_ = cmb->get_history_offset("dissipation_force");
            if (!dissipation_history_offset_)
                error->one(FLERR, "Internal error: Could not find dissipation history offset");
        }
        else
            fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy", "property/atom", "vector", 0, 0, "dissipated energy"));
        if (!fix_dissipated_)
            error->one(FLERR, "Surface model has not registered dissipated_energy fix");
      }
    }

    /* ---------------------------------------------------------------------- */

    void connectToProperties(PropertyRegistry & registry)
    {
      // check if fix_strenghtening_factor or fix_bond_random_id exist
      fix_strenghtening_factor_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("bond_strengthening_factor","property/atom","vector",0,0,"cohesion bond",false));
      fix_bond_random_id_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("bond_random_id","property/atom","scalar",0,0,"cohesion bond",false));

      registry.registerProperty("radiusMultiplierBond", &MODEL_PARAMS::createRadiusMultiplierBond);
      registry.registerProperty("normalBondStiffness", &MODEL_PARAMS::createNormalBondStiffness);
      registry.registerProperty("tangentialBondStiffness", &MODEL_PARAMS::createTangentialBondStiffness);

      registry.connect("radiusMultiplierBond", radiusMultiplier_,"cohesion bond");
      registry.connect("normalBondStiffness", normalBondStiffness_,"cohesion bond");
      registry.connect("tangentialBondStiffness", tangentialBondStiffness_,"cohesion bond");

      if (dampingsmooth_)
        dampingflag_ = true;

      if (!dampingflag_ && !dissipationflag_)
        error->one(FLERR,"Damping or dissipation has to be enabled.");

      if (dampingflag_)
      {
        registry.registerProperty("dampingNormalForceBond", &MODEL_PARAMS::createDampingNormalForceBond);
        registry.registerProperty("dampingTangentialForceBond", &MODEL_PARAMS::createDampingTangentialForceBond);
        registry.registerProperty("dampingNormalTorqueBond", &MODEL_PARAMS::createDampingNormalTorqueBond);
        registry.registerProperty("dampingTangentialTorqueBond", &MODEL_PARAMS::createDampingTangentialTorqueBond);

        registry.connect("dampingNormalForceBond", fn_damping_,"cohesion bond");
        registry.connect("dampingTangentialForceBond", ft_damping_,"cohesion bond");
        registry.connect("dampingNormalTorqueBond", tn_damping_,"cohesion bond");
        registry.connect("dampingTangentialTorqueBond", tt_damping_,"cohesion bond");
      }

      if (dissipationflag_)
      {
        registry.registerProperty("dissipationNormalForceBond", &MODEL_PARAMS::createDissipationNormalForceBond);
        registry.registerProperty("dissipationTangentialForceBond", &MODEL_PARAMS::createDissipationTangentialForceBond);
        registry.registerProperty("dissipationNormalTorqueBond", &MODEL_PARAMS::createDissipationNormalTorqueBond);
        registry.registerProperty("dissipationTangentialTorqueBond", &MODEL_PARAMS::createDissipationTangentialTorqueBond);

        registry.connect("dissipationNormalForceBond", fn_dissipation_,"cohesion bond");
        registry.connect("dissipationTangentialForceBond", ft_dissipation_,"cohesion bond");
        registry.connect("dissipationNormalTorqueBond", tn_dissipation_,"cohesion bond");
        registry.connect("dissipationTangentialTorqueBond", tt_dissipation_,"cohesion bond");
        if (dissipatedflag_)
        {
            const int max_type = registry.max_type();
            for (int i=1; i<max_type+1; i++)
            {
                for(int j=1;j<max_type+1;j++)
                {
                    if(!MathExtraLiggghts::compDouble(fn_dissipation_[i][j], ft_dissipation_[i][j]) ||
                       !MathExtraLiggghts::compDouble(fn_dissipation_[i][j], tn_dissipation_[i][j]) ||
                       !MathExtraLiggghts::compDouble(fn_dissipation_[i][j], tt_dissipation_[i][j])   )
                        error->one(FLERR, "Dissipated energy can only be tracked if all dissipation parameters are equal");
                }
            }
        }
      }

      // create bond settings
      if (createbond_always_flag_)
      {
        registry.registerProperty("createDistanceBond", &MODEL_PARAMS::createCreateDistanceBond);
        registry.connect("createDistanceBond", createDistanceBond_, "cohesion bond");
      }
      else
      {
        // required settings
        registry.registerProperty("tsCreateBond", &MODEL_PARAMS::createTsCreateBond);
        registry.registerProperty("createDistanceBond", &MODEL_PARAMS::createCreateDistanceBond);
        registry.connect("tsCreateBond", tsCreateBond_, "cohesion bond");
        registry.connect("createDistanceBond", createDistanceBond_, "cohesion bond");
      }

      // break up settings
      // TODO: Optimize the contactDistanceFactor; minrad as per-type property
      double cdf_all = 1.; // minimum value is one

      const int max_type = registry.max_type();
      const double minrad = registry.min_radius();

      if (minrad <= 0.)
        error->one(FLERR,"Bond settings: The minimum radius can't be <= 0!");

      if (!stressbreakflag_)
      {
        // required settings
        registry.registerProperty("maxDistanceBond", &MODEL_PARAMS::createMaxDistanceBond);
        registry.connect("maxDistanceBond", maxDist_,"cohesion bond");
        breakmode_ = BREAKSTYLE_SIMPLE;

        // calculate the contactDistanceFactor
        for (int i=1; i<max_type+1; i++)
        {
          for(int j=1;j<max_type+1;j++)
          {

            double cdf_one = 1.1 * 0.5*maxDist_[i][j]/minrad;

            cdf_all = cdf_one > cdf_all ? cdf_one : cdf_all;

          }
        }
      }
      else
      {
        // required settings
        registry.registerProperty("maxSigmaBond", &MODEL_PARAMS::createMaxSigmaBond);
        registry.registerProperty("maxTauBond", &MODEL_PARAMS::createMaxTauBond);
        registry.connect("maxSigmaBond", maxSigma_,"cohesion bond");
        registry.connect("maxTauBond", maxTau_,"cohesion bond");
        if (useRatioTC_)
        {
          registry.registerProperty("ratioTensionCompression", &MODEL_PARAMS::createRatioTensionCompression);
          registry.connect("ratioTensionCompression", ratioTensionCompression_,"cohesion bond");
        }
        if (druckerPrager_)
        {
          registry.registerProperty("frictionAngle", &MODEL_PARAMS::createFrictionAngle);
          registry.connect("frictionAngle", tanFrictionAngle_,"cohesion bond");
          for (int i=1; i<max_type+1; i++)
            for (int j=1; j<max_type+1; j++)
                tanFrictionAngle_[i][j] = tan(tanFrictionAngle_[i][j]*M_PI/180.0);
        }
        breakmode_ = BREAKSTYLE_STRESS;

        // calculate the contactDistanceFactor
        // it depends only on the normalStress due to normal distance!
        //  --> maxSigma_, normalBondStiffness_, minrad!
        for (int i=1; i<max_type+1; i++)
        {
          for (int j=1; j<max_type+1; j++)
          {

            double stress = maxSigma_[i][j];
            if (useRatioTC_)
                stress = fmax(stress, stress*ratioTensionCompression_[i][j]);
            if(normalBondStiffness_[i][j] <= 1e-15)
                error->one(FLERR,"Bond settings: In case of stress breakage, the normal bond stiffness can't be <= 0!");
            double cdf_one = 0.5*(1.1*createDistanceBond_[i][j]/minrad + 1.1 * stress / (normalBondStiffness_[i][j] * minrad));
            cdf_all = cdf_one > cdf_all ? cdf_one : cdf_all;

          }
        }
      }

      if (heatbreakflag_)
      {
        breakmode_ = BREAKSTYLE_STRESS_TEMP;
        error->all(FLERR,"Bond breakage due to temperature no implemented at the moment!");
      }

      //set neighbor contact_distance_factor here
      if(cdf_all > 10.)
        error->all(FLERR,"Maximum bond distance exceeding 10 x particle diameter, please reduce maxDistanceBond or maxSigmaBond/maxTauBond");
      neighbor->register_contact_dist_factor(cdf_all);
      compute_bond_counter_ = static_cast<ComputeBondCounter*>(modify->find_compute_style_strict("bond/counter", 0));
    }

    /* ---------------------------------------------------------------------- */

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&)
    {
        if(0 == neighbor->ago && fix_bond_random_id_)
            fix_bond_random_id_->do_forward_comm();
    }

    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    /* ---------------------------------------------------------------------- */

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      const bool update_history = scdata.computeflag && scdata.shearupdate;
      // initial settings for first contact are in surfaceTouch

      // check if cohesion model is working for this type-pair
      const int i = scdata.i;
      const int j = scdata.j;
      const int itype = scdata.itype;
      const int jtype = scdata.jtype;
      const double lambda = radiusMultiplier_[itype][jtype];
      if (lambda < SMALL_COHESION_MODEL_BOND)
        return; // do nothing in case of rb = 0.0

      // history values
      double * const contflag = &scdata.contact_history[history_offset_];
      double * const contact_pos = &scdata.contact_history[history_offset_+2]; // 3-elements

      // gather required properties
      // the distance r is the distance between the two particle centers
      // in case of a wall sqrt(scdata.rsq) is the distance to the intersection point on the wall, thus, assuming
      // that the wall is a particle of the same size, the distance of the centers is twice that distance
      double r_create = (scdata.is_wall ? 2.0 : 1.0) * sqrt(scdata.rsq);
      double r = (scdata.is_wall ? r_create/2. : r_create);

      // temporay vectors for intermediate results;
      double tmp1[3], tmp2[3];
      if (update_history)
      {
          // check if particles are already bonded
          if (contflag[0] < SMALL_COHESION_MODEL_BOND)
          {

            const bool create_bond_A =  (createbond_always_flag_ || (update->ntimestep == static_cast<int>(tsCreateBond_))) && (r_create < createDistanceBond_[itype][jtype]);

            if (create_bond_A)
                createBond(scdata, r);
            else if (fix_bond_random_id_)
            {

                const double * const bond_rand_id = fix_bond_random_id_->vector_atom;
                const bool create_bond_B = bond_rand_id && !scdata.is_wall && (update->ntimestep < bond_rand_id[i]) && MathExtraLiggghts::compDouble(bond_rand_id[i],bond_rand_id[j]) && (r_create < createDistanceBond_[itype][jtype]);

                if (create_bond_B)
                    createBond(scdata, r);
		else
                    return;
            }
            else
                return;
          }
          else
          {
              // if already bonded:
              // count bonds that already exist
              if (compute_bond_counter_)
                  compute_bond_counter_->bond_count(scdata.is_wall, scdata.i, scdata.j);
              // update contact_pos for moving walls const double dt = update->dt;
              if (scdata.is_wall)
              {
                  const double dt = update->dt;
                  vectorScalarMult3D(scdata.v_j,dt,tmp1);
                  vectorAdd3D(contact_pos,tmp1,contact_pos);
              }
          }
      }
      else if (contflag[0] < SMALL_COHESION_MODEL_BOND)
      {
          return;
      }

      double * const initialdist = &scdata.contact_history[history_offset_+1];
      double * const force_tang_ptr = &scdata.contact_history[history_offset_+5]; // 3-elements
      double * const torque_normal_ptr = &scdata.contact_history[history_offset_+8]; // 3-elements
      double * const torque_tang_ptr = &scdata.contact_history[history_offset_+11]; // 3-elements

      double force_tang[3], torque_normal[3], torque_tang[3];
      vectorCopy3D(force_tang_ptr, force_tang);
      vectorCopy3D(torque_normal_ptr, torque_normal);
      vectorCopy3D(torque_tang_ptr, torque_tang);

      // gather more properties
      const double radi = scdata.radi;
      const double radj = scdata.radj;
      double delta[3];
      if (scdata.is_wall)
      {
          // at walls we assume that the contact point remains the initial one
          vectorSubtract3D(atom->x[scdata.i],contact_pos,delta);
          r = vectorMag3D(delta);
      }
      else
      {
          vectorCopy3D(scdata.delta,delta);
      }

      // simple bond breakage
      if (breakmode_ == BREAKSTYLE_SIMPLE && r > maxDist_[itype][jtype] && update_history)
      {
        breakBond(scdata);
        return; // bond broken
      }

      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;

      const double rinv = 1./r;
      double en[3];
      vectorScalarMult3D(delta,rinv,en);
      const double * const vi = scdata.v_i;
      const double * const vj = scdata.v_j;
      const double * const omegai = atom->omega[i];
      double wall_omega[3] = {0., 0., 0.};
      double * const omegaj = scdata.is_wall ? wall_omega : atom->omega[j];
      const double kn_pb = normalBondStiffness_[itype][jtype];
      const double kt_pb = tangentialBondStiffness_[itype][jtype];
      const double rb = lambda* (scdata.is_wall ? radi : std::min(radi,radj));
      const double A = M_PI*rb*rb;
      const double J = 0.5*A*rb*rb;
      const double I = 0.5*J;
      const double dt = update->dt;

      const double radsuminv = 1./(radi+radj);
      const double cri = scdata.is_wall ? (r) : (r * radi * radsuminv);
      const double crj = r * radj * radsuminv;

      // calculate relative velocities at the contact point (for contact this is done in the surface_model)

      // relative translational velocity

      double vr[3];
      vectorSubtract3D(vi,vj,vr);

      // normal component

      double vn[3];
      vectorProject3D(vr,en,vn);

      // tangential component

      double vt[3];
      vectorSubtract3D(vr,vn,vt);

      // relative rotational velocity for shear

      double wr[3];
      if (scdata.is_wall)
      {
          vectorScalarMult3D(omegai,0.5,wr);
      }
      else
      {
        // wr = (radi*omegai + radj*omegaj) * rinv;
        vectorScalarMult3D(omegai,radi*radsuminv,tmp1);
        vectorScalarMult3D(omegaj,radj*radsuminv,tmp2);
        vectorAdd3D(tmp1,tmp2,wr);
      }

      // relative velocities for shear
      // v_rel = vi + (x_ic x omegai) - (vj + (-1)*(x_cj x omegaj)) = vi - vj + n_ij x (ri*omegai + rj*omegaj)
      // use only tangential part for velocities
      // cross product introduces only tangential velocity

      double vtr[3];
      vectorCross3D(delta,wr,tmp1);
      vectorAdd3D(vt,tmp1,vtr);

      // relative rotational velocity for torsion and bending

      vectorSubtract3D(omegai,omegaj,wr);

      // normal component

      double wn[3];
      vectorProject3D(wr,en,wn);

      // tangential component

      double wt[3];
      vectorSubtract3D(wr,wn,wt);

      // force calculation
      double nforce[3] = {0.,0.,0.}, dtforce[3] = {0.,0.,0.}, dntorque[3] = {0.,0.,0.}, dttorque[3] = {0.,0.,0.};
      double nforce_damped[3] = {0.,0.,0.}, tforce_damped[3] = {0.,0.,0.}, ntorque_damped[3] = {0.,0.,0.}, ttorque_damped[3] = {0.,0.,0.};

      const double displacement = (initialdist[0] - r);

      const double dissipate = dissipationflag_ ? std::min(dt*fn_dissipation_[itype][jtype],1.0) : 0.0; // dissipation_ already inverted during settings

      if (tensionflag_ || compressionflag_)
      {
        // TODO: distinguish for each type; make it time independent (better?)!
        if (dissipationflag_ && update_history)
        {
          // disspiate first - since explicit
          const double dissipation_delta = (r-initialdist[0])*dissipate;
          initialdist[0] += dissipation_delta; // relax normal spring
        }

        // tension or compression? enabled?
        if ( (tensionflag_ && displacement < -1.e-15) || (compressionflag_ && displacement > 1.e-15) )
        {
          // normal force - tracking distance explicitly
          const double frcmag = kn_pb * A * displacement;
          vectorScalarMult3D(en,frcmag,nforce);

          // no rotation since explicit
          // no history value since explicit

          // damping
          // for positive vn (approaching particles) the damping force has to be repulsive (negative!)
          if (dampingflag_)
          {
            const double minvel = 1e-5* (scdata.is_wall ? radi : fmin(radi,radj)) /dt;
            const double dvX = 0.01*nforce[0]*dt;
            const double dvY = 0.01*nforce[1]*dt;
            const double dvZ = 0.01*nforce[2]*dt;
            const double multiplierX = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vn[0]/fmax(minvel, dvX))) : MathExtraLiggghts::sgn(vn[0]);
            const double multiplierY = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vn[1]/fmax(minvel, dvY))) : MathExtraLiggghts::sgn(vn[1]);
            const double multiplierZ = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vn[2]/fmax(minvel, dvZ))) : MathExtraLiggghts::sgn(vn[2]);
            nforce_damped[0] = nforce[0] - fn_damping_[itype][jtype]*fabs(nforce[0])*multiplierX;
            nforce_damped[1] = nforce[1] - fn_damping_[itype][jtype]*fabs(nforce[1])*multiplierY;
            nforce_damped[2] = nforce[2] - fn_damping_[itype][jtype]*fabs(nforce[2])*multiplierZ;
          } else {
            vectorCopy3D(nforce,nforce_damped);
          }
        }
      }

      // all other forces - track difference
      if (shearflag_) {
        // calc change in shear forces
        vectorScalarMult3D(vtr,-kt_pb*A*dt,dtforce);

        // rotate old forces (only non-explicity calculated one)
        //rotate tangential force
        vectorProject3D(force_tang,en,tmp1);
        vectorSubtract3D(force_tang,tmp1,force_tang);

        // update history value and dissipate
        // TODO: distinguish for each type; make it time independent (better?)!
        if (dissipationflag_)
          vectorScalarMult3D(force_tang, 1.0 - fmin(dt*ft_dissipation_[itype][jtype], 1.0));

        vectorAdd3D(force_tang,dtforce,force_tang);

        // damping
        if (dampingflag_)
        {
            const double minvel = 1e-5* (scdata.is_wall ? radi : fmin(radi,radj)) /dt;
            const double dvX = 0.01*force_tang[0]*dt;
            const double dvY = 0.01*force_tang[1]*dt;
            const double dvZ = 0.01*force_tang[2]*dt;
            const double multiplierX = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vtr[0]/fmax(minvel, dvX))) : MathExtraLiggghts::sgn(vtr[0]);
            const double multiplierY = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vtr[1]/fmax(minvel, dvY))) : MathExtraLiggghts::sgn(vtr[1]);
            const double multiplierZ = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vtr[2]/fmax(minvel, dvZ))) : MathExtraLiggghts::sgn(vtr[2]);
            tforce_damped[0] = force_tang[0] - ft_damping_[itype][jtype]*fabs(force_tang[0])*multiplierX;
            tforce_damped[1] = force_tang[1] - ft_damping_[itype][jtype]*fabs(force_tang[1])*multiplierY;
            tforce_damped[2] = force_tang[2] - ft_damping_[itype][jtype]*fabs(force_tang[2])*multiplierZ;
        } else {
          vectorCopy3D(force_tang,tforce_damped);
        }
      }

      if (ntorqueflag_) {
        // calc change in normal torque
        vectorScalarMult3D(wn,-kt_pb*J*dt,dntorque);

        //rotate normal torque
        vectorProject3D(torque_normal,en,torque_normal);

        // update history value and dissipate
        // TODO: distinguish for each type; make it time independent (better?)!
        if (dissipationflag_)
          vectorScalarMult3D(torque_normal, 1.0 - fmin(dt*tn_dissipation_[itype][jtype], 1.0));

        vectorAdd3D(torque_normal,dntorque,torque_normal);

        // damping
        if (dampingflag_)
        {
          ntorque_damped[0] = torque_normal[0] - tn_damping_[itype][jtype]*fabs(torque_normal[0])*MathExtraLiggghts::sgn(wn[0]);
          ntorque_damped[1] = torque_normal[1] - tn_damping_[itype][jtype]*fabs(torque_normal[1])*MathExtraLiggghts::sgn(wn[1]);
          ntorque_damped[2] = torque_normal[2] - tn_damping_[itype][jtype]*fabs(torque_normal[2])*MathExtraLiggghts::sgn(wn[2]);
        } else {
          vectorCopy3D(torque_normal,ntorque_damped);
        }
      }

      if (ttorqueflag_)
      {
        const double wtsq = vectorMag3DSquared(wt);

        // the current projection works only if the tangential relative rotational velocity
        // is unequal zero!
        // thus the rotation of the torque vector is only possible in case of wtsq != 0.
        // TODO: Check this. Maybe there is a better solution?
        if (wtsq > 0)
        {
          // calc change in tang torque
          vectorScalarMult3D(wt,-kn_pb*I*dt,dttorque);

          //rotate tangential torque (project onto the relative tangential rotation)
          // TODO: Better version? Maybe rotation matrix?
          vectorProject3D(torque_tang,wt,torque_tang);

          // update history value and dissipate
          // TODO: distinguish for each type; make it time independent (better?)!
          if (dissipationflag_)
            vectorScalarMult3D(torque_tang, 1.0 - fmin(dt*tt_dissipation_[itype][jtype], 1.0));

          vectorAdd3D(torque_tang,dttorque,torque_tang);

          // damping
          if (dampingflag_)
          {
            ttorque_damped[0] = torque_tang[0] - tt_damping_[itype][jtype]*fabs(torque_tang[0])*MathExtraLiggghts::sgn(wt[0]);
            ttorque_damped[1] = torque_tang[1] - tt_damping_[itype][jtype]*fabs(torque_tang[1])*MathExtraLiggghts::sgn(wt[1]);
            ttorque_damped[2] = torque_tang[2] - tt_damping_[itype][jtype]*fabs(torque_tang[2])*MathExtraLiggghts::sgn(wt[2]);
          } else {
            vectorCopy3D(torque_tang,ttorque_damped);
          }
        }
      }

      // bond breakage due to stress (use un-damped forces/torques)
      if (breakmode_ == BREAKSTYLE_STRESS)
      {
        const double nforce_mag = vectorMag3D(nforce); //TODO: Potyondy and Cundall do not use the absolut value
        const double tforce_mag = vectorMag3D(force_tang);
        const double ntorque_mag = vectorMag3D(torque_normal);
        const double ttorque_mag = vectorMag3D(torque_tang);

        double s_factor = 1.;
        if (fix_strenghtening_factor_)
        {

           double **s_f =  fix_strenghtening_factor_->array_atom;

           if (scdata.is_wall && !MathExtraLiggghts::compDouble(s_f[i][0],0.))
            s_factor = s_f[i][0];
           else if(!MathExtraLiggghts::compDouble(s_f[i][0],0.) && !MathExtraLiggghts::compDouble(s_f[j][0],0.))
            s_factor = 2. * (s_f[i][0] * s_f[j][0]) / (s_f[i][0] + s_f[j][0]);

        }

        const double maxSigma = maxSigma_[itype][jtype] * ((useRatioTC_ && displacement < 1e-16) ? ratioTensionCompression_[itype][jtype] : 1.0);
        const bool nstress = ( s_factor * maxSigma ) < (nforce_mag/A + ttorque_mag*rb/I);

        const double maxTau = druckerPrager_ ?
            maxTau_[itype][jtype] + MathExtraLiggghts::sgn(displacement)*nforce_mag/A*tanFrictionAngle_[itype][jtype] :
            maxTau_[itype][jtype];
        const bool tstress = ( s_factor * maxTau ) < (tforce_mag/A + ntorque_mag*rb/J);

        if ((nstress || tstress) && update_history) // bond broken?
        {
          breakBond(scdata);
          return; // bond broken
        }
      }

      // calculate current torque: -(r x Ft) = Ft x r
      double tor[3] = {0., 0., 0.};
      vectorCross3D(tforce_damped,en,tor);

      // inform pair style to update forces!
      scdata.has_force_update = true;

      const double F[3] = {nforce_damped[0] + tforce_damped[0],
                           nforce_damped[1] + tforce_damped[1],
                           nforce_damped[2] + tforce_damped[2]};
      const double dTorque_i[3] = {cri*tor[0]+ntorque_damped[0]+ttorque_damped[0],
                                   cri*tor[1]+ntorque_damped[1]+ttorque_damped[1],
                                   cri*tor[2]+ntorque_damped[2]+ttorque_damped[2]};
      const double dTorque_j[3] = {crj*tor[0]-ntorque_damped[0]-ttorque_damped[0],
                                   crj*tor[1]-ntorque_damped[1]-ttorque_damped[1],
                                   crj*tor[2]-ntorque_damped[2]-ttorque_damped[2]};

      if (update_history)
      {
        vectorCopy3D(force_tang, force_tang_ptr);
        vectorCopy3D(torque_normal, torque_normal_ptr);
        vectorCopy3D(torque_tang, torque_tang_ptr);

        if (elasticpotflag_ || dissipatedflag_)
        {
            // elastic torque component
            double torEl[3];
            vectorCross3D(force_tang,en,torEl);
            const double dTorqueEl_i[3] = {cri*torEl[0]+torque_normal[0]+torque_tang[0],
                                           cri*torEl[1]+torque_normal[1]+torque_tang[1],
                                           cri*torEl[2]+torque_normal[2]+torque_tang[2]};
            const double dTorqueEl_j[3] = {crj*torEl[0]+torque_normal[0]+torque_tang[0],
                                           crj*torEl[1]+torque_normal[1]+torque_tang[1],
                                           crj*torEl[2]+torque_normal[2]+torque_tang[2]};
            // compute increment in elastic potential
            double * const elastic_pot = &scdata.contact_history[elastic_potential_offset_];
            const double dissipated_energy = dissipationflag_ ? (elastic_pot[0])*dissipate*dissipate : 0.0; // square as E_elastic = k x^2
            if (elasticpotflag_)
            {
                if (scdata.is_wall)
                {
                    double delta[3] = {0.,0.,0.};
                    vectorCopy3D(scdata.v_j,delta);
                    vectorScalarMult3D(delta, update->dt);
                    // -= because force is in opposite direction
                    // no *dt as delta is v*dt of the contact position

                    elastic_pot[0] -= (delta[0]*elastic_pot[1] +
                                       delta[1]*elastic_pot[2] +
                                       delta[2]*elastic_pot[3])*0.5
                                       // from previous half step
                                       + elastic_pot[10];

                    elastic_pot[10] = -(delta[0]*(nforce[0] + force_tang[0]) +
                                        delta[1]*(nforce[1] + force_tang[1]) +
                                        delta[2]*(nforce[2] + force_tang[2]))*0.5;
                }
                elastic_pot[1] = -(nforce[0] + force_tang[0]);
                elastic_pot[2] = -(nforce[1] + force_tang[1]);
                elastic_pot[3] = -(nforce[2] + force_tang[2]);
                elastic_pot[4] = -dTorqueEl_i[0];
                elastic_pot[5] = -dTorqueEl_i[1];
                elastic_pot[6] = -dTorqueEl_i[2];
                elastic_pot[7] = -dTorqueEl_j[0];
                elastic_pot[8] = -dTorqueEl_j[1];
                elastic_pot[9] = -dTorqueEl_j[2];
                if (dissipationflag_)
                    *elastic_pot -= dissipated_energy; // remove dissipated energy due to dissipation model
            }
            if (dissipatedflag_)
            {
                double * const * const dissipated = fix_dissipated_->array_atom;

                double * const dissipated_i = dissipated[scdata.i];
                double * const dissipated_j = dissipated[scdata.j];
                if (force->newton_pair || scdata.j < atom->nlocal || atom->tag[scdata.i] < atom->tag[scdata.j])
                    dissipated_i[0] += dissipated_energy;
                dissipated_i[1] += -(F[0] - (nforce[0] + force_tang[0]));
                dissipated_i[2] += -(F[1] - (nforce[1] + force_tang[1]));
                dissipated_i[3] += -(F[2] - (nforce[2] + force_tang[2]));
                const double dTorqueDamp_i[3] = {dTorque_i[0]-dTorqueEl_i[0], dTorque_i[1]-dTorqueEl_i[1], dTorque_i[2]-dTorqueEl_i[2]};
                dissipated_i[4] += -dTorqueDamp_i[0];
                dissipated_i[5] += -dTorqueDamp_i[1];
                dissipated_i[6] += -dTorqueDamp_i[2];
                if (scdata.j < atom->nlocal && !scdata.is_wall)
                {
                    dissipated_j[1] -= -(F[0] - (nforce[0] + force_tang[0]));
                    dissipated_j[2] -= -(F[1] - (nforce[1] + force_tang[1]));
                    dissipated_j[3] -= -(F[2] - (nforce[2] + force_tang[2]));
                    const double dTorqueDamp_j[3] = {dTorque_j[0]-dTorqueEl_j[0], dTorque_j[1]-dTorqueEl_j[1], dTorque_j[2]-dTorqueEl_j[2]};
                    dissipated_j[4] += -dTorqueDamp_j[0];
                    dissipated_j[5] += -dTorqueDamp_j[1];
                    dissipated_j[6] += -dTorqueDamp_j[2];
                }
                else if (scdata.is_wall)
                {
                    double * const diss_force = &scdata.contact_history[dissipation_history_offset_];
                    diss_force[0] += F[0] - (nforce[0] + force_tang[0]);
                    diss_force[1] += F[1] - (nforce[1] + force_tang[1]);
                    diss_force[2] += F[2] - (nforce[2] + force_tang[2]);
                }
            }
        }
      }

      // return resulting forces
      if(scdata.is_wall) {

        i_forces.delta_F[0] += F[0];
        i_forces.delta_F[1] += F[1];
        i_forces.delta_F[2] += F[2];
        i_forces.delta_torque[0] += dTorque_i[0];
        i_forces.delta_torque[1] += dTorque_i[1];
        i_forces.delta_torque[2] += dTorque_i[2];
      } else {
        i_forces.delta_F[0] += F[0];
        i_forces.delta_F[1] += F[1];
        i_forces.delta_F[2] += F[2];
        i_forces.delta_torque[0] += dTorque_i[0];
        i_forces.delta_torque[1] += dTorque_i[1];
        i_forces.delta_torque[2] += dTorque_i[2];

        j_forces.delta_F[0] -= F[0];
        j_forces.delta_F[1] -= F[1];
        j_forces.delta_F[2] -= F[2];
        j_forces.delta_torque[0] += dTorque_j[0];
        j_forces.delta_torque[1] += dTorque_j[1];
        j_forces.delta_torque[2] += dTorque_j[2];
      }
    }

    /* ---------------------------------------------------------------------- */

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      // check if cohesion model is working for this type-pair
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double lambda = radiusMultiplier_[itype][jtype];
      if (lambda < SMALL_COHESION_MODEL_BOND)
        return; // do nothing in case of rb=0.0

      //initial settings for first contact

      // the distance r is the distance between the two particle centers
      // in case of a wall sidata.r is the distance to the intersection point on the wall, thus, assuming
      // that the wall is a particle of the same size, the distance of the centers is twice that distance
      //  const double r = (sidata.is_wall ? 2.0 : 1.0) * sidata.r;
      //create bond
      // if createbondflag_ && correct timestep && within range

      //if( (update->ntimestep == static_cast<int>(tsCreateBond_) && createbondflag_ && r < createDistanceBond_[itype][jtype]) )
      //  createBond(sidata, r);

      // TODO create member function that does the calculation and gets the velocities as input parameter.
      //      thus we do not double calculate velocities in case of contact!
      surfacesClose(sidata, i_forces, j_forces);
    }

    void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}

    /* ---------------------------------------------------------------------- */

  private:
    void createBond(SurfacesCloseData & scdata, double dist)
    {
      double * const contflag = &scdata.contact_history[history_offset_];
      double * const initialdist = &scdata.contact_history[history_offset_+1];
      double * const contact_pos = &scdata.contact_history[history_offset_+2]; // 3-elements
      double * const ft_old = &scdata.contact_history[history_offset_+5]; // 3-elements
      double * const torque_normal = &scdata.contact_history[history_offset_+8]; // 3-elements
      double * const torque_tang = &scdata.contact_history[history_offset_+11]; // 3-elements

      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
      contflag[0] = 1.0;
      initialdist[0] = dist; // use current distance
      vectorSubtract3D(atom->x[scdata.i], scdata.delta, contact_pos);
      vectorZeroize3D(ft_old);
      vectorZeroize3D(torque_normal);
      vectorZeroize3D(torque_tang);

      if (compute_bond_counter_)
          compute_bond_counter_->bond_created(scdata.is_wall, scdata.i, scdata.j);
    }

    void breakBond(SurfacesCloseData & scdata)
    {
      double * const contflag = &scdata.contact_history[history_offset_];
      double * const initialdist = &scdata.contact_history[history_offset_+1];
      double * const elastic_pot = &scdata.contact_history[elastic_potential_offset_];

      if (scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
      contflag[0] = 0.;
      initialdist[0] = 0.;
      if (elasticpotflag_)
      {
        if (dissipatedflag_ && (scdata.is_wall || force->newton_pair || scdata.j < atom->nlocal || atom->tag[scdata.i] < atom->tag[scdata.j]))
        {
            // elastic energy is dissipated instantly
            double * const * const dissipated = fix_dissipated_->array_atom;
            double * const dissipated_i = dissipated[scdata.i];
            dissipated_i[0] += elastic_pot[0];
        }
        // reset elastic energy, forces & torque
        elastic_pot[0] = 0.0;
        elastic_pot[1] = 0.0;
        elastic_pot[2] = 0.0;
        elastic_pot[3] = 0.0;
        elastic_pot[4] = 0.0;
        elastic_pot[5] = 0.0;
        elastic_pot[6] = 0.0;
        elastic_pot[7] = 0.0;
        elastic_pot[8] = 0.0;
        elastic_pot[9] = 0.0;
      }
      if (compute_bond_counter_)
          compute_bond_counter_->bond_broken(scdata.is_wall, scdata.i, scdata.j);
    }

    /* ---------------------------------------------------------------------- */

  private:
    int history_offset_;
    int elastic_potential_offset_;
    int dissipation_history_offset_;
    class FixPropertyAtom *fix_dissipated_;
    double ** radiusMultiplier_;
    double ** normalBondStiffness_;
    double ** tangentialBondStiffness_;
    double ** fn_dissipation_;
    double ** ft_dissipation_;
    double ** tn_dissipation_;
    double ** tt_dissipation_;

    double ** fn_damping_;
    double ** ft_damping_;
    double ** tn_damping_;
    double ** tt_damping_;

    double ** maxDist_;
    double ** maxSigma_;
    double ** maxTau_;
    double ** ratioTensionCompression_;
    double ** tanFrictionAngle_;

    double tsCreateBond_;
    double ** createDistanceBond_;

  protected:
    bool stressbreakflag_; // false if max dist break, true if stress break
    bool heatbreakflag_;

    bool tensionflag_;
    bool compressionflag_;
    bool shearflag_;
    bool ntorqueflag_;
    bool ttorqueflag_;

    bool createbond_always_flag_;

    bool dampingflag_;
    bool dampingsmooth_;
    bool dissipationflag_;

    bool elasticpotflag_;
    bool dissipatedflag_;

    int breakmode_;
    bool useRatioTC_;
    bool druckerPrager_;

    // fix holding per-atom values for optional strengthening
    // of the bond
    class FixPropertyAtom *fix_strenghtening_factor_;
    class FixPropertyAtom *fix_bond_random_id_;

    // compute that keeps track of the number of bonds and created/destroyed ones
    class ComputeBondCounter *compute_bond_counter_;
  };
}
}
#endif // COHESION_MODEL_BOND_H_
#endif
