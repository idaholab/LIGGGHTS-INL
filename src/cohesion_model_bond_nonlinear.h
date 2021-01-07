/* ----------------------------------------------------------------------
   Copyright 2021, Battelle Energy Alliance, LLC  All Rights Reserved
-------------------------------------------------------------------------

   Contributing author and copyright for this file:

   Yidong Xia (Idaho National Laboratory)
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

    Yidong Xia (Idaho National Laboratory; yidong.xia@inl.gov)
    Yuan Guo (Clemson University & Idaho National Laboratory)
    Qiushi Chen (Clemson University)

    v1.0: Separate stiffness in each freedom, linear model
    v1.1: Quadratic model for compression, linear for tension
    v1.2: Squre root model for tension, linear for compression
    v1.3: Squre root model for compression, linear for tension
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_BOND_NONLINEAR,bond/nonlinear,15)
#else

#ifndef COHESION_MODEL_BOND_NONLINEAR_H_
#define COHESION_MODEL_BOND_NONLINEAR_H_

#include "contact_models.h"
#include "cohesion_model_base.h"
#include "compute_bond_counter.h"
#include <cmath>
#include <algorithm>
#include "error.h"
#include "math_extra_liggghts.h"
#include "fix_property_atom.h"
#include "neighbor.h"

// required model properties
namespace MODEL_PARAMS
{
  inline static MatrixProperty* createRadiusMultiplierBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "radiusMultiplierBondnonlinear", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createDampingNormalForceBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "dampingNormalForceBondnonlinear", caller, sanity_checks, 0.0, 1.0);}

  inline static MatrixProperty* createDampingTangentialForceBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "dampingTangentialForceBondnonlinear", caller, sanity_checks, 0.0, 1.0);}

  inline static MatrixProperty* createDampingNormalTorqueBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "dampingNormalTorqueBondnonlinear", caller, sanity_checks, 0.0, 1.0);}

  inline static MatrixProperty* createDampingTangentialTorqueBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "dampingTangentialTorqueBondnonlinear", caller, sanity_checks, 0.0, 1.0);}

  inline static MatrixProperty* createMaxDistanceBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "maxDistanceBondnonlinear", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createMaxSigmaBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "maxSigmaBondnonlinear", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createMaxTauBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "maxTauBondnonlinear", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createRatioTensionCompressionBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "ratioTensionCompressionBondnonlinear", caller, sanity_checks, 0.0);}

  inline static ScalarProperty* createTsCreateBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createScalarProperty(registry, "tsCreateBondnonlinear", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createCreateDistanceBondnonlinear(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "createDistanceBondnonlinear", caller, sanity_checks, 0.0);}

  // bond stiffness
  inline static MatrixProperty* createStiffnessK_fn1(PropertyRegistry & registry, const char * caller, bool sanity_checks)           // compression
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaK_fn1", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKu_fn1(PropertyRegistry & registry, const char * caller, bool sanity_checks)          // compression
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKu_fn1", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKc_fn1(PropertyRegistry & registry, const char * caller, bool sanity_checks)          // compression
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKc_fn1", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessK_fn2(PropertyRegistry & registry, const char * caller, bool sanity_checks)           // tension
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaK_fn2", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKu_fn2(PropertyRegistry & registry, const char * caller, bool sanity_checks)          // tension
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKu_fn2", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKc_fn2(PropertyRegistry & registry, const char * caller, bool sanity_checks)          // tension
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKc_fn2", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessK_ft(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaK_ft", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessK_tn(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaK_tn", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKu_tn(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKu_tn", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKc_tn(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKc_tn", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessK_tt(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaK_tt", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKu_tt(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKu_tt", caller, sanity_checks, 0.0);}

  inline static MatrixProperty* createStiffnessKc_tt(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {return createPerTypePairProperty(registry, "stiffnessPerUnitAreaKc_tt", caller, sanity_checks, 0.0);}
}

/* ---------------------------------------------------------------------- */
namespace LIGGGHTS {

namespace ContactModels {
  enum{
       BREAK_SIMPLE,
       BREAK_STRESS
      };

  template<>
  class CohesionModel<COHESION_BOND_NONLINEAR> : public CohesionModelBase {

  public:

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
      CohesionModelBase(lmp, hsetup, cmb),
      history_offset_(0),
      radiusMultiplierBondnonlinear_(NULL),
      stiffnessK_fn1_(NULL),
      stiffnessKu_fn1_(NULL),
      stiffnessKc_fn1_(NULL),
      stiffnessK_fn2_(NULL),
      stiffnessKu_fn2_(NULL),
      stiffnessKc_fn2_(NULL),
      stiffnessK_ft_(NULL),
      stiffnessK_tn_(NULL),
      stiffnessKu_tn_(NULL),
      stiffnessKc_tn_(NULL),
      stiffnessK_tt_(NULL),
      stiffnessKu_tt_(NULL),
      stiffnessKc_tt_(NULL),
      maxDist_(NULL),
      maxSigma_(NULL),
      maxTau_(NULL),
      ratioTensionCompressionBondnonlinear_(NULL),                                    // tensile/compression strength (positive value), sigma_max is compression strength
      tsCreateBondnonlinear_(0.0),                                                    // e.g., 0.5 corresponds to tensile strength equaling to half of compression strength
      createDistanceBondnonlinear_(NULL),
      stressbreakflag_(false),
      tensionflag_(true),
      compressionflag_(true),
      shearflag_(true),
      ntorqueflag_(true),
      ttorqueflag_(true),
      createBondNonlinear_always_flag_(false),
      dampingflag_(true),
      dampingsmooth_(false),
      breakmode_(BREAK_SIMPLE),
      useRatioTC_(false),
      fix_bondnonlinear_random_id_(NULL),
      compute_bond_counter_(NULL)
    {
      // by adding a history value with add_history_value, the second argument is the newtonflag
      // if newtonflag = 1, the history value changes its sign depending on particle interaction (i->j or j->i)
      // if newtonflag = 0, it does nothing
      history_offset_ = hsetup->add_history_value("bondFlag", "0");               // history_offset_
      hsetup->add_history_value("initial_dist", "0");                             // history_offset_+1
      hsetup->add_history_value("contactPosX", "0");                              // history_offset_+2
      hsetup->add_history_value("contactPosY", "0");
      hsetup->add_history_value("contactPosZ", "0");
      hsetup->add_history_value("ftx", "1");                                      // history_offset_+5
      hsetup->add_history_value("fty", "1");
      hsetup->add_history_value("ftz", "1");
      hsetup->add_history_value("thetanx", "1");                                  // history_offset_+8
      hsetup->add_history_value("thetany", "1");
      hsetup->add_history_value("thetanz", "1");
      hsetup->add_history_value("thetatx", "1");                                  // history_offset_+11
      hsetup->add_history_value("thetaty", "1");
      hsetup->add_history_value("thetatz", "1");
      hsetup->add_history_value("deltan_max", "0");                               // history_offset_+14, deltan_max
      hsetup->add_history_value("thetanx_max", "0");                              // history_offset_+15
      hsetup->add_history_value("thetany_max", "0");
      hsetup->add_history_value("thetanz_max", "0");
      hsetup->add_history_value("thetatx_max", "0");                              // history_offset_+18
      hsetup->add_history_value("thetaty_max", "0");
      hsetup->add_history_value("thetatz_max", "0");
      hsetup->add_history_value("thetanx_min", "0");                              // history_offset_+21
      hsetup->add_history_value("thetany_min", "0");
      hsetup->add_history_value("thetanz_min", "0");
      hsetup->add_history_value("thetatx_min", "0");                              // history_offset_+24
      hsetup->add_history_value("thetaty_min", "0");
      hsetup->add_history_value("thetatz_min", "0");
      hsetup->add_history_value("deltan_min", "0");                               // history_offset_+27, deltan_min
      cmb->add_history_offset("bondnonlinear_contactflag", history_offset_);
    }

    /* ---------------------------------------------------------------------- */
    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("stressBreak", stressbreakflag_, false);
      settings.registerOnOff("tensionStress", tensionflag_, true);
      settings.registerOnOff("compressionStress", compressionflag_, true);
      settings.registerOnOff("shearStress", shearflag_, true);
      settings.registerOnOff("normalTorqueStress", ntorqueflag_, true);
      settings.registerOnOff("shearTorqueStress", ttorqueflag_, true);
      settings.registerOnOff("createBondAlways", createBondNonlinear_always_flag_, false); // if true, bonds are created if within range, else at a certain timestep and within range
      settings.registerOnOff("dampingBond", dampingflag_, true);
      settings.registerOnOff("dampingBondSmooth", dampingsmooth_, false);
      settings.registerOnOff("ratioTensionCompressionBond", useRatioTC_, false);
    }

    /* ---------------------------------------------------------------------- */
    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb){}

    /* ---------------------------------------------------------------------- */
    void connectToProperties(PropertyRegistry & registry)
    {
      // check if fix_bondnonlinear_random_id exist
      fix_bondnonlinear_random_id_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("bondnonlinear_random_id","property/atom","scalar",0,0,"cohesion bond/nonlinear",false));

      registry.registerProperty("radiusMultiplierBondnonlinear", &MODEL_PARAMS::createRadiusMultiplierBondnonlinear);
      registry.registerProperty("stiffnessK_fn1", &MODEL_PARAMS::createStiffnessK_fn1);
      registry.registerProperty("stiffnessKu_fn1", &MODEL_PARAMS::createStiffnessKu_fn1);
      registry.registerProperty("stiffnessKc_fn1", &MODEL_PARAMS::createStiffnessKc_fn1);
      registry.registerProperty("stiffnessK_fn2", &MODEL_PARAMS::createStiffnessK_fn2);
      registry.registerProperty("stiffnessKu_fn2", &MODEL_PARAMS::createStiffnessKu_fn2);
      registry.registerProperty("stiffnessKc_fn2", &MODEL_PARAMS::createStiffnessKc_fn2);
      registry.registerProperty("stiffnessK_ft", &MODEL_PARAMS::createStiffnessK_ft);
      registry.registerProperty("stiffnessK_tn", &MODEL_PARAMS::createStiffnessK_tn);
      registry.registerProperty("stiffnessKu_tn", &MODEL_PARAMS::createStiffnessKu_tn);
      registry.registerProperty("stiffnessKc_tn", &MODEL_PARAMS::createStiffnessKc_tn);
      registry.registerProperty("stiffnessK_tt", &MODEL_PARAMS::createStiffnessK_tt);
      registry.registerProperty("stiffnessKu_tt", &MODEL_PARAMS::createStiffnessKu_tt);
      registry.registerProperty("stiffnessKc_tt", &MODEL_PARAMS::createStiffnessKc_tt);

      registry.connect("radiusMultiplierBondnonlinear", radiusMultiplierBondnonlinear_,"cohesion bond/nonlinear");
      registry.connect("stiffnessK_fn1", stiffnessK_fn1_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKu_fn1", stiffnessKu_fn1_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKc_fn1", stiffnessKc_fn1_,"cohesion bond/nonlinear");
      registry.connect("stiffnessK_fn2", stiffnessK_fn2_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKu_fn2", stiffnessKu_fn2_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKc_fn2", stiffnessKc_fn2_,"cohesion bond/nonlinear");
      registry.connect("stiffnessK_ft", stiffnessK_ft_,"cohesion bond/nonlinear");
      registry.connect("stiffnessK_tn", stiffnessK_tn_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKu_tn", stiffnessKu_tn_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKc_tn", stiffnessKc_tn_,"cohesion bond/nonlinear");
      registry.connect("stiffnessK_tt", stiffnessK_tt_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKu_tt", stiffnessKu_tt_,"cohesion bond/nonlinear");
      registry.connect("stiffnessKc_tt", stiffnessKc_tt_,"cohesion bond/nonlinear");

      // create damping settings
      if (dampingsmooth_) dampingflag_ = true;

      if (dampingflag_)
      {
        registry.registerProperty("dampingNormalForceBondnonlinear", &MODEL_PARAMS::createDampingNormalForceBondnonlinear);
        registry.registerProperty("dampingTangentialForceBondnonlinear", &MODEL_PARAMS::createDampingTangentialForceBondnonlinear);
        registry.registerProperty("dampingNormalTorqueBondnonlinear", &MODEL_PARAMS::createDampingNormalTorqueBondnonlinear);
        registry.registerProperty("dampingTangentialTorqueBondnonlinear", &MODEL_PARAMS::createDampingTangentialTorqueBondnonlinear);

        registry.connect("dampingNormalForceBondnonlinear", damping_fn_,"cohesion bond/nonlinear");
        registry.connect("dampingTangentialForceBondnonlinear", damping_ft_,"cohesion bond/nonlinear");
        registry.connect("dampingNormalTorqueBondnonlinear", damping_tn_,"cohesion bond/nonlinear");
        registry.connect("dampingTangentialTorqueBondnonlinear", damping_tt_,"cohesion bond/nonlinear");
      }

      // create contact settings
      if (createBondNonlinear_always_flag_)
      {
        registry.registerProperty("createDistanceBondnonlinear", &MODEL_PARAMS::createCreateDistanceBondnonlinear);
        registry.connect("createDistanceBondnonlinear", createDistanceBondnonlinear_, "cohesion bond/nonlinear");
      }
      else
      {
        registry.registerProperty("tsCreateBondnonlinear", &MODEL_PARAMS::createTsCreateBondnonlinear);
        registry.registerProperty("createDistanceBondnonlinear", &MODEL_PARAMS::createCreateDistanceBondnonlinear);

        registry.connect("tsCreateBondnonlinear", tsCreateBondnonlinear_, "cohesion bond/nonlinear");
        registry.connect("createDistanceBondnonlinear", createDistanceBondnonlinear_, "cohesion bond/nonlinear");
      }

      // break up settings
      // find cdf: contact distance factor
      double cdf_all = 1.;                                                        // minimum value is one
      const int max_type = registry.max_type();
      const double minrad = registry.min_radius();                                // min radius of all spheres

      if (minrad <= 0.)
        error->one(FLERR,"Bondnonlinear settings: The minimum radius can't be <= 0!");

      if (!stressbreakflag_)
      {
        // required settings
        registry.registerProperty("maxDistanceBondnonlinear", &MODEL_PARAMS::createMaxDistanceBondnonlinear);
        registry.connect("maxDistanceBondnonlinear", maxDist_,"cohesion bond/nonlinear");
        breakmode_ = BREAK_SIMPLE;

        // calculate the contactDistanceFactor
        for (int i=1; i<max_type+1; i++)
        {
          for(int j=1;j<max_type+1;j++)
          {
            double cdf_one = 1.1 * 0.5*maxDist_[i][j]/minrad;
            cdf_all = cdf_one > cdf_all ? cdf_one : cdf_all;                      // whichever is larger
          }
        }
      }
      else
      {
        // required settings
        registry.registerProperty("maxSigmaBondnonlinear", &MODEL_PARAMS::createMaxSigmaBondnonlinear);
        registry.registerProperty("maxTauBondnonlinear", &MODEL_PARAMS::createMaxTauBondnonlinear);
        registry.connect("maxSigmaBondnonlinear", maxSigma_,"cohesion bond/nonlinear");
        registry.connect("maxTauBondnonlinear", maxTau_,"cohesion bond/nonlinear");
        if (useRatioTC_)
        {
          registry.registerProperty("ratioTensionCompressionBondnonlinear", &MODEL_PARAMS::createRatioTensionCompressionBondnonlinear);
          registry.connect("ratioTensionCompressionBondnonlinear", ratioTensionCompressionBondnonlinear_,"cohesion bond/nonlinear");
        }
        breakmode_ = BREAK_STRESS;

        // calculate the cdf
        // it depends only on the normalStress due to normal distance, maxSigma_, stiffnessK_fn2_, minrad!
        for (int i=1; i<max_type+1; i++)
        {
          for (int j=1; j<max_type+1; j++)
          {
            double stress = maxSigma_[i][j];
            if (useRatioTC_)
                stress = fmax(stress, stress*ratioTensionCompressionBondnonlinear_[i][j]);
            if(stiffnessK_fn2_[i][j] <= 1e-15)
                error->one(FLERR,"Bondnonlinear settings: In case of stress breakage, the normal bond/nonlinear stiffness can't be <= 0!");
            double cdf_one = 0.5*1.1*createDistanceBondnonlinear_[i][j]/minrad + 0.5*1.1*stress/(stiffnessK_fn2_[i][j]*minrad); // Linear model for tension
            cdf_all = cdf_one > cdf_all ? cdf_one : cdf_all;
          }
        }
      }

      //set neighbor contact_distance_factor here
      if(cdf_all > 10.)
        error->all(FLERR,"Maximum bondnonlinear distance exceeding 10 x particle diameter, please reduce maxDistanceBondnonlinear or maxSigmaBondnonlinear/maxTauBondnonlinear");
      neighbor->register_contact_dist_factor(cdf_all);
      compute_bond_counter_ = static_cast<ComputeBondCounter*>(modify->find_compute_style_strict("bond/nonlinear/counter", 0));
    }

    /* ---------------------------------------------------------------------- */
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&)
    {
        if(0 == neighbor->ago && fix_bondnonlinear_random_id_)
            fix_bondnonlinear_random_id_->do_forward_comm();
    }

    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    /* ---------------------------------------------------------------------- */
    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      const bool update_history = scdata.computeflag && scdata.shearupdate;
      // initial settings for first contact are in surfaceTouch

      // gather particle properties
      const int i = scdata.i;
      const int j = scdata.j;
      const int itype = scdata.itype;
      const int jtype = scdata.jtype;
      const double radi = scdata.radi;
      const double radj = scdata.radj;
      const double lambda = radiusMultiplierBondnonlinear_[itype][jtype];
      if (lambda < 1.e-15)
        return; // do nothing in case of rb = 0.0

      // history values
      double * const bondFlag = &scdata.contact_history[history_offset_];
      double * const initialdist = &scdata.contact_history[history_offset_+1];
      double * const contact_pos = &scdata.contact_history[history_offset_+2];       // 3-elements
      double * const force_tang_ptr = &scdata.contact_history[history_offset_+5];    // 3-elements
      double * const theta_normal_ptr = &scdata.contact_history[history_offset_+8];  // 3-elements
      double * const theta_tang_ptr = &scdata.contact_history[history_offset_+11];   // 3-elements

      // the distance r is the distance between the two particle centers
      // in case of a wall sqrt(scdata.rsq) is the distance to the intersection point on the wall, thus, assuming
      // that the wall is a particle of the same size, the distance of the centers is twice that distance
      double r_centers = (scdata.is_wall ? 2.0 : 1.0) * sqrt(scdata.rsq);
      double r = (scdata.is_wall ? r_centers/2. : r_centers);

      // temporay vectors for intermediate results
      double tmp1[3], tmp2[3];
      if (update_history)
      {
          // check if particles are already bonded
          if (bondFlag[0] < 1.e-15)
          //if unbonded:
          {
            const bool create_bondnonlinear_A =  (createBondNonlinear_always_flag_ || (update->ntimestep == static_cast<int>(tsCreateBondnonlinear_))) && (r_centers < createDistanceBondnonlinear_[itype][jtype]);
            if (create_bondnonlinear_A)
                createBondnonlinear(scdata, r);
            else if (fix_bondnonlinear_random_id_)
            {
                const double * const bondnonlinear_rand_id = fix_bondnonlinear_random_id_->vector_atom;
                const bool create_bondnonlinear_B = bondnonlinear_rand_id && !scdata.is_wall && (update->ntimestep < bondnonlinear_rand_id[i]) && MathExtraLiggghts::compDouble(bondnonlinear_rand_id[i],bondnonlinear_rand_id[j]) && (r_centers < createDistanceBondnonlinear_[itype][jtype]);
                if (create_bondnonlinear_B)
                    createBondnonlinear(scdata, r);
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
      else if (bondFlag[0] < 1.e-15)
      {
          return;
      }

      double delta[3], force_tang[3], theta_normal[3], theta_tang[3];
      vectorCopy3D(force_tang_ptr, force_tang);
      vectorCopy3D(theta_normal_ptr, theta_normal);
      vectorCopy3D(theta_tang_ptr, theta_tang);

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
      if (breakmode_ == BREAK_SIMPLE && r > maxDist_[itype][jtype] && update_history)
      {
        breakBondnonlinear(scdata);
        return;                                                                   // bond broken
      }

      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;

      // gather more particle properties
      const double rinv = 1./r;
      double en[3];                                                               // normal force unit vector
      vectorScalarMult3D(delta,rinv,en);
      const double * const vi = scdata.v_i;
      const double * const vj = scdata.v_j;
      const double * const omegai = atom->omega[i];
      double wall_omega[3] = {0., 0., 0.};
      double * const omegaj = scdata.is_wall ? wall_omega : atom->omega[j];
      const double k_fn1 = stiffnessK_fn1_[itype][jtype];
      const double ku_fn1 = stiffnessKu_fn1_[itype][jtype];
      const double kc_fn1 = stiffnessKc_fn1_[itype][jtype];
      const double k_fn2 = stiffnessK_fn2_[itype][jtype];
      const double ku_fn2 = stiffnessKu_fn2_[itype][jtype];
      const double kc_fn2 = stiffnessKc_fn2_[itype][jtype];
      const double k_ft = stiffnessK_ft_[itype][jtype];
      const double k_tn = stiffnessK_tn_[itype][jtype];
      const double ku_tn = stiffnessKu_tn_[itype][jtype];
      const double kc_tn = stiffnessKc_tn_[itype][jtype];
      const double k_tt = stiffnessK_tt_[itype][jtype];
      const double ku_tt = stiffnessKu_tt_[itype][jtype];
      const double kc_tt = stiffnessKc_tt_[itype][jtype];
      const double rb = lambda* (scdata.is_wall ? radi : std::min(radi,radj));
      const double A = M_PI*rb*rb;
      const double J = 0.5*A*rb*rb;
      const double I = 0.5*J;
      const double dt = update->dt;
      const double radsuminv = 1./(radi+radj);
      const double cri = scdata.is_wall ? (r) : (r * radi * radsuminv);           // corrected radius for i & j
      const double crj = r * radj * radsuminv;

      // calculate relative velocities at the contact point (for contact this is done in the surface_model)
      // relative translational velocity of i with respect to j
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
      {vectorScalarMult3D(omegai,0.5,wr);}
      else
      {
        // corrected for delta, wr = (radi*omegai + radj*omegaj) * rinv;
        vectorScalarMult3D(omegai,radi*radsuminv,tmp1);
        vectorScalarMult3D(omegaj,radj*radsuminv,tmp2);
        vectorAdd3D(tmp1,tmp2,wr);
      }
      // total relative velocities for shear, vtr = tangential translational velocity, vt + rotational velocity for shear, wr
      double vtr[3];
      vectorCross3D(delta,wr,tmp1);
      vectorAdd3D(vt,tmp1,vtr);

      // relative rotational velocity for torsion and bending
      vectorSubtract3D(omegai,omegaj,wr);
      // normal component for torsion
      double wn[3];
      vectorProject3D(wr,en,wn);
      // tangential component for bending
      double wt[3];
      vectorSubtract3D(wr,wn,wt);

      // -------------------- normal force calculation --------------------
      double nforce[3] = {0.,0.,0.}, dtforce[3] = {0.,0.,0.};
      double nforce_damped[3] = {0.,0.,0.}, tforce_damped[3] = {0.,0.,0.}, ntorque_damped[3] = {0.,0.,0.}, ttorque_damped[3] = {0.,0.,0.};
      const double displacement = (initialdist[0] - r);                           //negative: detached; positive: contacted

      // update history value
      double displacement_max, displacement_min;
      if (displacement > scdata.contact_history[history_offset_+14]){
        scdata.contact_history[history_offset_+14] = displacement;
        displacement_max = displacement;}
      else{
        displacement_max = scdata.contact_history[history_offset_+14];}
      if (displacement < scdata.contact_history[history_offset_+27]){
        scdata.contact_history[history_offset_+27] = displacement;
        displacement_min = displacement;}
      else{
        displacement_min = scdata.contact_history[history_offset_+27];}
      // set to zero for next loading-unloading cycle
      if (displacement < 0.0)
      {scdata.contact_history[history_offset_+14] = 0.0;}
      if (displacement > 0.0)
      {scdata.contact_history[history_offset_+27] = 0.0;}

      if (tensionflag_ || compressionflag_)
      {
        // tension or compression enabled
        if ( (tensionflag_ && displacement < -1.e-15) || (compressionflag_ && displacement > 1.e-15) )
        {
          //normal force - tracking distance explicitly
          double displacement_c1 = pow((ku_fn1-k_fn1)/(ku_fn1+kc_fn1),2.0)*displacement_max;
          double displacement_c2 = pow((ku_fn2-k_fn2)/(ku_fn2+kc_fn2),1.0)*displacement_min;
          double frcmag;

          if (displacement < displacement_min)
          {frcmag = k_fn2*A*displacement;}
          else
          {
              if (displacement < displacement_c2)
              {frcmag = ku_fn2*A*displacement+(k_fn2-ku_fn2)*A*displacement_min;}
              else
              {
                  if (displacement < 0.0)
                  {frcmag = -kc_fn2*A*displacement;}
                  else
                  {
                      if (displacement < displacement_c1)
                      {frcmag = -kc_fn1*pow(fabs(displacement),0.5);}
                      else
                      {
                          if (displacement < displacement_max)
                          {frcmag = ku_fn1*pow(fabs(displacement),0.5)+(k_fn1-ku_fn1)*pow(fabs(displacement_max),0.5);}
                          else
                          {frcmag = k_fn1*pow(fabs(displacement),0.5);}
                      }
                  }
              }
          }
          vectorScalarMult3D(en,frcmag,nforce);

          // damping
          // for positive vn (approaching particles) the damping force has to be repulsive (negative)
          if (dampingflag_)
          {
            const double minvel = 1e-5* (scdata.is_wall ? radi : fmin(radi,radj)) /dt;
            const double dvX = 0.01*nforce[0]*dt;
            const double dvY = 0.01*nforce[1]*dt;
            const double dvZ = 0.01*nforce[2]*dt;
            const double multiplierX = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vn[0]/fmax(minvel, dvX))) : MathExtraLiggghts::sgn(vn[0]);
            const double multiplierY = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vn[1]/fmax(minvel, dvY))) : MathExtraLiggghts::sgn(vn[1]);
            const double multiplierZ = dampingsmooth_ ? fmin(1.0, fmax(-1.0, vn[2]/fmax(minvel, dvZ))) : MathExtraLiggghts::sgn(vn[2]);
            nforce_damped[0] = nforce[0] - damping_fn_[itype][jtype]*fabs(nforce[0])*multiplierX;
            nforce_damped[1] = nforce[1] - damping_fn_[itype][jtype]*fabs(nforce[1])*multiplierY;
            nforce_damped[2] = nforce[2] - damping_fn_[itype][jtype]*fabs(nforce[2])*multiplierZ;
          } else {
            vectorCopy3D(nforce,nforce_damped);
          }
        }
      }

      // -------------------- shear force calculation --------------------
      if (shearflag_) {
        // calc change in shear forces
        vectorScalarMult3D(vtr,-k_ft*A*dt,dtforce);
        // rotate old forces (only non-explicity calculated one)
        // rotate tangential force
        vectorProject3D(force_tang,en,tmp1);
        vectorSubtract3D(force_tang,tmp1,force_tang);
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
            tforce_damped[0] = force_tang[0] - damping_ft_[itype][jtype]*fabs(force_tang[0])*multiplierX;
            tforce_damped[1] = force_tang[1] - damping_ft_[itype][jtype]*fabs(force_tang[1])*multiplierY;
            tforce_damped[2] = force_tang[2] - damping_ft_[itype][jtype]*fabs(force_tang[2])*multiplierZ;
        } else {
          vectorCopy3D(force_tang,tforce_damped);
        }
      }

      // -------------------- normal torque for twisting --------------------
      double theta_c1, theta_c2;
      double dtheta_normal[3] = {0.,0.,0.}, torque_normal[3] = {0.,0.,0.};
      if (ntorqueflag_)
      {
          vectorScalarMult3D(wn,dt,dtheta_normal);
          vectorAdd3D(theta_normal,dtheta_normal,theta_normal);
          vectorScalarMult3D(theta_normal,-k_tn*J,torque_normal);
          vectorProject3D(torque_normal,wn,torque_normal);

          if (theta_normal[0] > 0.0)
          {scdata.contact_history[history_offset_+21] = 0.0;}                     //theta_nx_min
          else{
            if (theta_normal[0] < 0.0)
            {scdata.contact_history[history_offset_+15] = 0.0;}                   //theta_nx_max
          }
          if (theta_normal[1] > 0.0)
          {scdata.contact_history[history_offset_+22] = 0.0;}                     //theta_ny_min
          else{
            if (theta_normal[1] < 0.0)
            {scdata.contact_history[history_offset_+16] = 0.0;}                   //theta_ny_max
          }
          if (theta_normal[2] > 0.0)
          {scdata.contact_history[history_offset_+23] = 0.0;}                     //theta_nz_min
          else{
            if (theta_normal[2] < 0.0)
            {scdata.contact_history[history_offset_+17] = 0.0;}                   //theta_nz_max
          }

          // normal X
          theta_c1 = (ku_tn-k_tn)/(ku_tn+kc_tn)*scdata.contact_history[history_offset_+15];
          theta_c2 = (ku_tn-k_tn)/(ku_tn+kc_tn)*scdata.contact_history[history_offset_+21];

          if ((theta_normal[0] >= scdata.contact_history[history_offset_+15]) || (theta_normal[0] <= scdata.contact_history[history_offset_+21]))
          {torque_normal[0] = -k_tn*J*theta_normal[0];}
          else{
              if (theta_normal[0] > theta_c1)
              {torque_normal[0] = -ku_tn*J*theta_normal[0]+(ku_tn-k_tn)*J*scdata.contact_history[history_offset_+15];}
               else{
                  if (theta_normal[0] >= theta_c2)
                  {torque_normal[0] = kc_tn*J*theta_normal[0];}
                  else
                  {torque_normal[0] = -ku_tn*J*theta_normal[0]+(ku_tn-k_tn)*J*scdata.contact_history[history_offset_+21];}
              }
          }
          // record theta_tx_max, theta_tx_min
          if (theta_normal[0] > scdata.contact_history[history_offset_+15])
            {scdata.contact_history[history_offset_+15] = theta_normal[0];}
          if (theta_normal[0] < scdata.contact_history[history_offset_+21])
            {scdata.contact_history[history_offset_+21] = theta_normal[0];}

          // normal Y
          theta_c1 = (ku_tn-k_tn)/(ku_tn+kc_tn)*scdata.contact_history[history_offset_+16];
          theta_c2 = (ku_tn-k_tn)/(ku_tn+kc_tn)*scdata.contact_history[history_offset_+22];

          if ((theta_normal[1] >= scdata.contact_history[history_offset_+16]) || (theta_normal[1] <= scdata.contact_history[history_offset_+22]))
          {torque_normal[1] = -k_tn*J*theta_normal[1];}
          else{
              if (theta_normal[1] > theta_c1)
              {torque_normal[1] = -ku_tn*J*theta_normal[1]+(ku_tn-k_tn)*J*scdata.contact_history[history_offset_+16];}
               else{
                  if (theta_normal[1] >= theta_c2)
                  {torque_normal[1] = kc_tn*J*theta_normal[1];}
                  else
                  {torque_normal[1] = -ku_tn*J*theta_normal[1]+(ku_tn-k_tn)*J*scdata.contact_history[history_offset_+22];}
              }
          }
          // record theta_tx_max, theta_tx_min
          if (theta_normal[1] > scdata.contact_history[history_offset_+16])
            {scdata.contact_history[history_offset_+16] = theta_normal[1];}
          if (theta_normal[1] < scdata.contact_history[history_offset_+22])
            {scdata.contact_history[history_offset_+22] = theta_normal[1];}

          // normal Z
          theta_c1 = (ku_tn-k_tn)/(ku_tn+kc_tn)*scdata.contact_history[history_offset_+17];
          theta_c2 = (ku_tn-k_tn)/(ku_tn+kc_tn)*scdata.contact_history[history_offset_+23];

          if ((theta_normal[2] >= scdata.contact_history[history_offset_+17]) || (theta_normal[2] <= scdata.contact_history[history_offset_+23]))
          {torque_normal[2] = -k_tn*J*theta_normal[2];}
          else{
              if (theta_normal[2] > theta_c1)
              {torque_normal[2] = -ku_tn*J*theta_normal[2]+(ku_tn-k_tn)*J*scdata.contact_history[history_offset_+17];}
               else{
                  if (theta_normal[2] >= theta_c2)
                  {torque_normal[2] = kc_tn*J*theta_normal[2];}
                  else
                  {torque_normal[2] = -ku_tn*J*theta_normal[2]+(ku_tn-k_tn)*J*scdata.contact_history[history_offset_+23];}
              }
          }
          // record theta_tx_max, theta_tx_min
          if (theta_normal[2] > scdata.contact_history[history_offset_+17])
            {scdata.contact_history[history_offset_+17] = theta_normal[2];}
          if (theta_normal[2] < scdata.contact_history[history_offset_+23])
            {scdata.contact_history[history_offset_+23] = theta_normal[2];}

          // damping
          if (dampingflag_)
          {
              ntorque_damped[0] = torque_normal[0] - damping_tn_[itype][jtype]*fabs(torque_normal[0])*MathExtraLiggghts::sgn(wn[0]);
              ntorque_damped[1] = torque_normal[1] - damping_tn_[itype][jtype]*fabs(torque_normal[1])*MathExtraLiggghts::sgn(wn[1]);
              ntorque_damped[2] = torque_normal[2] - damping_tn_[itype][jtype]*fabs(torque_normal[2])*MathExtraLiggghts::sgn(wn[2]);
          }
          else
          {vectorCopy3D(torque_normal,ntorque_damped);}
      }

      // -------------------- tangential torque for bending --------------------
      double dtheta_tang[3] = {0.,0.,0.}, torque_tang[3] = {0.,0.,0.};
      if (ttorqueflag_)
      {
          vectorScalarMult3D(wt,dt,dtheta_tang);
          vectorAdd3D(theta_tang,dtheta_tang,theta_tang);
          vectorScalarMult3D(theta_tang,-k_tt*I,torque_tang);
          vectorProject3D(torque_tang,wt,torque_tang);

          if (theta_tang[0] > 0.0)
          {scdata.contact_history[history_offset_+24] = 0.0;}                     //theta_tx_min
          else{
            if (theta_tang[0] < 0.0)
            {scdata.contact_history[history_offset_+18] = 0.0;}                   //theta_tx_max
          }
          if (theta_tang[1] > 0.0)
          {scdata.contact_history[history_offset_+25] = 0.0;}                     //theta_ty_min
          else{
            if (theta_tang[1] < 0.0)
            {scdata.contact_history[history_offset_+19] = 0.0;}                   //theta_ty_max
          }
          if (theta_tang[2] > 0.0)
          {scdata.contact_history[history_offset_+26] = 0.0;}                     //theta_tz_min
          else{
            if (theta_tang[2] < 0.0)
            {scdata.contact_history[history_offset_+20] = 0.0;}                   //theta_tz_max
          }

          // tangential torque X
          theta_c1 = (ku_tt-k_tt)/(ku_tt+kc_tt)*scdata.contact_history[history_offset_+18];
          theta_c2 = (ku_tt-k_tt)/(ku_tt+kc_tt)*scdata.contact_history[history_offset_+24];
          if ((theta_tang[0] >= scdata.contact_history[history_offset_+18]) || (theta_tang[0] <= scdata.contact_history[history_offset_+24]))
          {torque_tang[0] = -k_tt*I*theta_tang[0];}
          else{
              if (theta_tang[0] > theta_c1)
              {torque_tang[0] = -ku_tt*I*theta_tang[0]+(ku_tt-k_tt)*I*scdata.contact_history[history_offset_+18];}
               else{
                  if (theta_tang[0] >= theta_c2)
                  {torque_tang[0] = kc_tt*I*theta_tang[0];}
                  else
                  {torque_tang[0] = -ku_tt*I*theta_tang[0]+(ku_tt-k_tt)*I*scdata.contact_history[history_offset_+24];}
              }
          }
          // record theta_tx_max, theta_tx_min
          if (theta_tang[0] > scdata.contact_history[history_offset_+18])
            {scdata.contact_history[history_offset_+18] = theta_tang[0];}
          if (theta_tang[0] < scdata.contact_history[history_offset_+24])
            {scdata.contact_history[history_offset_+24] = theta_tang[0];}

          // tangential torque Y
          theta_c1 = (ku_tt-k_tt)/(ku_tt+kc_tt)*scdata.contact_history[history_offset_+19];
          theta_c2 = (ku_tt-k_tt)/(ku_tt+kc_tt)*scdata.contact_history[history_offset_+25];
          if ((theta_tang[1] >= scdata.contact_history[history_offset_+19]) || (theta_tang[1] <= scdata.contact_history[history_offset_+25]))
          {torque_tang[1] = -k_tt*I*theta_tang[1];}
          else{
              if (theta_tang[1] > theta_c1)
              {torque_tang[1] = -ku_tt*I*theta_tang[1]+(ku_tt-k_tt)*I*scdata.contact_history[history_offset_+19];}
               else{
                  if (theta_tang[1] >= theta_c2)
                  {torque_tang[1] = kc_tt*I*theta_tang[1];}
                  else
                  {torque_tang[1] = -ku_tt*I*theta_tang[1]+(ku_tt-k_tt)*I*scdata.contact_history[history_offset_+25];}
              }
          }
          // record theta_tx_max, theta_tx_min
          if (theta_tang[1] > scdata.contact_history[history_offset_+19])
            {scdata.contact_history[history_offset_+19] = theta_tang[1];}
          if (theta_tang[1] < scdata.contact_history[history_offset_+25])
            {scdata.contact_history[history_offset_+25] = theta_tang[1];}

          // tangential torque Z
          theta_c1 = (ku_tt-k_tt)/(ku_tt+kc_tt)*scdata.contact_history[history_offset_+20];
          theta_c2 = (ku_tt-k_tt)/(ku_tt+kc_tt)*scdata.contact_history[history_offset_+26];
          if ((theta_tang[2] >= scdata.contact_history[history_offset_+20]) || (theta_tang[2] <= scdata.contact_history[history_offset_+26]))
          {torque_tang[2] = -k_tt*I*theta_tang[2];}
          else{
              if (theta_tang[2] > theta_c1)
              {torque_tang[2] = -ku_tt*I*theta_tang[2]+(ku_tt-k_tt)*I*scdata.contact_history[history_offset_+20];}
               else{
                  if (theta_tang[2] >= theta_c2)
                  {torque_tang[2] = kc_tt*I*theta_tang[2];}
                  else
                  {torque_tang[2] = -ku_tt*I*theta_tang[2]+(ku_tt-k_tt)*I*scdata.contact_history[history_offset_+26];}
              }
          }
          // record theta_tx_max, theta_tx_min
          if (theta_tang[2] > scdata.contact_history[history_offset_+20])
            {scdata.contact_history[history_offset_+20] = theta_tang[2];}
          if (theta_tang[2] < scdata.contact_history[history_offset_+26])
            {scdata.contact_history[history_offset_+26] = theta_tang[2];}

          // damping
          if (dampingflag_)
          {
              ttorque_damped[0] = torque_tang[0] - damping_tt_[itype][jtype]*fabs(torque_tang[0])*MathExtraLiggghts::sgn(wt[0]);
              ttorque_damped[1] = torque_tang[1] - damping_tt_[itype][jtype]*fabs(torque_tang[1])*MathExtraLiggghts::sgn(wt[1]);
              ttorque_damped[2] = torque_tang[2] - damping_tt_[itype][jtype]*fabs(torque_tang[2])*MathExtraLiggghts::sgn(wt[2]);
          }
          else
          {vectorCopy3D(torque_tang,ttorque_damped);}
      }
      // -----------------------------------------------------------------------

      // bond breakage due to stress (use un-damped forces/torques)
      if (breakmode_ == BREAK_STRESS)
      {
        const double nforce_mag = vectorMag3D(nforce);                            //TODO: Potyondy and Cundall do not use the absolut value
        const double tforce_mag = vectorMag3D(force_tang);
        const double ntorque_mag = vectorMag3D(torque_normal);
        const double ttorque_mag = vectorMag3D(torque_tang);
        const double maxSigma = maxSigma_[itype][jtype] * ((useRatioTC_ && displacement < -1.e-15) ? ratioTensionCompressionBondnonlinear_[itype][jtype] : 1.0);
        const bool nstress = maxSigma < (nforce_mag/A + ttorque_mag*rb/I);
        const double maxTau = maxTau_[itype][jtype];
        const bool tstress = maxTau < (tforce_mag/A + ntorque_mag*rb/J);
        if ((nstress || tstress) && update_history)
        {
          breakBondnonlinear(scdata);
          return;                                                                 // bond broken
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
        vectorCopy3D(theta_normal, theta_normal_ptr);
        vectorCopy3D(theta_tang, theta_tang_ptr);
      }

      // return resulting forces (contact model overriden)
      if(scdata.is_wall) {
        i_forces.delta_F[0] = F[0];
        i_forces.delta_F[1] = F[1];
        i_forces.delta_F[2] = F[2];
        i_forces.delta_torque[0] = dTorque_i[0];
        i_forces.delta_torque[1] = dTorque_i[1];
        i_forces.delta_torque[2] = dTorque_i[2];
      } else {
        i_forces.delta_F[0] = F[0];
        i_forces.delta_F[1] = F[1];
        i_forces.delta_F[2] = F[2];
        i_forces.delta_torque[0] = dTorque_i[0];
        i_forces.delta_torque[1] = dTorque_i[1];
        i_forces.delta_torque[2] = dTorque_i[2];

        j_forces.delta_F[0] = -F[0];
        j_forces.delta_F[1] = -F[1];
        j_forces.delta_F[2] = -F[2];
        j_forces.delta_torque[0] = dTorque_j[0];
        j_forces.delta_torque[1] = dTorque_j[1];
        j_forces.delta_torque[2] = dTorque_j[2];
      }
    }

    /* ---------------------------------------------------------------------- */
    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      // check if cohesion model is working for this type-pair
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double lambda = radiusMultiplierBondnonlinear_[itype][jtype];
      if (lambda < 1.e-15)
        return;                                                                   // do nothing in case of rb=0.0
      surfacesClose(sidata, i_forces, j_forces);
    }

    void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}

  /* ====================================================================== */
  private:
    void createBondnonlinear(SurfacesCloseData & scdata, double dist)
    {
      double * const bondFlag = &scdata.contact_history[history_offset_];
      double * const initialdist = &scdata.contact_history[history_offset_+1];
      double * const contact_pos = &scdata.contact_history[history_offset_+2];    // 3-elements
      double * const ft_old = &scdata.contact_history[history_offset_+5];         // 3-elements
      double * const theta_normal = &scdata.contact_history[history_offset_+8];   // 3-elements
      double * const theta_tang = &scdata.contact_history[history_offset_+11];    // 3-elements
      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
      bondFlag[0] = 1.0;
      initialdist[0] = dist;                                                      // use current distance
      vectorSubtract3D(atom->x[scdata.i], scdata.delta, contact_pos);
      vectorZeroize3D(ft_old);
      vectorZeroize3D(theta_normal);
      vectorZeroize3D(theta_tang);
      if (compute_bond_counter_)
          compute_bond_counter_->bond_created(scdata.is_wall, scdata.i, scdata.j);
    }

    void breakBondnonlinear(SurfacesCloseData & scdata)
    {
      double * const bondFlag = &scdata.contact_history[history_offset_];
      double * const initialdist = &scdata.contact_history[history_offset_+1];

      if (scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
      bondFlag[0] = 0.;
      initialdist[0] = 0.;

      if (compute_bond_counter_)
          compute_bond_counter_->bond_broken(scdata.is_wall, scdata.i, scdata.j);
    }

    int history_offset_;
    double ** radiusMultiplierBondnonlinear_;
    double ** stiffnessK_fn1_;
    double ** stiffnessKu_fn1_;
    double ** stiffnessKc_fn1_;
    double ** stiffnessK_fn2_;
    double ** stiffnessKu_fn2_;
    double ** stiffnessKc_fn2_;
    double ** stiffnessK_ft_;
    double ** stiffnessK_tn_;
    double ** stiffnessKu_tn_;
    double ** stiffnessKc_tn_;
    double ** stiffnessK_tt_;
    double ** stiffnessKu_tt_;
    double ** stiffnessKc_tt_;
    double ** damping_fn_;
    double ** damping_ft_;
    double ** damping_tn_;
    double ** damping_tt_;
    double ** maxDist_;
    double ** maxSigma_;
    double ** maxTau_;
    double ** ratioTensionCompressionBondnonlinear_;
    double tsCreateBondnonlinear_;
    double ** createDistanceBondnonlinear_;

  /* ====================================================================== */
  protected:
    bool stressbreakflag_; // false if max dist break, true if stress break
    bool tensionflag_;
    bool compressionflag_;
    bool shearflag_;
    bool ntorqueflag_;
    bool ttorqueflag_;
    bool createBondNonlinear_always_flag_;
    bool dampingflag_;
    bool dampingsmooth_;
    int breakmode_;
    bool useRatioTC_;
    class FixPropertyAtom *fix_bondnonlinear_random_id_;
    class ComputeBondCounter *compute_bond_counter_; // compute that keeps track of the number of bonds and created/destroyed ones
  };
}
}
#endif // COHESION_MODEL_BOND_NONLINEAR_H_
#endif
