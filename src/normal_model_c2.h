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

    Rahul Mohanty (University of Edinburgh, P&G)
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(C2,c2,16)
#else
#ifndef NORMAL_MODEL_C2_H_
#define NORMAL_MODEL_C2_H_
#include "contact_models.h"
#include "normal_model_base.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"

// namespace MODEL_PARAMS
// {
//     inline static ScalarProperty* createCoeffAlphaScalarC2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
//     {
//       ScalarProperty* CoeffAlphaScalar = MODEL_PARAMS::createScalarProperty(registry, "Alpha", caller);
//       return CoeffAlphaScalar;
//     }
//     inline static ScalarProperty* createCoeffCinScalarC2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
//     {
//       ScalarProperty* CoeffCinScalar = MODEL_PARAMS::createScalarProperty(registry, "Cin", caller);
//       return CoeffCinScalar;
//     }
//     inline static ScalarProperty* createCoeffA1ScalarC2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
//     {
//       ScalarProperty* CoeffA1Scalar = MODEL_PARAMS::createScalarProperty(registry, "A1", caller);
//       return CoeffA1Scalar;
//     }
//     inline static ScalarProperty* createCoeffA2ScalarC2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
//     {
//       ScalarProperty* CoeffA2Scalar = MODEL_PARAMS::createScalarProperty(registry, "A2", caller);
//       return CoeffA2Scalar;
//     }
//     inline static ScalarProperty* createCoeffA3ScalarC2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
//     {
//       ScalarProperty* CoeffA3Scalar = MODEL_PARAMS::createScalarProperty(registry, "A3", caller);
//       return CoeffA3Scalar;
//     }
// }

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<C2> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      K_elastic(NULL),
      CoeffRestLog(NULL),
      kn2k1(NULL),
      kn2kc(NULL),
      phiF(NULL),
      f_adh(NULL),
//      Alpha(NULL),
//      Cin(NULL),
//      A1(NULL),
//      A2(NULL),
//      A3(NULL),
      limitForce(false)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0"); //history[0]
      hsetup->add_history_value("deltaZero", "0"); //history[1]
      hsetup->add_history_value("k1", "0"); //history[2]
      hsetup->add_history_value("deltaZero_old","0"); //history[3]
      hsetup->add_history_value("k1_old","0"); //history[4]
      hsetup->add_history_value("delta_old","0"); //history[5]
      kc_offset = hsetup->add_history_value("kc", "1");
      fo_offset = hsetup->add_history_value("fo", "1");
      c->add_history_offset("kc_offset", kc_offset);
      c->add_history_offset("fo_offset", fo_offset);
    }

    void registerSettings(Settings & settings){
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce, true);
    }
    void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("K_elastic", &MODEL_PARAMS::createLoadingStiffness,"model c2"); //initial k1 in 1st cyclicloading
      registry.registerProperty("CoeffRestLog", &MODEL_PARAMS::createCoeffRestLog,"model c2");
      registry.registerProperty("kn2k1", &MODEL_PARAMS::createUnloadingStiffness,"model c2"); //k2 = k1 * kn2k1
      registry.registerProperty("kn2kc", &MODEL_PARAMS::createCoeffAdhesionStiffness,"model c2"); //kc = k1 * kn2kc
      registry.registerProperty("phiF", &MODEL_PARAMS::createCoeffPlasticityDepth,"model c2");
      registry.registerProperty("f_adh", &MODEL_PARAMS::createPullOffForce,"model c2");
//      registry.registerProperty("Alpha", &MODEL_PARAMS::createCoeffAlphaC2,"model c2"); //model parameter a
//      registry.registerProperty("Cin", &MODEL_PARAMS::createCoeffCinC2,"model c2"); //model parameter C
//      registry.registerProperty("A1", &MODEL_PARAMS::createCoeffA1C2,"model c2"); //model parameter A1
//      registry.registerProperty("A2", &MODEL_PARAMS::createCoeffA2C2,"model c2"); //model parameter A2
//      registry.registerProperty("A3", &MODEL_PARAMS::createCoeffA3C2,"model c2"); //model parameter A3

      registry.connect("K_elastic", K_elastic,"model c2");
      registry.connect("CoeffRestLog", CoeffRestLog,"model c2");
      registry.connect("kn2kc", kn2kc,"model c2");
      registry.connect("kn2k1", kn2k1, "model c2");
      registry.connect("phiF", phiF,"model c2");
      registry.connect("f_adh", f_adh,"model c2");
//      registry.connect("Alpha", Alpha, "model c2");
//      registry.connect("Cin", Cin, "model c2");
//      registry.connect("A1", A1, "model c2");
//      registry.connect("A2", A2, "model c2");
//      registry.connect("A3", A3, "model c2");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model c2");
    }

    // effective exponent for stress-strain relationship

    double stressStrainExponent()
    {
      return 1.;
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double deltan = sidata.deltan; //read overlap distance
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : (ri*rj/(ri+rj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif
      double meff=sidata.meff;
      double kn = K_elastic[itype][jtype];
      double kt = kn;

      // input paramters: Alpha, Cin, A1, A2, A3
      double Alpha = 20.0;
      double Cin = 1.0e-7;
      double A1 = 5.0e7;
      double A2 = 1.0e6;
      double A3 = 2.5;
      double CoeffGn = 10.5*exp(0.6*CoeffRestLog[itype][jtype]);

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      // k2Max, kc f_0 will not be used in this implemented model
      const double k2Max = kn * kn2k1[itype][jtype];
      const double kc = kn2kc[itype][jtype] * kn;
      const double f_0 = f_adh[itype][jtype];

      // calculate damping coefficient
      double gamman, gammat;
      // gamman = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));
      // gammat = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));
      // gammat = 0.0;

      if (!tangential_damping) gammat = 0.0;

      // get the history value -- 0. deltaMax, 1. deltaZero, 2. k_1, 3. deltaZero_old, 4. k1_old, 5. deltan_old
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      double * const kc_history = &sidata.contact_history[kc_offset];
      double * const fo_history = &sidata.contact_history[fo_offset];
      // the 4th value of the history array is deltaMax
      double deltaMax;
      if (deltan >= history[0])
         {
         history[0] = deltan;
         deltaMax = deltan;
         }
      else
         {
         deltaMax = history[0];
         }
      double deltaZero = history[1];
      double k1 = history[2];
      double deltaZero_old = history[3];
      double k1_old = history[4];
      double delta_old = history[5]; //overlap for previous timestep

      double k2 = A3*k1;
      double fHys;

       int tag_status;
       if (deltan >= delta_old)
       {
          if (sidata.vn > 0)
          {
             tag_status = 2; //unloading
          }
          else
          {
             tag_status = 1; //loading
          }
       }
       else
       {
          tag_status = 2; //unloading
       }

       if (tag_status == 1)//(deltan >= delta_old) //loading part
       {
         if (deltaZero == 0) //initial loading
         {
            k1 = A2;
            history[2] = A2;
         }
         else
         {
            k1 = history[2];
         }
         if (deltan <= deltaZero)
         {
            fHys = 0.0;
         }
         else
         {
            fHys = Alpha*k1*pow(deltan-deltaZero,2);
         }
         history[1] = deltaZero;
         history[3] = deltaZero;
         history[4] = k1;
         history[5] = deltan;
         gamman = sqrt(5/4)*sqrt(4.*meff/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])))*sqrt(Alpha*k1)*pow(deltan,0.25);
         //gamman = CoeffGn*sqrt(4.*meff/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])))*sqrt(Alpha*k1)*pow(deltan,0.5);
         gammat = gamman;
       }
       else // deltan < deltaMax, unloading part
       {
         k2 = A3*(A1*deltaMax+A2);
         deltaZero_old = history[3];
         k1_old = history[4];
         double beta = log(Alpha*k1_old/Cin/k2*pow(deltaMax-deltaZero_old,2)+1)/(deltaMax-(1-k1_old/k2)*deltaMax);
         deltaZero = (1-k1_old/k2)*deltaMax;
         if (deltan <= deltaZero)
         {
            fHys = 0.0;
         }
         else
         {
            fHys = Cin*k2*(exp(beta*(deltan-deltaZero))-1.0);
         }
         k1 = A1*deltaMax+A2;
         history[1] = deltaZero;
         history[2] = k1;
         history[5] = deltan;
         gamman = sqrt(5/4)*sqrt(4.*meff/(1.+(M_PI/CoeffRestLog[itype][jtype])/(M_PI/CoeffRestLog[itype][jtype])))*sqrt(Alpha*k1_old)*pow(deltan,0.25);
         //gamman = 5*CoeffGn*sqrt(4.*meff/(1.+(M_PI/CoeffRestLog[itype][jtype])/(M_PI/CoeffRestLog[itype][jtype])))*sqrt(k2)*pow(deltan,0.5);
         gammat = gamman;
       }

      const double Fn_damping = -gamman*sidata.vn;
      double Fn = fHys + Fn_damping + f_0;

      if(limitForce && (Fn<0.0) && kc == 0 && f_0 == 0.0){
          Fn = 0.0;
      }
      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      kc_history[0] = kc;
      fo_history[0] = f_0;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

#ifdef NONSPHERICAL_ACTIVE_FLAG
      double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
      double torque_i[3] = {0.0, 0.0, 0.0}; //initialized here with zeros to avoid compiler warnings
      if(sidata.is_non_spherical) {
        double xci[3];
        vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
        vectorCross3D(xci, Fn_i, torque_i);
      }
#endif

      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] = Fn_ * sidata.en[0];
        i_forces.delta_F[1] = Fn_ * sidata.en[1];
        i_forces.delta_F[2] = Fn_ * sidata.en[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
        if(sidata.is_non_spherical) {
          //for non-spherical particles normal force can produce torque!
          i_forces.delta_torque[0] += sidata.area_ratio*torque_i[0];
          i_forces.delta_torque[1] += sidata.area_ratio*torque_i[1];
          i_forces.delta_torque[2] += sidata.area_ratio*torque_i[2];
        }
        #endif
      } else {
        i_forces.delta_F[0] = sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] = sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] = sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
#ifdef NONSPHERICAL_ACTIVE_FLAG
        if(sidata.is_non_spherical) {
          //for non-spherical particles normal force can produce torque!
          double xcj[3], torque_j[3];
          double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
          vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
          vectorCross3D(xcj, Fn_j, torque_j);

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

    void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_NORMAL_MODEL;
      double * const history = &scdata.contact_history[history_offset];
      history[0] = 0.0;
      history[1] = 0.0;
      history[2] = 0.0;
      history[3] = 0.0;
      history[4] = 0.0;
      history[5] = 0.0;
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double **K_elastic;
    double **CoeffRestLog;
    double **kn2k1;
    double **kn2kc;
    double **phiF;
    double **f_adh;

    int history_offset;
    int kc_offset;
    int fo_offset;

    bool tangential_damping;
    bool limitForce;
  };
}
}
#endif // NORMAL_MODEL_TEST_H_
#endif
