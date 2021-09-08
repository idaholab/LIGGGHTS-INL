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

------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(HYSTERETIC_NONLINEAR2,hysteretic/nonlinear2,18)
#else
#ifndef NORMAL_MODEL_HYSTERETIC_NONLINEAR2_H_
#define NORMAL_MODEL_HYSTERETIC_NONLINEAR2_H_
#include "contact_models.h"
#include "normal_model_base.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"
#include <algorithm>
#include "math_extra_liggghts.h"
#include <math.h>

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HYSTERETIC_NONLINEAR2> : public NormalModelBase
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
      Alpha(NULL),
      Cin(NULL),
      A1(NULL),
      A2(NULL),
      A3(NULL),
      kcin(NULL),
      limitForce(false)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0");
      hsetup->add_history_value("deltaZero","0"); //history[1]
      hsetup->add_history_value("k1", "0"); //history[2]
      hsetup->add_history_value("deltaZero_old","0"); //history[3]
      hsetup->add_history_value("k1_old","0"); //history[4]
      hsetup->add_history_value("delta_old","0"); //history[5]
      hsetup->add_history_value("deltaMin","0"); //history[6]
      hsetup->add_history_value("betan","0"); //history[7]
      hsetup->add_history_value("f0","0"); //history[8]
      hsetup->add_history_value("f0_old","0"); //history[9]
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
      registry.registerProperty("K_elastic", &MODEL_PARAMS::createLoadingStiffness,"model hysteretic/nonlinear2");
      registry.registerProperty("CoeffRestLog", &MODEL_PARAMS::createCoeffRestLog,"model hysteretic/nonlinear2");
      registry.registerProperty("kn2k1", &MODEL_PARAMS::createUnloadingStiffness,"model hysteretic/nonlinear2");
      registry.registerProperty("kn2kc", &MODEL_PARAMS::createCoeffAdhesionStiffness,"model hysteretic/nonlinear2");
      registry.registerProperty("phiF", &MODEL_PARAMS::createCoeffPlasticityDepth,"model hysteretic/nonlinear2");
      registry.registerProperty("f_adh", &MODEL_PARAMS::createPullOffForce,"model hysteretic/nonlinear2");
      //hysteretic nonlinear2 model properties
      registry.registerProperty("Alpha", &MODEL_PARAMS::createAlphaCustom,"model hysteretic/nonlinear2");
      registry.registerProperty("Cin", &MODEL_PARAMS::createCinCustom,"model hysteretic/nonlinear2");
      registry.registerProperty("A1", &MODEL_PARAMS::createAoneCustom,"model hysteretic/nonlinear2");
      registry.registerProperty("A2", &MODEL_PARAMS::createAtwoCustom,"model hysteretic/nonlinear2");
      registry.registerProperty("A3", &MODEL_PARAMS::createAthreeCustom,"model hysteretic/nonlinear2");
      registry.registerProperty("kcin", &MODEL_PARAMS::createKcinCustom,"model hysteretic/nonlinear2");

      //

      registry.connect("K_elastic", K_elastic,"model hysteretic/nonlinear2");
      registry.connect("CoeffRestLog", CoeffRestLog,"model hysteretic/nonlinear2");
      registry.connect("kn2kc", kn2kc,"model hysteretic/nonlinear2");
      registry.connect("kn2k1", kn2k1, "model hysteretic/nonlinear2");
      registry.connect("phiF", phiF,"model hysteretic/nonlinear2");
      registry.connect("f_adh", f_adh,"model hysteretic/nonlinear2");
      //hysteretic nonlinear2 model properties
      registry.connect("Alpha", Alpha,"model hysteretic/nonlinear2");
      registry.connect("Cin", Cin,"model hysteretic/nonlinear2");
      registry.connect("A1", A1,"model hysteretic/nonlinear2");
      registry.connect("A2", A2,"model hysteretic/nonlinear2");
      registry.connect("A3", A3,"model hysteretic/nonlinear2");
      registry.connect("kcin", kcin,"model hysteretic/nonlinear2");

      //

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hysteretic/nonlinear2");
    }

    // effective exponent for stress-strain relationship

    double stressStrainExponent()
    {
      return 1.;
    }

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int id_i = sidata.i;
      const int id_j = sidata.j;
      double * vi = sidata.v_i;
      double * vj = sidata.v_j;
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double deltan = sidata.deltan;
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : (ri*rj/(ri+rj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
#endif
      double meff=sidata.meff;
      double kn = K_elastic[itype][jtype];
      // double kt = kn;

      // input parameters for hysteretic/nonlinear2 model: Alpha, Cin, A1, A2, A3, CoeffGn, k_c
      //double Alpha_in = 20;
      //double Cin = 1.0e-7;
      //double A1 = 7.0e8; // 1.0e9; // 5.0e7; - stiffness increasing rate
      //double A2 = 4e4; // 3e4; //1.0e6; - initial loading stiffness
      //double A3 = 5.0;

      double Alpha_in = Alpha[itype][jtype];
      double Cin_in = Cin[itype][jtype];
      double A1_in = A1[itype][jtype];
      double A2_in = A2[itype][jtype];
      double A3_in = A3[itype][jtype];
      double kcin_in = kcin[itype][jtype];

      double CoeffGn = 10.5*exp(0.6*CoeffRestLog[itype][jtype]);
      double k_c = kcin_in * A2_in;
      double fcc = 1.0;

      // convert Kn and Kt from pressure units to force/distance^2
      // kn /= force->nktv2p;
      // kt /= force->nktv2p;

      //const double k1 = kn;
      const double k2Max = kn * kn2k1[itype][jtype];

      const double kc = kn2kc[itype][jtype] * kn;
      const double f_0 = f_adh[itype][jtype];

      double gamman, gammat;

      //gamman = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));
      //gammat = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));

      //if (!tangential_damping) gammat = 0.0;

      // get the history value -- 0.deltaMax, 1.deltaZero, 2.k1, 3.deltaZero_old, 4.k1_old, 5.deltan_old, 6.deltaMin, 7.betan, 8.f0, 9.f0_old
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      double * const kc_history = &sidata.contact_history[kc_offset];
      double * const fo_history = &sidata.contact_history[fo_offset];
      double deltaMax; // the 4th value of the history array is deltaMax
      // initial values
      if (deltan > history[0]) {
          history[0] = deltan;
          deltaMax = deltan;
      } else{
        deltaMax = history[0];
      }
      double deltaZero = history[1];
      double k1 = history[2];
      //double k1 = A2;
      //if (history[2] == 0)
      //{
      //   k1 = A2;
      //   history[2] = k1;
      //}
      //else
      //{
      //   k1 = history[2];
      //}
      double deltaZero_old = history[3];
      double k1_old = history[4];
      //double k1_old = A2;
      //if (k1 >= history[4])
      //{
      //   k1_old = history[4];
      //}
      //else
      //{
      //   k1_old = k1;
      //   history[4] = k1_old;
      //}
      double delta_old = history[5];
      double deltaMin = history[6];
      double betan = history[7];
      double f0 = history[8];
      double f0_old = history[9];


      // k2 dependent on the maximum overlap
      // this accounts for an increasing stiffness with deformation - to capture nonlinearity
      const double deltaMaxLim =(k2Max/(k2Max-k1))*phiF[itype][jtype]*2*reff;

      double k2, fHys;
      k2 = A3_in*k1;

// - contact status, loading - 1, unloading - 2
      int tag_status;
      if (deltan >= delta_old)
      {
         if (delta_old == 0)
         {
            tag_status = 1;
         }
         else
         {
         if (sidata.vn > 0)
            {
               tag_status = 2;
            }
         else
            {
               tag_status = 1;
            }
         }
      }
      else
      {
         tag_status =2;
      }

// - contact force calculation beginning
      if (tag_status == 1 ) //(deltan>=delta_old) //loading part
      {
         if (deltaZero == 0) //initial loading
         {
            k1 = A2_in;
         }
         else
         {
            k1 = history[2];
         }
         deltaMin = history[6];
         f0_old = f0;
         if (deltan <= deltaZero)
         {
            //if (deltan <= deltaMin)
            //   {
            //       fHys = -k_c*deltan;
            //   }
            //else
            //   {
            //       fHys = betan*(deltan-deltaZero);
            //   }
	    fHys = Cin_in*k2*(exp(betan*(deltan-deltaZero))-1);
         }
         else
         {
            fHys = Alpha_in*k1*pow(deltan-deltaZero,2)+f0;
         }
         //calculate damping coefficient
         gamman = sqrt(5/4)*sqrt(4.*meff*Alpha_in*k1/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])))*(pow(deltan,0.5)+pow(deltaZero,0.5));
         gammat = sqrt(5/4)*sqrt(4.*meff*Alpha_in*k1/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])))*(pow(deltan,0.5)+pow(deltaZero,0.5));
         //update historic values
         history[1] = deltaZero;
         history[2] = k1;
         history[3] = deltaZero;
         history[4] = k1;
         history[5] = deltan;
         history[6] = deltaMin;
         history[7] = betan;
         history[8] = f0;
         history[9] = f0_old;
      }
      else //deltan < delta_old, unloading part
      {
         k2 = A3_in*(A1_in*deltaMax+A2_in); //unloading stiffness
         deltaZero = fcc*(1-k1_old/k2)*deltaMax; //plastic deformation
	 deltaZero_old = history[3];
         betan = log(Alpha_in*k1_old/Cin_in/k2*pow(deltaMax-deltaZero_old,2)+1)/(deltaMax-fcc*(1-k1_old/k2)*deltaMax); //calculor
         deltaMin = betan*(k2-k1_old)/(betan*k2+k_c)*deltaMax;
         k1 = deltaMax*A1_in+A2_in; //updated loading stiffness

	 //if (deltan >= deltaMin)
         //   {
               if (deltan >= deltaZero)
               {
                   fHys = Cin_in*k2*(exp(betan*(deltan-deltaZero))-1) +(deltan-deltaZero)*f0_old/(deltaMax-deltaZero);
                   f0 = fHys - Alpha_in*k1*pow(deltan-deltaZero,2);
		   //std::cout << std::scientific;
		   //std::cout << "delta " << deltan << " fHys " << fHys << std::endl;
               }
               else
               {
                   fHys = Cin_in*k2*(exp(betan*(deltan-deltaZero))-1); 
                   f0 = 0;
               }
         //   }
         //else
         //   {
         //      fHys = -k_c*deltan;
         //      f0 = 0;
         //   }

         //calculate damping coefficient
         gamman = 1*0.001*sqrt(5/4)*sqrt(4.*meff*Alpha_in*k1_old/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])))*(pow(deltan,-0.25));
         gammat = 1*0.001*sqrt(5/4)*sqrt(4.*meff*Alpha_in*k1_old/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])))*(pow(deltan,-0.25));
         //update historic values
         history[1] = deltaZero;
         history[2] = k1;
         history[3] = deltaZero_old;
         history[4] = k1_old;
         history[5] = deltan;
         history[6] = deltaMin;
         history[7] = betan;
         history[8] = f0;
      }

// - stiffness modulues calculate
      kn = k1;
      double kt = kn;
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

// - contact force calculation end

      const double Fn_damping = -gamman*sidata.vn;
      double Fn = fHys + Fn_damping + f_0;

      //std::cout << std::scientific;
      //std::cout << "id_i " << id_i << " id_j " << id_j << " d_n " << deltan << " Fhys " << fHys << " f_damping " << Fn_damping << " gamman " << gamman << " vj_3 " << vj[2] << std::endl;


      if(limitForce && (Fn<0.0) && kc == 0 && f_0 == 0.0){
          Fn = 0.0;
      }
      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.deltaZero = deltaZero;
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
      history[6] = 0.0;
      history[7] = 0.0;
      history[8] = 0.0;
      history[9] = 0.0;
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
    double **A1;
    double **A2;
    double **A3;
    double **Alpha;
    double **Cin;
    double **kcin;

    int history_offset;
    int kc_offset;
    int fo_offset;

    bool tangential_damping;
    bool limitForce;
  };
}
}
#endif // NORMAL_MODEL_HYSTERETIC_NONLINEAR2_H_
#endif
