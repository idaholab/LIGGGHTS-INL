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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Christoph Kloss (DCS Computing GmbH, JKU Linz)
    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef MESHMODULE_CLASS

MeshModuleStyle(6dof,MeshModuleStress6DOF)
// do not allow the combination of the following
MeshModuleRestrict(6dof,servo)

#else

#ifndef LMP_MESH_MODULE_STRESS_6DOF_H
#define LMP_MESH_MODULE_STRESS_6DOF_H

#include "fix.h"
#include "input.h"
#include <cmath>
#include "mesh_module_stress.h"
#include "mesh_module.h"
#include "fix_gravity.h"

namespace LAMMPS_NS
{

class MeshModuleStress6DOF : public MeshModule
{

    public:

      MeshModuleStress6DOF(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh);
      virtual ~MeshModuleStress6DOF();

      virtual void post_create();

      void init();
      virtual void setup(int vflag);
      virtual void setup_pre_force(int vflag);

      int setmask();
      void initial_integrate(int vflag);
      void final_integrate();

      double compute_vector(int n);

      inline int get_num_vector_components() const
      { return 10; }

      void get_xcm(double *xcm)
      { vectorCopy3D(xcm_(0),xcm); }
      void get_dX_dx(double *dX, double *dx)
      { vectorCopy3D(dX_,dX); vectorCopy3D(dx_,dx); }
      void get_dQ_dq(double *dQ, double *dq)
      { vectorCopy4D(dQ_,dQ); vectorCopy4D(dq_,dq); }

    private:

      void init_defaults();
      void error_checks();

      void init_rotation_props();
      void calc_displace();

      void set_vel();
      void add_gravity();
      void add_suspension_force();
      void add_external_force();
      void rot_flip(); 

      int modify_param(int, char **);

      // properties of 6 dof body

      VectorContainer<double,3> &xcm_;
      VectorContainer<double,4> &quat_;
      VectorContainer<double,3> &vcm_;
      VectorContainer<double,3> &omega_;
      VectorContainer<double,3> &angmom_;
      ScalarContainer<double> &mass_;
      MultiVectorContainer<double,3,3> &moi_;
      VectorContainer<double,3> &inertia_;
      VectorContainer<double,3> &ex_space_;
      VectorContainer<double,3> &ey_space_;
      VectorContainer<double,3> &ez_space_;

      // original position and orientation
      
      VectorContainer<double,3> &xcm_orig_;
      VectorContainer<double,4> &quat_orig_;
      
      VectorContainer<double,3> &xcm_orig_0_;
      VectorContainer<double,4> &quat_orig_0_;

      // timesteps and flags for integration

      double dtf_,dtv_,dtfm_,dtq_;
      VectorContainer<bool,3> &fflag_;
      VectorContainer<bool,3> &tflag_;

      // displace for each mesh node

      MultiVectorContainer<double,3,3> &displace_;

      // velocity for each node

      MultiVectorContainer<double,3,3> &v_;

      // additional constraint
      bool suspension_flag_;
      double k_t_, c_t_, k_r_, c_r_;

      bool externalForce_flag_;
      double xvalue, yvalue, zvalue;
      char *xstr, *ystr, *zstr;
      int xvar, yvar, zvar, xstyle, ystyle, zstyle;

      bool rot_flip_flag_;
      double rot_flip_angle_;

      FixGravity *fix_gravity_;
      bool gravity_set_;
      double limit_vel_;
      double axis_[3];
      bool limit_movement_axis_;

      MeshModuleStress *mm_stress;

      // made these class variables because their values are
      // used by MeshMoverFollow6DOF
      double dX_[3], dx_[3], dQ_[4], dq_[4], qOld_[4];

}; //end class

}

#endif
#endif
