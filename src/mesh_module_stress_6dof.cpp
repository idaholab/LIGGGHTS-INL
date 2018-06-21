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

#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include "mesh_module_stress_6dof.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "fix_rigid.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include <mpi.h>
#include "vector_liggghts.h"
#include "fix_property_global.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define EPSILON 1.0e-7

/* ---------------------------------------------------------------------- */

MeshModuleStress6DOF::MeshModuleStress6DOF(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh) :
    MeshModule(lmp, iarg_, narg, arg, fix_mesh),

  xcm_(      *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm","comm_none","frame_invariant","restart_yes",3)),
  quat_(     *mesh->prop().addGlobalProperty< VectorContainer<double,4> > ("quat","comm_none","frame_invariant","restart_yes",1)),
  vcm_(      *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("vcm","comm_none","frame_invariant","restart_yes",1)),
  omega_(    *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("omega","comm_none","frame_invariant","restart_yes",1)),
  angmom_(   *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("angmom","comm_none","frame_invariant","restart_yes",1)),
  mass_(     *mesh->prop().addGlobalProperty< ScalarContainer<double> >   ("mass","comm_none","frame_invariant","restart_yes",3)),
  moi_(      *mesh->prop().addGlobalProperty< MultiVectorContainer<double,3,3> > ("moi","comm_none","frame_invariant","restart_yes",2)),
  inertia_(  *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("inertia","comm_none","frame_invariant","restart_yes",2)),
  ex_space_( *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("ex_space","comm_none","frame_invariant","restart_yes",2)),
  ey_space_( *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("ey_space","comm_none","frame_invariant","restart_yes",2)),
  ez_space_( *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("ez_space","comm_none","frame_invariant","restart_yes",2)),

  xcm_orig_(  *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig","comm_none","frame_invariant","restart_yes",3)),
  quat_orig_( *mesh->prop().addGlobalProperty< VectorContainer<double,4> > ("quat_orig","comm_none","frame_invariant","restart_yes",1)),

  xcm_orig_0_(  *mesh->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig_0","comm_none","frame_invariant","restart_yes",3)),
  quat_orig_0_( *mesh->prop().addGlobalProperty< VectorContainer<double,4> > ("quat_orig_0","comm_none","frame_invariant","restart_yes",1)),

  fflag_(    *mesh->prop().addGlobalProperty< VectorContainer<bool,3> > ("fflag","comm_none","frame_invariant","restart_yes",1)),
  tflag_(    *mesh->prop().addGlobalProperty< VectorContainer<bool,3> > ("tflag","comm_none","frame_invariant","restart_yes",1)),

  displace_( *mesh->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("displace","comm_none","frame_invariant","restart_yes",3)),

  v_(        *mesh->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("v","comm_exchange_borders","frame_invariant","restart_yes",1)),

  suspension_flag_(false),
  k_t_(0.),
  c_t_(0.),
  k_r_(0.),
  c_r_(0.),
  externalForce_flag_(false),
  rot_flip_flag_(false),
  rot_flip_angle_(0.),
  fix_gravity_(NULL),
  gravity_set_(false),
  limit_vel_(-1.),
  limit_movement_axis_(false)
{
    mm_stress = static_cast<MeshModuleStress*>(fix_mesh->get_module("stress"));
    if (!mm_stress)
        error->one(FLERR,"Mesh module \"6dof\" requires mesh modules \"stress\"");

    if(!mm_stress->trackStress())
        error->one(FLERR,"stress = 'on' required");

    if(fix_mesh->manipulated())
        error->warning(FLERR,"Mesh has been scaled, moved, or rotated.\n"
                             "Please note that values for 'com', 'vel', 'mass', 'moi' refer to the scaled, moved, or rotated configuration");

    // set defaults

    init_defaults();

    xstr = ystr = zstr = NULL;

    // parse further args
    
    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
      hasargs = false;
      if(strcmp(arg[iarg_],"com") == 0) {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'com'");
          iarg_++;
          double _com[3];
          _com[0] = force->numeric(FLERR,arg[iarg_++]);
          _com[1] = force->numeric(FLERR,arg[iarg_++]);
          _com[2] = force->numeric(FLERR,arg[iarg_++]);
          xcm_.add(_com);
          mm_stress->set_p_ref(xcm_(0));
          hasargs = true;
      }
      else if(strcmp(arg[iarg_],"limit_vel") == 0)
      {
          if (narg < iarg_+2)
              error->one(FLERR,"not enough arguments for 'limit_vel'");
          limit_vel_ = force->numeric(FLERR, arg[iarg_+1]);
          iarg_ += 2;
          hasargs = true;
      }
      else if(strcmp(arg[iarg_],"movement_axis") == 0)
      {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'movement_axis'");
          iarg_++;
          axis_[0] = force->numeric(FLERR,arg[iarg_++]);
          axis_[1] = force->numeric(FLERR,arg[iarg_++]);
          axis_[2] = force->numeric(FLERR,arg[iarg_++]);
          vectorNormalize3D(axis_);
          limit_movement_axis_ = true;
          hasargs = true;
      }
      else if(strcmp(arg[iarg_],"vel") == 0) {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'vel'");
          iarg_++;
          double _vel[3];
          _vel[0] = force->numeric(FLERR,arg[iarg_++]);
          _vel[1] = force->numeric(FLERR,arg[iarg_++]);
          _vel[2] = force->numeric(FLERR,arg[iarg_++]);
          vcm_.set(0,_vel);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"angmom") == 0) {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'angmom'");
          iarg_++;
          double _angmom[3];
          _angmom[0] = force->numeric(FLERR,arg[iarg_++]);
          _angmom[1] = force->numeric(FLERR,arg[iarg_++]);
          _angmom[2] = force->numeric(FLERR,arg[iarg_++]);
          angmom_.set(0,_angmom);
          hasargs = true;
      } else if (strcmp(arg[iarg_],"mass") == 0) {
          if (narg < iarg_+2) error->one(FLERR,"not enough arguments");
          iarg_++;
          double _mass = force->numeric(FLERR,arg[iarg_++]);
          if(_mass <= 0.)
            error->one(FLERR,"mass > 0 required");
          mass_.add(_mass);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"moi") == 0) {
          if (narg < iarg_+7) error->one(FLERR,"not enough arguments for 'moi'");
          iarg_++;
          double **_moi;
          memory->create<double>(_moi,3,3,"6dof:_moi");
          _moi[0][0] = force->numeric(FLERR,arg[iarg_++]);
          _moi[1][1] = force->numeric(FLERR,arg[iarg_++]);
          _moi[2][2] = force->numeric(FLERR,arg[iarg_++]);
          _moi[0][1] = force->numeric(FLERR,arg[iarg_++]);
          _moi[0][2] = force->numeric(FLERR,arg[iarg_++]);
          _moi[1][2] = force->numeric(FLERR,arg[iarg_++]);
          _moi[1][0] = _moi[0][1];
          _moi[2][0] = _moi[0][2];
          _moi[2][1] = _moi[1][2];
          moi_.add(_moi);
          memory->destroy<double>(_moi);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"forceflags") == 0) {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'forceflags'");
          iarg_++;
          bool flags[3];
          if(strcmp("0",arg[iarg_]) == 0)       flags[0] = false;
          else if(strcmp("1",arg[iarg_]) == 0)  flags[0] = true;
          else error->one(FLERR,"wrong arguments for 'forceflags', needs 0 or 1");
          if(strcmp("0",arg[iarg_+1]) == 0)       flags[1] = false;
          else if(strcmp("1",arg[iarg_+1]) == 0)  flags[1] = true;
          else error->one(FLERR,"wrong arguments for 'forceflags', needs 0 or 1");
          if(strcmp("0",arg[iarg_+2]) == 0)       flags[2] = false;
          else if(strcmp("1",arg[iarg_+2]) == 0)  flags[2] = true;
          else error->one(FLERR,"wrong arguments for 'forceflags', needs 0 or 1");
          fflag_.set(0,flags);
          iarg_ = iarg_+3;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"torqueflags") == 0) {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'torqueflags'");
          iarg_++;
          bool flags[3];
          if(strcmp("0",arg[iarg_]) == 0)       flags[0] = false;
          else if(strcmp("1",arg[iarg_]) == 0)  flags[0] = true;
          else error->one(FLERR,"wrong arguments for 'torqueflags', needs 0 or 1");
          if(strcmp("0",arg[iarg_+1]) == 0)       flags[1] = false;
          else if(strcmp("1",arg[iarg_+1]) == 0)  flags[1] = true;
          else error->one(FLERR,"wrong arguments for 'torqueflags', needs 0 or 1");
          if(strcmp("0",arg[iarg_+2]) == 0)       flags[2] = false;
          else if(strcmp("1",arg[iarg_+2]) == 0)  flags[2] = true;
          else error->one(FLERR,"wrong arguments for 'torqueflags', needs 0 or 1");
          tflag_.set(0,flags);
          iarg_ = iarg_+3;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"rot_flip") == 0) {
          if (narg < iarg_+3) error->one(FLERR,"not enough arguments for 'rot_flip'");
          iarg_++;
          rot_flip_flag_ = true;
          if(strcmp(arg[iarg_++], "angle"))
            error->one(FLERR,"expecting keyword 'angle'");
          else
            rot_flip_angle_ = force->numeric(FLERR,arg[iarg_++]);
          if(rot_flip_angle_ < 0. || rot_flip_angle_ > 90.)
            error->one(FLERR,"0° < angle < 90° required");
          rot_flip_angle_ = rot_flip_angle_ * M_PI / 180.;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"suspension") == 0) {
          if (narg < iarg_+9) error->one(FLERR,"not enough arguments for 'suspension'");
          iarg_++;
          suspension_flag_ = true;
          if(strcmp(arg[iarg_++], "k_t"))
            error->one(FLERR,"expecting keyword 'k_t'");
          k_t_ = force->numeric(FLERR,arg[iarg_++]);
          if(k_t_ < 0.)
            error->one(FLERR,"k_t >= 0 required");
          if(strcmp(arg[iarg_++], "c_t"))
            error->one(FLERR,"expecting keyword 'c_t'");
          c_t_ = force->numeric(FLERR,arg[iarg_++]);
          if(c_t_ < 0.)
            error->one(FLERR,"c_t >= 0 required");
          if(strcmp(arg[iarg_++], "k_r"))
            error->one(FLERR,"expecting keyword 'k_r'");
          k_r_ = force->numeric(FLERR,arg[iarg_++]);
          if(k_r_ < 0.)
            error->one(FLERR,"k_r >= 0 required");
          if(strcmp(arg[iarg_++], "c_r"))
            error->one(FLERR,"expecting keyword 'c_r'");
          c_r_ = force->numeric(FLERR,arg[iarg_++]);
          if(c_r_ < 0.)
            error->one(FLERR,"c_r >= 0 required");
          hasargs = true;
      } else if(strcmp(arg[iarg_],"externalForce") == 0) {
          if (narg < iarg_+4) error->one(FLERR,"not enough arguments for 'externalForce'");
          iarg_++;
          externalForce_flag_ = true;
          // check if forces are given as variables:
          if(strstr(arg[iarg_],"v_") == arg[iarg_]){
              int n = strlen(&arg[iarg_][2]) + 1;
              xstr = new char[n];
              strcpy(xstr,&arg[iarg_][2]);
              iarg_++;
          } else {
              xvalue = force->numeric(FLERR,arg[iarg_++]);
              xstyle = CONSTANT;
          }
          if(strstr(arg[iarg_],"v_") == arg[iarg_]){
              int n = strlen(&arg[iarg_][2]) + 1;
              ystr = new char[n];
              strcpy(ystr,&arg[iarg_][2]);
              iarg_++;
          } else {
              yvalue = force->numeric(FLERR,arg[iarg_++]);
              ystyle = CONSTANT;
          }
          if(strstr(arg[iarg_],"v_") == arg[iarg_]){
              int n = strlen(&arg[iarg_][2]) + 1;
              zstr = new char[n];
              strcpy(zstr,&arg[iarg_][2]);
              iarg_++;
          } else {
              zvalue = force->numeric(FLERR,arg[iarg_++]);
              zstyle = CONSTANT;
          }
          hasargs = true;
      }
    }

    error_checks();

    // store original position and rotation state
    xcm_orig_.add(xcm_(0));
    quat_orig_.add(quat_(0));
    xcm_orig_0_.add(xcm_(0));
    quat_orig_0_.add(quat_(0));
}

/* ---------------------------------------------------------------------- */

MeshModuleStress6DOF::~MeshModuleStress6DOF()
{
    delete [] xstr;
    delete [] ystr;
    delete [] zstr;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::post_create()
{
    
    init_rotation_props();

    calc_displace();

    mesh->registerMove(false,true,true);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::init_defaults()
{
    bool truevec[] ={true,true,true};
    fflag_.add(truevec);
    tflag_.add(truevec);
    double unitquat[4];
    quatIdentity4D(unitquat);
    quat_.add(unitquat);
    double zerovec[3];
    vectorZeroize3D(zerovec);
    vcm_.add(zerovec);
    angmom_.add(zerovec);
    omega_.add(zerovec);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::error_checks()
{
    
    if(!xcm_.size())
        error->one(FLERR,"please define 'com' for the 6dof mesh");
    if(!mass_.size())
        error->one(FLERR,"please define 'mass' for the 6dof mesh");
    if(!moi_.size())
        error->one(FLERR,"please define 'moi' for the 6dof mesh");

    if(mesh->nMove() > 1)
      error->one(FLERR,"this fix does not allow superposition with moving mesh fixes");
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::init_rotation_props()
{
  double **evectors;
  memory->create(evectors,3,3,"MeshModuleStress6DOF:evectors");
  
  int ierror = MathExtra::jacobi(moi_(0),inertia_(0),evectors);
  if (ierror) error->one(FLERR,"Insufficient Jacobi rotations for rigid body");

  ex_space_(0)[0] = evectors[0][0];  ex_space_(0)[1] = evectors[1][0];  ex_space_(0)[2] = evectors[2][0];
  ey_space_(0)[0] = evectors[0][1];  ey_space_(0)[1] = evectors[1][1];  ey_space_(0)[2] = evectors[2][1];
  ez_space_(0)[0] = evectors[0][2];  ez_space_(0)[1] = evectors[1][2];  ez_space_(0)[2] = evectors[2][2];

  // if any principal moment < scaled EPSILON, set to 0.0
  
  double max;
  max = std::max(inertia_(0)[0],inertia_(0)[1]);
  max = std::max(max,inertia_(0)[2]);

  if (inertia_(0)[0] < EPSILON*max) inertia_(0)[0] = 0.0;
  if (inertia_(0)[1] < EPSILON*max) inertia_(0)[1] = 0.0;
  if (inertia_(0)[2] < EPSILON*max) inertia_(0)[2] = 0.0;

  double ez0 = ex_space_(0)[1]*ey_space_(0)[2] -ex_space_(0)[2]*ey_space_(0)[1];
  double ez1 = ex_space_(0)[2]*ey_space_(0)[0] -ex_space_(0)[0]*ey_space_(0)[2];
  double ez2 = ex_space_(0)[0]*ey_space_(0)[1] -ex_space_(0)[1]*ey_space_(0)[0];

  if (ez0*ez_space_(0)[0] + ez1*ez_space_(0)[1] + ez2*ez_space_(0)[2] < 0.0)
  {
      ez_space_(0)[0] = -ez_space_(0)[0];
      ez_space_(0)[1] = -ez_space_(0)[1];
      ez_space_(0)[2] = -ez_space_(0)[2];
  }

  // create initial quaternion
  MathExtra::exyz_to_q(ex_space_(0),ey_space_(0),ez_space_(0),quat_(0));

  memory->destroy(evectors);

}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::calc_displace()
{
    
    int nall = mesh->size();
    int nnodes = mesh->numNodes();
    double nodepos[3],tmp[3];

    double **dsplc;
    memory->create<double>(dsplc,nnodes,3,"6dof:dsplc");

    for(int i = 0; i < nall; i++)
    {
        for(int j = 0; j < nnodes; j++)
        {
            mesh->node_slow(i,j,nodepos);
            vectorSubtract3D(nodepos,xcm_(0),tmp);
            dsplc[j][0] = vectorDot3D(tmp,ex_space_(0));
            dsplc[j][1] = vectorDot3D(tmp,ey_space_(0));
            dsplc[j][2] = vectorDot3D(tmp,ez_space_(0));
        }
        displace_.set(i,dsplc);
    }

    memory->destroy<double>(dsplc);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::init()
{
    dtv_ = update->dt;
    dtf_ = 0.5 * update->dt * force->ftm2v;
    dtfm_ = dtf_ / mass_(0);

    dtq_ = 0.5 * update->dt;

    if (strcmp(update->integrate_style,"respa") == 0)
        error->one(FLERR,"not respa-compatible");

    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0)
        error->all(FLERR,"Variable name for fix addforce does not exist");
      if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
      else error->all(FLERR,"Variable for fix addforce is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0)
        error->all(FLERR,"Variable name for fix addforce does not exist");
      if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
      else error->all(FLERR,"Variable for fix addforce is invalid style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0)
        error->all(FLERR,"Variable name for fix addforce does not exist");
      if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
      else error->all(FLERR,"Variable for fix addforce is invalid style");
    }

    const int nfix_gravity = modify->n_fixes_style_strict("gravity");
    if(nfix_gravity == 1)
        fix_gravity_ = static_cast<FixGravity*>(modify->find_fix_style_strict("gravity",0));
    else if (nfix_gravity > 1)
        error->all(FLERR, "mesh module 6dof can only handle one gravity fix");
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::setup_pre_force(int vflag)
{
    mm_stress->set_p_ref(xcm_(0));
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::setup(int vflag)
{
    vectorCopy3D(xcm_(0),xcm_orig_(0));
    vectorCopy4D(quat_(0),quat_orig_(0));
}

/* ---------------------------------------------------------------------- */

int MeshModuleStress6DOF::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::initial_integrate(int vflag)
{
    // double dX[3],dx[3], qOld[4],dQ[4], dq[4];

    if (!gravity_set_)
        add_gravity();

    // update vcm by 1/2 step

    if(fflag_(0)[0]) vcm_(0)[0] += dtfm_ * mm_stress->f_total(0);
    if(fflag_(0)[1]) vcm_(0)[1] += dtfm_ * mm_stress->f_total(1);
    if(fflag_(0)[2]) vcm_(0)[2] += dtfm_ * mm_stress->f_total(2);

    if (limit_movement_axis_)
    {
        const double dot = vectorDot3D(vcm_(0), axis_);
        vectorScalarMult3D(axis_, dot, vcm_(0));
    }

    const double vel = vectorMag3D(vcm_(0));
    if (limit_vel_ > 0 && vel > limit_vel_)
        vectorScalarMult3D(vcm_(0), limit_vel_/vel);

    // update xcm by full step

    dx_[0] = dtv_ * vcm_(0)[0];
    dx_[1] = dtv_ * vcm_(0)[1];
    dx_[2] = dtv_ * vcm_(0)[2];
    vectorAdd3D(xcm_(0),dx_,xcm_(0));
    vectorSubtract3D(xcm_(0),xcm_orig_(0),dX_);

    // update angular momentum by 1/2 step

    if(tflag_(0)[0]) angmom_(0)[0] += dtf_ * mm_stress->torque_total(0);
    if(tflag_(0)[1]) angmom_(0)[1] += dtf_ * mm_stress->torque_total(1);
    if(tflag_(0)[2]) angmom_(0)[2] += dtf_ * mm_stress->torque_total(2);

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    vectorCopy4D(quat_(0),qOld_);
    MathExtra::angmom_to_omega(angmom_(0),ex_space_(0),ey_space_(0),ez_space_(0),inertia_(0),omega_(0));
    MathExtra::richardson(quat_(0),angmom_(0),omega_(0),inertia_(0),dtq_);
    MathExtra::q_to_exyz(quat_(0),ex_space_(0),ey_space_(0),ez_space_(0));

    // calculate quaternion difference
    
    MathExtraLiggghts::quat_diff(quat_(0),qOld_,dq_);
    MathExtraLiggghts::quat_diff(quat_(0),quat_orig_(0),dQ_);

    // move and rotate mesh, set velocity

    mesh->move(dX_,dx_);
    mesh->rotate(dQ_,dq_,xcm_(0));
    set_vel();

    // update reference point to COM
    
    mm_stress->set_p_ref(xcm_(0));
    
    gravity_set_ = false;

}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::final_integrate()
{
    // add extra forces

    add_gravity();

    if(suspension_flag_)
        add_suspension_force();

    if(externalForce_flag_)
        add_external_force();

    if(rot_flip_flag_)
        rot_flip();

    // update vcm by 1/2 step

    if(fflag_(0)[0]) vcm_(0)[0] += dtfm_ * mm_stress->f_total(0);
    if(fflag_(0)[1]) vcm_(0)[1] += dtfm_ * mm_stress->f_total(1);
    if(fflag_(0)[2]) vcm_(0)[2] += dtfm_ * mm_stress->f_total(2);

    if (limit_movement_axis_)
    {
        const double dot = vectorDot3D(vcm_(0), axis_);
        vectorScalarMult3D(axis_, dot, vcm_(0));
    }

    const double vel = vectorMag3D(vcm_(0));
    if (limit_vel_ > 0. && vel > limit_vel_)
        vectorScalarMult3D(vcm_(0), limit_vel_/vel);

    // update angular momentum by 1/2 step

    if(tflag_(0)[0]) angmom_(0)[0] += dtf_ * mm_stress->torque_total(0);
    if(tflag_(0)[1]) angmom_(0)[1] += dtf_ * mm_stress->torque_total(1);
    if(tflag_(0)[2]) angmom_(0)[2] += dtf_ * mm_stress->torque_total(2);

    MathExtra::angmom_to_omega(angmom_(0),ex_space_(0),ey_space_(0),ez_space_(0),inertia_(0),omega_(0));

    set_vel();
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::add_gravity()
{
    double gravity[3],f_grav[3];

    vectorZeroize3D(gravity);
    if (fix_gravity_)
        fix_gravity_->get_gravity(gravity);

    vectorScalarMult3D(gravity,mass_(0),f_grav);
    
    mm_stress->add_global_external_contribution(f_grav);
    gravity_set_ = true;
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::rot_flip()
{
    double dangle;

    dangle = 2.*acos(quat_(0)[0]) - rot_flip_angle_;

    if(dangle > 0.)
    {
        vectorFlip3D(vcm_(0));
        vectorFlip3D(angmom_(0));
        vectorFlip3D(omega_(0));
    }
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::add_suspension_force()
{
    double dX[3], dQ[4], axis[3], angle, force[3], torque[3];
    double force_spring[3], force_damper[3], torque_spring[3], torque_damper[3];

    vectorSubtract3D(xcm_(0),xcm_orig_0_(0),dX);
    MathExtraLiggghts::quat_diff(quat_(0),quat_orig_0_(0),dQ);

    // force components

    vectorScalarMult3D(dX,-k_t_,force_spring);
    vectorScalarMult3D(vcm_(0),-c_t_,force_damper);

    // torque components
    
    angle = 2.*acos(dQ[0]);
    MathExtraLiggghts::vec_from_quat(dQ,axis);
    vectorNormalize3D(axis);
    vectorScalarMult3D(axis,-k_r_*angle,torque_spring);
    vectorScalarMult3D(omega_(0),-c_r_,torque_damper);

    // add forces and contribute them

    vectorAdd3D(force_spring,force_damper,force);
    vectorAdd3D(torque_spring,torque_damper,torque);
    
    mm_stress->add_global_external_contribution(force,torque);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::add_external_force()
{
    double force[3];

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);

    modify->addstep_compute(update->ntimestep + 1);

    if(fflag_(0)[0]) force[0] = xvalue;
    if(fflag_(0)[1]) force[1] = yvalue;
    if(fflag_(0)[2]) force[2] = zvalue;

    mm_stress->add_global_external_contribution(force);
}

/* ---------------------------------------------------------------------- */

void MeshModuleStress6DOF::set_vel()
{
    const int nall = mesh->size();
    const int nnodes = mesh->numNodes();
    double dx,dy,dz;
    double **vNodes;

    memory->create<double>(vNodes,nnodes,3,"6dof:vNodes");

    for(int i = 0; i < nall; i++)
    {
        for(int j = 0; j < nnodes; j++)
        {
            dx = ex_space_(0)[0]*displace_(i)[j][0] + ey_space_(0)[0]*displace_(i)[j][1] +  ez_space_(0)[0]*displace_(i)[j][2];
            dy = ex_space_(0)[1]*displace_(i)[j][0] + ey_space_(0)[1]*displace_(i)[j][1] +  ez_space_(0)[1]*displace_(i)[j][2];
            dz = ex_space_(0)[2]*displace_(i)[j][0] + ey_space_(0)[2]*displace_(i)[j][1] +  ez_space_(0)[2]*displace_(i)[j][2];
            vNodes[j][0] = omega_(0)[1]*dz - omega_(0)[2]*dy + vcm_(0)[0];
            vNodes[j][1] = omega_(0)[2]*dx - omega_(0)[0]*dz + vcm_(0)[1];
            vNodes[j][2] = omega_(0)[0]*dy - omega_(0)[1]*dx + vcm_(0)[2];
        }
        v_.set(i,vNodes);
    }

    memory->destroy<double>(vNodes);
}

/* ----------------------------------------------------------------------
   Compute eulerian angle from quaternion
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   return total force or torque component on body or xcm
------------------------------------------------------------------------- */

double MeshModuleStress6DOF::compute_vector(int n)
{
    double dQ[4];
    MathExtraLiggghts::quat_diff(quat_(0),quat_orig_(0),dQ);
    double myOmega[3];
    MathExtraLiggghts::quat_to_EulerAngle(dQ,myOmega);

    myOmega[0] = myOmega[0]*180/M_PI;
    myOmega[1] = myOmega[1]*180/M_PI;
    myOmega[2] = myOmega[2]*180/M_PI;

    if (n < 3)
        return xcm_(0)[n];
    else if (n < 7)
        return quat_(0)[n-3];
    else 
        return myOmega[n-7];

}

/* ---------------------------------------------------------------------- */

int MeshModuleStress6DOF::modify_param(int narg, char **arg)
{
    std::string arg0(arg[0]);
    std::size_t slash = arg0.find_first_of('/');
    std::string command(arg0.substr(slash == std::string::npos ? 0 : slash+1));

    if (command.compare("forceflags") == 0)
    {
        if (narg < 4)
            error->one(FLERR,"not enough arguments for fix_modify 'forceflags'");
        bool flags[3];
        if(strcmp("0",arg[1]) == 0)
            flags[0] = false;
        else if(strcmp("1",arg[1]) == 0)
            flags[0] = true;
        else
            error->one(FLERR,"wrong arguments for 'forceflags', needs 0 or 1");

        if(strcmp("0",arg[2]) == 0)
            flags[1] = false;
        else if(strcmp("1",arg[2]) == 0)
            flags[1] = true;
        else
            error->one(FLERR,"wrong arguments for 'forceflags', needs 0 or 1");

        if(strcmp("0",arg[3]) == 0)
            flags[2] = false;
        else if(strcmp("1",arg[3]) == 0)
            flags[2] = true;
        else
            error->one(FLERR,"wrong arguments for 'forceflags', needs 0 or 1");

        fflag_.set(0,flags);
        return 4;
    }
    else if(command.compare("torqueflags") == 0)
    {
        if (narg < 4)
            error->one(FLERR,"not enough arguments for fix_modify 'torqueflags'");
        bool flags[3];
        if(strcmp("0",arg[1]) == 0)
            flags[0] = false;
        else if(strcmp("1",arg[1]) == 0)
            flags[0] = true;
        else
            error->one(FLERR,"wrong arguments for 'torqueflags', needs 0 or 1");

        if(strcmp("0",arg[2]) == 0)
            flags[1] = false;
        else if(strcmp("1",arg[2]) == 0)
            flags[1] = true;
        else
            error->one(FLERR,"wrong arguments for 'torqueflags', needs 0 or 1");

        if(strcmp("0",arg[3]) == 0)
            flags[2] = false;
        else if(strcmp("1",arg[3]) == 0)
            flags[2] = true;
        else
            error->one(FLERR,"wrong arguments for 'torqueflags', needs 0 or 1");

        tflag_.set(0,flags);
        return 4;
    }
    else if(command.compare("vel") == 0)
    {
        if (narg < 4)
            error->one(FLERR,"not enough arguments for fix_modify 'vel'");
        double _vel[3];
        _vel[0] = force->numeric(FLERR,arg[1]);
        _vel[1] = force->numeric(FLERR,arg[2]);
        _vel[2] = force->numeric(FLERR,arg[3]);
        vcm_.set(0,_vel);

        return 4;
    }
    else if(command.compare("move_mesh") == 0)
    {
        if (narg < 4)
            error->one(FLERR,"not enough arguments for fix_modify 'move_mesh'");
        double dx[3];
        dx[0] = force->numeric(FLERR,arg[1]);
        dx[1] = force->numeric(FLERR,arg[2]);
        dx[2] = force->numeric(FLERR,arg[3]);
        fix_mesh->moveMesh(dx[0], dx[1], dx[2]);
        vectorAdd3D(xcm_(0),dx,xcm_(0));

        return 4;
    }
    else if(command.compare("angmom") == 0)
    {
        if (narg < 4)
            error->one(FLERR,"not enough arguments for fix_modify 'angmom'");
        double _angmom[3];
        _angmom[0] = force->numeric(FLERR,arg[1]);
        _angmom[1] = force->numeric(FLERR,arg[2]);
        _angmom[2] = force->numeric(FLERR,arg[3]);
        angmom_.set(0,_angmom);

        return 4;
    }
    else if(command.compare("omega") == 0)
    {
        if (narg < 4)
            error->one(FLERR,"not enough arguments for fix_modify 'omega'");
        // read the omega values and convert them to the angmom
        double _omega[3];
        _omega[0] = force->numeric(FLERR,arg[1]);
        _omega[1] = force->numeric(FLERR,arg[2]);
        _omega[2] = force->numeric(FLERR,arg[3]);
        double _angmom[3];
        //MathExtra::angmom_to_omega(_omega,ex_space_(0),ey_space_(0),ez_space_(0),inertia_(0),_angmom);
        MathExtra::omega_to_angmom(_omega,ex_space_(0),ey_space_(0),ez_space_(0),inertia_(0),_angmom);
        angmom_.set(0,_angmom);

        return 4;
    }

    return 0;
}
