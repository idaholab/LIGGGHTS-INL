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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MULTISPHERE_H
#define LMP_MULTISPHERE_H

#include "pointers.h"
#include "custom_value_tracker.h"
#include "mpi_liggghts.h"
#include "update.h"
#include "math_extra.h"
#include "parallel_base.h"
#include <vector>
#include <cmath>

namespace LAMMPS_NS {

  class Multisphere : protected Pointers, public ParallelBase {

    friend class FixMultisphere;
    friend class FixChangeSizeMultisphere;
    friend class SetMultisphere;
    friend class FixMoveMultisphere;

    public:

      void add_body(int nspheres, double *xcm_ins, double *xcm_to_xbound_ins,
                    double r_bound_ins, double *v_ins, double *omega_ins,
                    double mass_ins, double dens_ins, int atomtype_ins, int type_ins,
                    double *inertia_ins, double *ex_space_ins, double *ey_space_ins, double *ez_space_ins,
                    double **displace_ins,bool *fflag, bool *tflag, int start_step_ins = -1, double *v_integrate_ins = NULL);

      void grow_arrays_per_body_local(int);
      void grow_arrays_per_body_global(int);

      void remove_body(int ilocal);
      void copy_body(int from_local, int to_local);

      void remap_bodies(int *body);

      void clear_map();
      void generate_map();
      void id_extend_body_extend(int *body);

      void calc_nbody_all();
      bool check_lost_atoms(int *body, double *atom_delflag,double *body_existflag, double *volumeweight);

      int calc_n_steps(int iatom,int body,double *p_ref,double *normalvec,double *v_normal);
      void recalc_n_steps(double dt_ratio);
      void release(int iatom,int body,double *v_toInsert,double *omega_toInsert);

      double max_r_bound();

      virtual void exchange() {}
      virtual void restart_serial(double *);
      virtual void write_restart_parallel(FILE *fp);
      virtual void restart_parallel(double *list);

      void reset_forces(bool extflag);

      void* extract(const char *name, int &, int &);

      double *extract_double_scalar(const char *name);
      double **extract_double_vector(const char *name);

      double extract_ke();
      double extract_rke();
      double extract_vave();
      double extract_omega_ave();

      // inline access functions

      inline int n_body() const
      { return nbody_; }

      inline int n_body_all() const
      { return nbody_all_; }

      inline int tag_max_body() const
      { return mapTagMax_; }

      inline int map(int ibody_local) const
      { return mapArray_?mapArray_[ibody_local]:-1; }

      inline int tag(int ibody_local) const
      { return id_(ibody_local); }

      inline bool has_tag(int _tag) const
      { return mapArray_[_tag] == -1 ? false : true;}

      inline int atomtype(int ibody_local) const
      { return atomtype_(ibody_local); }

      inline int nrigid(int ibody_local) const
      { return nrigid_(ibody_local); }

      inline void xcm(double *x_cm,int ibody_local) const
      { vectorCopy3D(xcm_(ibody_local),x_cm); }

      inline void vcm(double *v_cm,int ibody_local) const
      { vectorCopy3D(vcm_(ibody_local),v_cm); }

      inline void omega(double * const omega, const int ibody_local) const
      { vectorCopy3D(omega_(ibody_local), omega); }

      inline void angmom(double * const angmom, const int ibody_local) const
      { vectorCopy3D(angmom_(ibody_local), angmom); }

      inline void quat(double *quat,int ibody_local) const
      { vectorCopy4D(quat_(ibody_local),quat); }

      inline void add_external_force(double *frc,int ibody_local)
      { vectorAdd3D(fcm_(ibody_local),frc,fcm_(ibody_local)); }

      inline void x_bound(double *x_bnd,int ibody_local) const
      {
        vectorZeroize3D(x_bnd);
        MathExtraLiggghts::local_coosys_to_cartesian(x_bnd,xcm_to_xbound_(ibody_local),
                            ex_space_(ibody_local),ey_space_(ibody_local),ez_space_(ibody_local));
        vectorAdd3D(xcm_(ibody_local),x_bnd,x_bnd);
      }

      inline double r_bound(int ibody_local) const
      { return r_bound_(ibody_local); }

      inline double mass(int ibody_local) const
      { return masstotal_(ibody_local); }

      inline double density(int ibody_local) const
      { return density_(ibody_local); }

      inline double volume(int ibody_local) const
      { return masstotal_(ibody_local)/density_(ibody_local); }

      inline void set_v_body(int ibody_local,double *vel)
      { vcm_.set(ibody_local,vel); }

      inline void set_omega_body(int ibody_local,double *omega)
      { omega_.set(ibody_local,omega); }

      inline void set_angmom_via_omega_body(int ibody_local,double *omega)
      {
        double angmom[3];
        MathExtra::omega_to_angmom(omega, ex_space_(ibody_local), ey_space_(ibody_local), ez_space_(ibody_local), inertia_(ibody_local), angmom);
        omega_.set(ibody_local,omega);
        angmom_.set(ibody_local,angmom);
      }

      inline void set_angmom_body(int ibody_local,double *angmom)
      { angmom_.set(ibody_local,angmom); }

      void set_fflag(int ibody_local,bool *fflag)
      { fflag_.set(ibody_local, fflag); }

      void set_tflag(int ibody_local,bool *tflag)
      { tflag_.set(ibody_local, tflag); }

      inline class CustomValueTracker& prop()
      { return customValues_; }

      void tagReset(int ibody_local)
      { id_(ibody_local)=-1; }

      void nrigidReset(int ibody_local, int resetValue)
      { nrigid_(ibody_local) = resetValue; }

      // implementation of the ParallelBase interface
      
      void clear_exchange() {}
      int pack_exchange(const int /*dim*/) { return 0; }
      void unpack_exchange(const int /*nrecv*/, const int /*dim*/) {}
      void post_exchange() {}

      int get_ntotal() const
      { return nbody_; } // multisphere has no ghosts
      int get_nlocal() const
      { return nbody_; }
      double get_pos(const int i, const int dim) const
      { return xcm_(i)[dim]; }
      bool has_radius() const
      { return true; }
      double get_radius(const int i) const
      { return r_bound_(i); }
      int get_type(const int i) const
      { return type_(i); }

      double *const *const get_xcm_ptr()
      { return xcm_.begin(); }
      double *const *const get_vcm_ptr()
      { return vcm_.begin(); }
      double *const *const get_fcm_ptr()
      { return fcm_.begin(); }
      double *const *const get_torquecm_ptr()
      { return torquecm_.begin(); }
      double *const *const get_ex_space_ptr()
      { return ex_space_.begin(); }
      double *const *const get_ey_space_ptr()
      { return ey_space_.begin(); }
      double *const *const get_ez_space_ptr()
      { return ez_space_.begin(); }
      double *const *const get_angmom_ptr()
      { return angmom_.begin(); }
      double *const *const get_omega_ptr()
      { return omega_.begin(); }
      double *const *const get_quat_ptr()
      { return quat_.begin(); }
      double *const *const get_inertia_ptr()
      { return inertia_.begin(); }
      double *const get_masstotal_ptr()
      { return masstotal_.begin(); }
      double *const get_density_ptr()
      { return density_.begin(); }
      int *const get_start_step_ptr()
      { return start_step_.begin(); }
      double *const *const get_v_integrate_ptr()
      { return v_integrate_.begin(); }
      bool *const *const get_fflag_ptr()
      { return fflag_.begin(); }
      bool *const *const get_tflag_ptr()
      { return tflag_.begin(); }

      virtual int size_restart() const;
      virtual int pack_restart(double *const buf) const;
      virtual int unpack_restart(double *const buf);
      virtual void finalize_restart();

    protected:

      Multisphere(LAMMPS *lmp);
      virtual ~Multisphere();

      // class holding fields
      CustomValueTracker &customValues_;

      // # of local rigid bodies, and global # bodies on all procs
      int nbody_, nbody_all_;

      // global-local lookup
      int mapTagMax_;
      int *mapArray_;

      // ID of rigid body
      
      ScalarContainer<int> &id_;

      // basic body properties

      // center-of-mass coords, vels, forces, torques of each rigid body
      // extra (external) force on center-of-mass of each
      VectorContainer<double,3> &xcm_;
      VectorContainer<double,3> &vcm_;
      VectorContainer<double,3> &fcm_;
      VectorContainer<double,3> &torquecm_;

      // angular momentum of each in space coords
      // angular velocity of each in space coords
      // quaternion of each rigid body
      VectorContainer<double,3> &angmom_;
      VectorContainer<double,3> &omega_;
      VectorContainer<double,4> &quat_;

      // density and total mass of each rigid body
      // 3 principal components of inertia of each
      // principal axes of each in space coords
      ScalarContainer<int> &atomtype_; 
      ScalarContainer<int> &type_;     
      ScalarContainer<double> &density_;
      ScalarContainer<double> &masstotal_;
      VectorContainer<double,3> &inertia_;
      VectorContainer<double,3> &ex_space_;
      VectorContainer<double,3> &ey_space_;
      VectorContainer<double,3> &ez_space_;

      // # of atoms in each rigid body
      ScalarContainer<int> &nrigid_;

      // image flags of xcm of each rigid body
      ScalarContainer<int> &imagebody_;
      VectorContainer<int,4> &remapflag_;

      // flag for on/off of center-of-mass force, torque
      VectorContainer<bool,3> &fflag_;
      VectorContainer<bool,3> &tflag_;

      // step to start from for integration
      ScalarContainer<int> &start_step_;
      VectorContainer<double,3> &v_integrate_;

      // bounding radius for each body
      // vector from xcm to center of bound sphere
      ScalarContainer<double> &r_bound_;
      VectorContainer<double,3> &xcm_to_xbound_;

      // temperature and buffer for each body
      ScalarContainer<double> &temp_;
      ScalarContainer<double> &temp_old_;
  };

  // *************************************
  #include "multisphere_I.h"
  // *************************************

} //Namespace

#endif
