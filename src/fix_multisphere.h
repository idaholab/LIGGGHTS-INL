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

#ifdef FIX_CLASS

FixStyle(multisphere,FixMultisphere)
FixStyle(multisphere/nointegration,FixMultisphere)
FixStyle(concave,FixMultisphere)

#else

#ifndef LMP_FIX_MULTISPHERE_H
#define LMP_FIX_MULTISPHERE_H

#include "fix.h"
#include "math_extra_liggghts.h"
#include "multisphere_parallel.h"
#include "fix_property_atom.h"
#include "fix_remove.h"
#include "fix_heat_gran.h"
#include "atom.h"
#include "comm.h"

namespace LAMMPS_NS {

enum
{
    MS_COMM_UNDEFINED,
    MS_COMM_FW_BODY,
    MS_COMM_FW_IMAGE_DISPLACE,
    MS_COMM_FW_V_OMEGA,
    MS_COMM_FW_F_TORQUE,
    MS_COMM_FW_TEMP,
    MS_COMM_REV_X_V_OMEGA,
    MS_COMM_REV_V_OMEGA,
    MS_COMM_REV_IMAGE,
    MS_COMM_REV_DISPLACE,
    MS_COMM_REV_TEMP
};

class FixMultisphere : public Fix
{
    friend class SetMultisphere;
    friend class FixMoveMultisphere;

public:

    FixMultisphere(class LAMMPS *, int, char **);
    virtual ~FixMultisphere();

    void post_create();
    void pre_delete(bool unfixflag);
    virtual int setmask();
    virtual void init();

    virtual void setup(int);
    virtual void setup_pre_force(int);
    virtual void setup_pre_exchange();
    virtual void setup_pre_neighbor();

    virtual double extend_cut_ghost();

    virtual void pre_force(int);
    virtual void pre_final_integrate();
    void comm_correct_force(bool setupflag);
    virtual void calc_force(bool setupflag);

    void add_body_finalize();

    double memory_usage();
    void grow_arrays(int);
    void copy_arrays(int i, int j, int delflag);

    void pre_exchange();
    void pre_neighbor();
    void set_arrays(int);

    int size_restart(int nlocal);
    int maxsize_restart();

    void write_restart(FILE *fp);
    void restart(char *buf, const Version &ver);

    void reset_dt();

    int modify_param(int narg, char **arg);

    // *************************************
    #include "fix_multisphere_comm_I.h"
    // *************************************

    int dof(int);
    double ** get_dump_ref(int &nb, int &nprop, char* prop);
    double max_r_bound();

    void add_remove_callback(FixRemove *ptr);

    void *extract_ms(const char *name, int &len1, int &len2);

    // public inline access

    class Multisphere& data()
    { return multisphere_;}

    class FixPropertyAtom* fix_delflag()
    { return fix_delflag_; }

    int belongs_to(int i) const
    { return body_[i]; }

    int n_body()
    { return data().n_body(); }

    int n_body_all()
    { return data().n_body_all(); }

    int tag_max_body()
    { return data().tag_max_body(); }

    int ntypes()
    { return ntypes_; }

    double* vclump()
    { return Vclump_; }

    double extract_ke()
    { return data().extract_ke(); }

    double extract_rke()
    { return data().extract_rke(); }

    double extract_vave()
    { return data().extract_vave(); }

    double extract_omega_ave()
    { return data().extract_omega_ave(); }

    void set_v_body_from_atom_index(int iatom,double *vel)
    {if(body_[iatom] >= 0) multisphere_.set_v_body(body_[iatom],vel); }

    void set_omega_body_from_atom_index(int iatom,double *omega)
    {if(body_[iatom] >= 0) multisphere_.set_omega_body(body_[iatom],omega); }

    void set_body_displace(int i,double *_displace,int body_id,double volume_weight)
    {
      body_[i] = body_id; vectorCopy3D(_displace,displace_[i]);
      if(fix_volumeweight_ms_)
          fix_volumeweight_ms_->vector_atom[i] = volume_weight;
    }

    int calc_n_steps(int iatom,double *p_ref,double *normalvec,double *v_normal)
    { return multisphere_.calc_n_steps(iatom,body_[iatom],p_ref,normalvec,v_normal); }

    void release(int iatom,double *v_toInsert,double *omega_toInsert)
    { return multisphere_.release(iatom,body_[iatom],v_toInsert,omega_toInsert); }

    bool allow_group_and_set()
    { return allow_group_and_set_; }

    bool use_volumeweight()
    { return use_volumeweight_ms_; }

    const FixPropertyAtom *get_volumeweight() const
    { return fix_volumeweight_ms_; }

    void scale_displace(int i, double factor)
    { vectorScalarMult3D(displace_[i],factor); }

    void rev_comm_displace()
    {
      rev_comm_flag_ = MS_COMM_REV_DISPLACE;
      reverse_comm();
    }

    void rev_comm(int flag)
    {
      rev_comm_flag_ = flag;
      reverse_comm();
    }

    void fw_comm(int flag)
    {
      fw_comm_flag_ = flag;
      forward_comm();
    }

    class ParallelBase* get_parallel_base_ptr() const
    { return &multisphere_; }

    void set_xv();
    void set_xv(int);
    void set_v();
    void set_v(int);

 protected:

    void ms_error(const char * file, int line,char const *errmsg);

    int map(int i)
    { return data().map(i); }

    int tag(int i)
    { return data().tag(i); }

    bool do_modify_body_forces_torques_;
    virtual void modify_body_forces_torques() {}

    class Multisphere &multisphere_;
    class FixPropertyAtom *fix_corner_ghost_;
    class FixPropertyAtom *fix_delflag_;
    class FixPropertyAtom *fix_existflag_;
    class FixPropertyAtom *fix_volumeweight_ms_; 
    bool use_volumeweight_ms_;
    class FixGravity *fix_gravity_;
    FixHeatGran *fix_heat_;

    //int comm_di_;
    int fw_comm_flag_;
    int rev_comm_flag_;

    // per-atom properties handled by this fix
    int *body_;                // which body each atom is part of (-1 if none)
    double **displace_;        // displacement of each atom in body coords

    double dtv,dtf,dtq;

    std::vector<FixRemove*> fix_remove_;

    // MS communication
    int ntypes_;
    double *Vclump_;

    bool allow_group_and_set_;
    bool allow_heatsource_;

    bool concave_;    

    int getMask(int ibody)
    { return 1; }

};

}

#endif
#endif
