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
    Alexander Podlozhnyuk, DCS Computing GmbH, Linz
    Christoph Kloss, DCS Computing GmbH, Linz

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <string.h>
#include "fix_nve_non_spherical.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h"
#include "math_extra_liggghts_nonspherical.h"
#include "fix_wall_gran.h"
#include "fix_contact_property_atom_wall.h"
#include "fix_cfd_coupling_force.h"

enum shapes {SUPERQUADRIC, MULTISPHERE, CONCAVE, CONVEX};
enum solver {RICHARDSON, SYMPLECTIC, PREDICTOR_CORRECTOR, WOODEM};

/* ---------------------------------------------------------------------- */

FixNVENonSpherical::FixNVENonSpherical(LAMMPS *lmp, int narg, char **arg) :
    FixNVE(lmp, narg, arg),
    concave_(false),
    multisphere_(0),
    fix_multisphere_(NULL),
    CAddRhoFluid_(0.),
    particle_type(-1),
    integrate_concave(true),
    integrate_superquadric(true),
    integrate_multisphere(true),
    integrate_convexhull(true)
{
    if (narg < 3)
        error->all(FLERR,"Illegal fix nve/nonspherical command");

    if (strcmp(style, "nve/convexhull") == 0)
        error->warning(FLERR, "Fix nve/convexhull no longer exists. For backward compatibility this fix was replaced by a fix nve/nonspherical, please update your input script as this will no longer work in future versions");
    if (strcmp(style, "nve/superquadric") == 0)
        error->warning(FLERR, "Fix nve/superquadric no longer exists. For backward compatibility this fix was replaced by a fix nve/nonspherical, please update your input script as this will no longer work in future versions");

    time_integrate = 1;

    // process extra keywords

    integration_scheme = SYMPLECTIC;

    int iarg = 3;
    while (iarg < narg)
    {
        if (strcmp(arg[iarg],"integration_scheme") == 0)
        {
            int integration_scheme_ = force->numeric(FLERR,arg[iarg+1]);
            if(integration_scheme_ == 0)
                integration_scheme = RICHARDSON;
            else if(integration_scheme_ == 1)
                integration_scheme = SYMPLECTIC;
            else if(integration_scheme_ == 2)
                integration_scheme = PREDICTOR_CORRECTOR;
            else if(integration_scheme_ == 3)
                integration_scheme = WOODEM;
            else
                error->one(FLERR,"Invalid integration scheme! Please choose between 0, 1 (default), 2 or 3!");
            iarg += 2;
        }
        else if (strcmp(arg[iarg],"integrate_superquadric") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'integrate_superquadric'");
            iarg ++;
            if (strcmp(arg[iarg],"yes") == 0)
              integrate_superquadric = true;
            else if (strcmp(arg[iarg],"no") == 0)
              integrate_superquadric = false;
            else
                error->one(FLERR,"Expected yes or no for 'integrate_superquadric'");
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_convexhull") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'integrate_convexhull'");
            iarg ++;
            if (strcmp(arg[iarg],"yes") == 0)
              integrate_convexhull = true;
            else if (strcmp(arg[iarg],"no") == 0)
              integrate_convexhull = false;
            else
                error->one(FLERR,"Expected yes or no for 'integrate_convexhull'");
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_multisphere") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'integrate_multisphere'");
            iarg ++;
            if (strcmp(arg[iarg],"yes") == 0)
              integrate_multisphere = true;
            else if (strcmp(arg[iarg],"no") == 0)
              integrate_multisphere = false;
            else
                error->one(FLERR,"Expected yes or no for 'integrate_multisphere'");
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_concave") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'integrate_concave'");
            iarg ++;
            if (strcmp(arg[iarg],"yes") == 0)
              integrate_concave = true;
            else if (strcmp(arg[iarg],"no") == 0)
              integrate_concave = false;
            else
                error->one(FLERR,"Expected yes or no for 'integrate_concave'");
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_superquadric_only") == 0)
        {
            integrate_superquadric = true;
            integrate_convexhull = false;
            integrate_multisphere = false;
            integrate_concave = false;
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_convexhull_only") == 0)
        {
            integrate_superquadric = false;
            integrate_convexhull = true;
            integrate_multisphere = false;
            integrate_concave = false;
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_multisphere_only") == 0)
        {
            integrate_superquadric = false;
            integrate_convexhull = false;
            integrate_multisphere = true;
            integrate_concave = false;
            iarg ++;
        }
        else if (strcmp(arg[iarg],"integrate_concave_only") == 0)
        {
            integrate_superquadric = false;
            integrate_convexhull = false;
            integrate_multisphere = false;
            integrate_concave = true;
            iarg ++;
        }

        else error->fix_error(FLERR,this,"unknown keyword");
    }

    global_freq = 1;
    extarray = 0;

    reset_dt();

}

FixNVENonSpherical::~FixNVENonSpherical()
{
    std::list<ContactForceCorrector*>::iterator it;
    for (it = correctors_.begin(); it != correctors_.end(); it++)
        delete *it;
    correctors_.clear();
}

void FixNVENonSpherical::convexHullInit()
{
    while(correctors_.size() > 0)
    {
        std::list<ContactForceCorrector*>::iterator it;
        for (it = correctors_.begin(); it != correctors_.end(); it++)
            delete *it;
        correctors_.clear();
    }

    int nwalls = modify->n_fixes_style("wall/gran");
    for(int iwall = 0; iwall < nwalls; iwall++)
    {
        FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",iwall));

        if(!fwg->store_force_contact())
            error->fix_error(FLERR,this,"requires to use 'store_force_contact yes' with granular walls");

        if(fwg->is_mesh_wall())
        {
            
            int n_meshes = fwg->n_meshes();
            for(int imesh = 0; imesh < n_meshes; imesh++)
            {
                char fixid[200];
                sprintf(fixid,"contactforces_%s",(fwg->mesh_list())[imesh]->id);

                FixContactPropertyAtomWall *fix_contact_wall = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
                if(!fix_contact_wall)
                    error->fix_error(FLERR,this,"Internal error: need fix contactforces");

                sprintf(fixid,"%s",(fwg->mesh_list())[imesh]->id);
                correctors_.push_back(new ContactForceCorrector(lmp,fix_contact_wall,NULL,NULL,true,fixid));
            }
        }
        else
        {
            
            char fixid[200];
            sprintf(fixid,"contactforces_%s",fwg->id);

            FixContactPropertyAtomWall *fix_contact_wall = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
            if(!fix_contact_wall)
                error->fix_error(FLERR,this,"Internal error: need fix contactforces");

            sprintf(fixid,"%s",fwg->id);
            correctors_.push_back(new ContactForceCorrector(lmp,fix_contact_wall,NULL,NULL,true,fixid));
        }
    }

    if(modify->my_index(this) <  modify->index_last_fix_of_style("wall/gran"))
        error->fix_error(FLERR,this,"needs to be after all fixes of type wall/gran");
}

/* ---------------------------------------------------------------------- */

void FixNVENonSpherical::init()
{
    FixNVE::init();

    fix_multisphere_ = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
    if(fix_multisphere_)
        multisphere_ = & fix_multisphere_->data();

    //    class FixCfdCouplingForce *fix_coupling_force_ = static_cast<FixCfdCouplingForce*>(modify->find_fix_style("couple/cfd/force",0));
    class FixCfdCouplingForce *fix_coupling_force_ = (FixCfdCouplingForce*)modify->find_fix_style("couple/cfd/force",0);
    if (fix_coupling_force_)
    {
        CAddRhoFluid_ = fix_coupling_force_->getCAddRhoFluid();
        if (comm->me == 0) fprintf(screen,"CAddRhoFluid is %e\n",CAddRhoFluid_);
    }

    if(!atom->sphere_flag)
    {
        integrate_multisphere = false;
    }

    if(!atom->superquadric_flag)
    {
        integrate_superquadric = false;
    }

    if(!atom->convex_hull_flag)
    {
        integrate_convexhull = false;
        integrate_concave = false;
    }

    if(!fix_multisphere_)
    {
        integrate_multisphere = false;
        integrate_concave = false;
    }

    if(integrate_convexhull) convexHullInit();

    reset_dt();

    // error checks might go here
}

/* ---------------------------------------------------------------------- */

//dynamic Euler equations for angular velocity in body's principal axes
void FixNVENonSpherical::dynamic_euler(double *wbody, double *tbody, double *inertia, double *result)
{
    result[0] = tbody[0] / inertia[0] + wbody[1]*wbody[2]*((inertia[1] - inertia[2]) / inertia[0]);
    result[1] = tbody[1] / inertia[1] + wbody[2]*wbody[0]*((inertia[2] - inertia[0]) / inertia[1]);
    result[2] = tbody[2] / inertia[2] + wbody[0]*wbody[1]*((inertia[0] - inertia[1]) / inertia[2]);
}

/* ---------------------------------------------------------------------- */

void FixNVENonSpherical::integrate_dynamic_euler(double dt, double *wbody, double *tbody, double *inertia)
{
    double omega_der[3];
    double tol = 1e-12;
    if(LAMMPS_NS::vectorMag3D(wbody)*dt > 1.0)
        error->one(FLERR, "Timestep is too big for rotation integration!");
    dynamic_euler(wbody, tbody, inertia, omega_der);

    double omega_half_prev[3], delta[3];
    double omega_half[] = {0.0, 0.0, 0.0};
    while(1)
    {
        LAMMPS_NS::vectorCopy3D(omega_half, omega_half_prev);
        omega_half[0] = wbody[0] + 0.5*dt*omega_der[0];
        omega_half[1] = wbody[1] + 0.5*dt*omega_der[1];
        omega_half[2] = wbody[2] + 0.5*dt*omega_der[2];
        LAMMPS_NS::vectorSubtract3D(omega_half_prev, omega_half, delta);
        double omega_half_mag = LAMMPS_NS::vectorMag3D(omega_half);
        if(omega_half_mag > 0.0)
        {
            double eps = LAMMPS_NS::vectorMag3D(delta) / omega_half_mag;
            if(eps < tol) break;
            dynamic_euler(omega_half, tbody, inertia, omega_der);
        }
        else break;
    }
    wbody[0] += dt*omega_der[0];
    wbody[1] += dt*omega_der[1];
    wbody[2] += dt*omega_der[2];
}

/* ---------------------------------------------------------------------- */

void FixNVENonSpherical::integrate_quaternion(double *quat, double *wbody)
{
    double q2[4];
    q2[0] = 1.0;
    q2[1] = -wbody[0] * dtq;
    q2[2] = -wbody[1] * dtq;
    q2[3] = -wbody[2] * dtq;
    MathExtraLiggghtsNonspherical::invquat(q2); //using implicit scheme

    double q_temp[4];
    MathExtra::quatquat(quat, q2, q_temp);
    vectorCopy4D(q_temp, quat);
    MathExtra::qnormalize(quat);
}

/* ---------------------------------------------------------------------- */

void FixNVENonSpherical::initial_integrate(int vflag)
{
    if(integrate_convexhull || integrate_superquadric) initial_integrate_single(vflag);
    if(integrate_concave || integrate_multisphere) initial_integrate_clump(vflag);
}

/* ---------------------------------------------------------------------- */
void FixNVENonSpherical::final_integrate()
{
    if(integrate_convexhull || integrate_superquadric) final_integrate_single();
    if(integrate_concave || integrate_multisphere) final_integrate_clump();
}

void FixNVENonSpherical::initial_integrate_single(int vflag)
{
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    double *rmass = atom->rmass;
    double **angmom = atom->angmom;
    double **quat = atom->quaternion;
    double **inertia = atom->inertia;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    double rotation_matrix[9];
    const double dtf2 = 2.0 * dtf;

    for (int i = 0; i < nlocal; i++)
    {
        bool isSingleParticle = true;
        if(fix_multisphere_ && fix_multisphere_->belongs_to(i) >= 0)
          isSingleParticle = false;
        if ((mask[i] & groupbit) && isSingleParticle) // do not integrate if atom belongs to multisphere or concave particle
        {
            MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix);
            const double dtfm = dtf / rmass[i];

            //update velocity by step t+1/2dt
            v[i][0] += dtfm * f[i][0];
            v[i][1] += dtfm * f[i][1];
            v[i][2] += dtfm * f[i][2];

            //update position by step t+dt
            x[i][0] += dtv * v[i][0];
            x[i][1] += dtv * v[i][1];
            x[i][2] += dtv * v[i][2];

            if(integration_scheme == RICHARDSON)
            {
               //update angular moment by step t+1/2dt
                bool tflag[3] = {true, true, true};
                double ex_space[3],ey_space[3],ez_space[3];
                richardson_initial_integrate(dtf, angmom[i], omega[i], torque[i], quat[i], inertia[i],
                                                                  tflag, ex_space, ey_space, ez_space);
            }
            else if(integration_scheme == SYMPLECTIC)           //symplectic scheme
                symplectic_initial_integrate(dtf2, angmom[i], quat[i], omega[i], torque[i], inertia[i]);
            else if (integration_scheme == PREDICTOR_CORRECTOR) //direct integration, 2nd order predictor-corrector
                predictor_corrector_initial_integrate(dtf, angmom[i], quat[i], omega[i], torque[i], inertia[i]);
            else if(integration_scheme == WOODEM)               //woodem scheme
                woodem_initial_integrate(dtf2, angmom[i], quat[i], omega[i], torque[i], inertia[i]);
            else
                error->one(FLERR,"Internal error");
         }
    }
}

void FixNVENonSpherical::final_integrate_single()
{
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    double *rmass = atom->rmass;
    double **angmom = atom->angmom;
    double **quat = atom->quaternion;
    double **inertia = atom->inertia;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    double rotation_matrix[9];

    double dtf2 = dtf * 2.0;

    for (int i = 0; i < nlocal; i++)
    {
        bool isSingleParticle = true;
        if(fix_multisphere_ && fix_multisphere_->belongs_to(i) >= 0)
          isSingleParticle = false;
        if ((mask[i] & groupbit) && isSingleParticle) // do not integrate if atom belongs to multisphere or concave particle
        {
            const double dtfm = dtf / rmass[i];
            MathExtraLiggghtsNonspherical::quat_to_mat(quat[i], rotation_matrix);

            //update velocity by step t+dt
            v[i][0] += dtfm * f[i][0];
            v[i][1] += dtfm * f[i][1];
            v[i][2] += dtfm * f[i][2];

            if(integration_scheme == RICHARDSON)
            {
                //update angular moment by step t+dt
                bool tflag[3] = {true, true, true};
                double ex_space[3],ey_space[3],ez_space[3];
                richardson_final_integrate(dtf, angmom[i], omega[i], torque[i], quat[i], inertia[i],
                                                  tflag, ex_space, ey_space, ez_space);
            }
            else if(integration_scheme == SYMPLECTIC)
                symplectic_final_integrate(dtf2, angmom[i], quat[i], omega[i], torque[i], inertia[i]);
            else if (integration_scheme == PREDICTOR_CORRECTOR)
                predictor_corrector_final_integrate(dtf, angmom[i], quat[i], omega[i], torque[i], inertia[i]);
            else if(integration_scheme == WOODEM)
            {
                //do nothing
            }
            else
                error->one(FLERR,"Internal error");
        }
    }
}

void FixNVENonSpherical::initial_integrate_clump(int vflag)
{
    double *const *const xcm = multisphere_->get_xcm_ptr();
    double *const *const vcm = multisphere_->get_vcm_ptr();
    double *const *const fcm = multisphere_->get_fcm_ptr();
    double *const *const torquecm = multisphere_->get_torquecm_ptr();
    double *const *const ex_space = multisphere_->get_ex_space_ptr();
    double *const *const ey_space = multisphere_->get_ey_space_ptr();
    double *const *const ez_space = multisphere_->get_ez_space_ptr();
    double *const *const angmom = multisphere_->get_angmom_ptr();
    double *const *const omega = multisphere_->get_omega_ptr();
    double *const *const quat = multisphere_->get_quat_ptr();
    double *const *const inertia = multisphere_->get_inertia_ptr();
    double *const masstotal = multisphere_->get_masstotal_ptr();
    double *const density = multisphere_->get_density_ptr();
    int *const start_step = multisphere_->get_start_step_ptr();
    double *const *const v_integrate = multisphere_->get_v_integrate_ptr();
    bool *const *const fflag = multisphere_->get_fflag_ptr();
    bool *const *const tflag = multisphere_->get_tflag_ptr();

    const int timestep = update->ntimestep;
    const int nbody = multisphere_->n_body();

    const int n_stream = modify->n_fixes_style("insert/stream");
    const bool has_stream = n_stream > 0;

    for (int ibody = 0; ibody < nbody; ibody++)
    {
        if(has_stream && timestep < start_step[ibody])
        {
            vectorCopy3D(v_integrate[ibody],vcm[ibody]);

            // update xcm by full step
            xcm[ibody][0] += dtv * vcm[ibody][0];
            xcm[ibody][1] += dtv * vcm[ibody][1];
            xcm[ibody][2] += dtv * vcm[ibody][2];
            continue;
        }

        // update vcm by 1/2 step

        const double addMassTerm = 1.0+CAddRhoFluid_/density[ibody];
        const double dtfm = dtf / (masstotal[ibody] * addMassTerm );

        if(fflag[ibody][0])
            vcm[ibody][0] += dtfm * fcm[ibody][0];
        if(fflag[ibody][1])
            vcm[ibody][1] += dtfm * fcm[ibody][1];
        if(fflag[ibody][2])
            vcm[ibody][2] += dtfm * fcm[ibody][2];

        // update xcm by full step

        xcm[ibody][0] += dtv * vcm[ibody][0];
        xcm[ibody][1] += dtv * vcm[ibody][1];
        xcm[ibody][2] += dtv * vcm[ibody][2];

        // update angular momentum by 1/2 step
        const double dtt = dtf / addMassTerm;
        richardson_initial_integrate(dtt, angmom[ibody], omega[ibody], torquecm[ibody], quat[ibody], inertia[ibody],
                                          tflag[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody]);
    }

    // virial setup before call to set_xv

    if (vflag)
        v_setup(vflag);
    else
        evflag = 0;

    // set coords/orient and velocity/rotation of atoms in rigid bodies
    // from quarternion and omega

    fix_multisphere_->set_xv();

    //    rev_comm_flag_ = MS_COMM_REV_X_V_OMEGA;
    //    reverse_comm();
    fix_multisphere_->rev_comm(MS_COMM_REV_X_V_OMEGA);

}

void FixNVENonSpherical::final_integrate_clump()
{
    double *const *const vcm = multisphere_->get_vcm_ptr();
    double *const *const fcm = multisphere_->get_fcm_ptr();
    double *const *const torquecm = multisphere_->get_torquecm_ptr();
    double *const *const ex_space = multisphere_->get_ex_space_ptr();
    double *const *const ey_space = multisphere_->get_ey_space_ptr();
    double *const *const ez_space = multisphere_->get_ez_space_ptr();
    double *const *const angmom = multisphere_->get_angmom_ptr();
    double *const *const omega = multisphere_->get_omega_ptr();
    double *const *const quat = multisphere_->get_quat_ptr();
    double *const *const inertia = multisphere_->get_inertia_ptr();
    double *const masstotal = multisphere_->get_masstotal_ptr();
    double *const density = multisphere_->get_density_ptr();
    int *const start_step = multisphere_->get_start_step_ptr();
    bool *const *const fflag = multisphere_->get_fflag_ptr();
    bool *const *const tflag = multisphere_->get_tflag_ptr();

    const int timestep = update->ntimestep;
    const int nbody = multisphere_->n_body();

    // calculate forces and torques on body

    //  calc_force(false); // done now in pre_final integrate

    const int n_stream = modify->n_fixes_style("insert/stream");
    const bool has_stream = n_stream > 0;

    // resume integration
    for (int ibody = 0; ibody < nbody; ibody++)
    {
        /*
        if (!(getMask(ibody) & groupbit))
          continue; */

        if(has_stream && timestep < start_step[ibody])
            continue;

        // update vcm by 1/2 step
        const double addMassTerm = 1.0+CAddRhoFluid_/density[ibody];
        const double dtfm = dtf / ( masstotal[ibody] * addMassTerm );
        if(fflag[ibody][0])
            vcm[ibody][0] += dtfm * fcm[ibody][0];
        if(fflag[ibody][1])
            vcm[ibody][1] += dtfm * fcm[ibody][1];
        if(fflag[ibody][2])
            vcm[ibody][2] += dtfm * fcm[ibody][2];

        // update angular momentum by 1/2 step
        const double dtt = dtf / addMassTerm;
        richardson_final_integrate(dtt, angmom[ibody], omega[ibody], torquecm[ibody], quat[ibody],
            inertia[ibody], tflag[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody]);
    }

    fix_multisphere_->set_v();

//  rev_comm_flag_ = MS_COMM_REV_V_OMEGA;
    fix_multisphere_->rev_comm(MS_COMM_REV_V_OMEGA);

//  fw_comm_flag_ = MS_COMM_FW_V_OMEGA;
//  forward_comm();
    fix_multisphere_->fw_comm(MS_COMM_FW_V_OMEGA);
}

/* ---------------------------------------------------------------------- */

int FixNVENonSpherical::setmask()
{
    int mask = FixNVE::setmask();
    mask |= FixConst::PRE_FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVENonSpherical::setup(int dummy)
{
    if(integrate_convexhull) modify_body_forces_torques();

    FixNVE::setup(dummy);
}

void FixNVENonSpherical::pre_final_integrate()
{
    if(integrate_convexhull) modify_body_forces_torques();
}

void FixNVENonSpherical::modify_body_forces_torques()
{
    int nall = atom->nlocal + atom->nghost;

    // correct forces and torques
    for (int i = 0; i < nall; i++)
    {
        std::list<ContactForceCorrector*>::iterator it;
        for (it = correctors_.begin(); it != correctors_.end(); it++)
            (*it)->modify_body_forces_torques_wall(i);
    }
}

/* ---------------------------------------------------------------------- */

void FixNVENonSpherical::reset_dt()
{
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
    dtq = 0.5 * update->dt;
}

void FixNVENonSpherical::richardson_initial_integrate(double dtf, double *angmom, double *omega, double *torque, double *quat, double *inertia,
                                  bool *tflag, double *ex_space, double *ey_space, double *ez_space)
{
    if(tflag[0])
        angmom[0] += dtf * torque[0];
    if(tflag[1])
        angmom[1] += dtf * torque[1];
    if(tflag[2])
        angmom[2] += dtf * torque[2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::mq_to_omega(angmom,quat,inertia,omega);
    MathExtra::richardson(quat,angmom,omega, inertia, dtq);
    MathExtra::q_to_exyz(quat, ex_space,ey_space,ez_space);
}

void FixNVENonSpherical::richardson_final_integrate(double dt, double *angmom, double *omega, double *torque,  double *quat, double *inertia,
                                  bool *tflag, double *ex_space, double *ey_space, double *ez_space)
{
    if(tflag[0])
        angmom[0] += dt * torque[0];
    if(tflag[1])
        angmom[1] += dt * torque[1];
    if(tflag[2])
        angmom[2] += dt * torque[2];

    MathExtra::mq_to_omega(angmom,quat,inertia,omega);
}

void FixNVENonSpherical::symplectic_initial_integrate(double dtf2, double *angmom, double *quat, double *omega, double *torque, double *inertia)
{
    double tbody[3], rotation_matrix[9];
    double fquat[4], mbody[3], wbody[3], conjqm[4];
    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix);

    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom,mbody);
    MathExtraLiggghtsNonspherical::calc_conjqm(quat,mbody,conjqm);

    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque, tbody);
    MathExtra::quatvec(quat, tbody, fquat);

    conjqm[0] += dtf2 * fquat[0];
    conjqm[1] += dtf2 * fquat[1];
    conjqm[2] += dtf2 * fquat[2];
    conjqm[3] += dtf2 * fquat[3];

    MathExtraLiggghtsNonspherical::no_squish_rotate(3,conjqm,quat,inertia,dtq);
    MathExtraLiggghtsNonspherical::no_squish_rotate(2,conjqm,quat,inertia,dtq);
    MathExtraLiggghtsNonspherical::no_squish_rotate(1,conjqm,quat,inertia,dtv);
    MathExtraLiggghtsNonspherical::no_squish_rotate(2,conjqm,quat,inertia,dtq);
    MathExtraLiggghtsNonspherical::no_squish_rotate(3,conjqm,quat,inertia,dtq);

    MathExtra::qnormalize(quat);
    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix);

    MathExtra::invquatvec(quat,conjqm,mbody);
    mbody[0] *= 0.5;
    mbody[1] *= 0.5;
    mbody[2] *= 0.5;

    wbody[0] = mbody[0] / inertia[0];
    wbody[1] = mbody[1] / inertia[1];
    wbody[2] = mbody[2] / inertia[2];

    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom);
    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega);
}

void FixNVENonSpherical::symplectic_final_integrate(double dtf2, double *angmom, double *quat, double *omega, double *torque, double *inertia)
{
    double tbody[3], rotation_matrix[9];
    double fquat[4], mbody[3], wbody[3];
    double conjqm[4];
    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix);
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque, tbody);
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom, mbody);
    MathExtraLiggghtsNonspherical::calc_conjqm(quat, mbody, conjqm);

    MathExtra::quatvec(quat, tbody, fquat);

    conjqm[0] += dtf2 * fquat[0];
    conjqm[1] += dtf2 * fquat[1];
    conjqm[2] += dtf2 * fquat[2];
    conjqm[3] += dtf2 * fquat[3];

    MathExtra::invquatvec(quat, conjqm, mbody);
    mbody[0] *= 0.5;
    mbody[1] *= 0.5;
    mbody[2] *= 0.5;

    wbody[0] = mbody[0] / inertia[0];
    wbody[1] = mbody[1] / inertia[1];
    wbody[2] = mbody[2] / inertia[2];

    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom);
    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega);
}

void FixNVENonSpherical::predictor_corrector_initial_integrate(double dtf, double *angmom, double *quat, double *omega, double *torque, double *inertia)
{
    double omega_half[3],rotation_matrix[9];
    double mbody[3], wbody[3], tbody[3];
    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix); //calculate rotation matrix from quaternion
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, omega,wbody); //angular velocity in body principal axes
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque, tbody); //torque in body principal axes

    omega_half[0] = wbody[0];
    omega_half[1] = wbody[1];
    omega_half[2] = wbody[2];

    integrate_dynamic_euler(dtf, wbody, tbody, inertia); //calculate angular velocity at step t+dt

    omega_half[0] = 0.5*(wbody[0] + omega_half[0]); //angular velocity at step t+dt/2
    omega_half[1] = 0.5*(wbody[1] + omega_half[1]);
    omega_half[2] = 0.5*(wbody[2] + omega_half[2]);

    integrate_quaternion(quat, wbody); //calculate quaternion at step t+dt

    mbody[0] = inertia[0]*wbody[0]; //calculate angular momentum at step t+dt from angular velocity in body's principal axes
    mbody[1] = inertia[1]*wbody[1];
    mbody[2] = inertia[2]*wbody[2];

    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom); // angular momentum to global frame
    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega); // angular velocity to global frame
}

void FixNVENonSpherical::predictor_corrector_final_integrate(double dtf, double *angmom, double *quat, double *omega, double *torque, double *inertia)
{
    double rotation_matrix[9];
    double mbody[3], wbody[3], tbody[3];
    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix); //calculate rotation matrix from quaternion
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, omega,wbody); //angular velocity in body principal axes
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, torque, tbody); //torque in body principal axes

    integrate_dynamic_euler(dtf, wbody, tbody, inertia); //calculate angular velocity at step t+dt

    mbody[0] = inertia[0]*wbody[0]; //calculate angular momentum at step t+dt from angular velocity in body's principal axes
    mbody[1] = inertia[1]*wbody[1];
    mbody[2] = inertia[2]*wbody[2];

    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, mbody, angmom); // angular momentum to global frame
    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, wbody, omega); // angular velocity to global frame
}

void FixNVENonSpherical::woodem_initial_integrate(double dtf2, double *angmom, double *quat, double *omega, double *torque, double *inertia)
{
    double angmom_half[3], angmom_half_local[3], angmom_next_local[3], rotation_matrix[9];
    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix);
    for(int j = 0; j < 3; j++)
        angmom_half[j] = angmom[j] + torque[j]*dtf;
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom_half, angmom_half_local);
    for(int j = 0; j < 3; j++)
        angmom[j] += torque[j]*dtf2;
    MathExtraLiggghtsNonspherical::transpose_matvec(rotation_matrix, angmom, angmom_next_local);

    double omega_half_local[3], omega_next_local[3], q_der[4], quat_half[4];
    for(int j = 0; j < 3; j++)
    {
        omega_half_local[j] = angmom_half_local[j] / inertia[j];
        omega_next_local[j] = angmom_next_local[j] / inertia[j];
    }

    MathExtra::quatvec(quat, omega_half_local, q_der);
    for(int j = 0; j < 4; j++)
        quat_half[j] = quat[j] + q_der[j]*dtq*0.5;
    MathExtra::qnormalize(quat_half);

    MathExtra::quatvec(quat_half, omega_next_local, q_der);
    for(int j = 0; j < 4; j++)
        quat[j] += q_der[j]*dtq;
    MathExtra::qnormalize(quat);

    MathExtraLiggghtsNonspherical::quat_to_mat(quat, rotation_matrix);
    MathExtraLiggghtsNonspherical::matvec(rotation_matrix, omega_next_local, omega);
}
