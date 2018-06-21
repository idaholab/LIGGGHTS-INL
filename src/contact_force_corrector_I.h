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

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_CONTACT_FORCE_CORRECTOR_I_H
#define LMP_CONTACT_FORCE_CORRECTOR_I_H

/* ---------------------------------------------------------------------- */

inline void ContactForceCorrector::modify_body_forces_torques_pair(int i)
{
    if (!fix_ms_)
        error->all(FLERR, "Contact force corrector only implemented for multispheres");

    double correction_factor, f_torque_contact[6];
    int n_atom_contacts = fix_contact_->n_partner(i);
    int j,j_tag,j_body_tag, ibody, ibody_contact;
    int i_body_tag = fix_ms_->belongs_to(i);

    if(i_body_tag < 0)
        return;

    ibody = multisphere_->map(i_body_tag);

    if(ibody < 0)
        return;

    for(int iatomcontact = 0; iatomcontact < n_atom_contacts; iatomcontact++)
    {
        
        j_tag = fix_contact_->partner(i,iatomcontact);
        j = atom->map(j_tag);
        j_body_tag = fix_ms_->belongs_to(j);

        if(j_body_tag > -1)
        {
            ibody_contact = 0;

            const int n_contact = (*n_ibody_contact_)(ibody);
            while(ibody_contact < n_contact)
            {
                if(j_body_tag == jbodytags_[ibody][ibody_contact])
                        break;
                ibody_contact++;
            }
            
            if(ibody_contact == n_contact)
            {
                
                error->one(FLERR,"internal error");
            }
            if(0 == ncontacts_[ibody][ibody_contact])
                error->one(FLERR,"internal error");

            correction_factor = 1.-1./static_cast<double>(ncontacts_[ibody][ibody_contact]);

            fix_contact_->contacthistory(i,iatomcontact,f_torque_contact);

            atom->f[i][0] -= correction_factor * f_torque_contact[0];
            atom->f[i][1] -= correction_factor * f_torque_contact[1];
            atom->f[i][2] -= correction_factor * f_torque_contact[2];
            atom->torque[i][0] -= correction_factor * f_torque_contact[3];
            atom->torque[i][1] -= correction_factor * f_torque_contact[4];
            atom->torque[i][2] -= correction_factor * f_torque_contact[5];
            
        }
    }
}

/* ---------------------------------------------------------------------- */

inline void ContactForceCorrector::modify_body_forces_torques_wall(int i)
{
    double correction_factor, f_torque_contact[6];
    int n_wall_contacts = fix_contact_->n_partner(i);

    if (n_wall_contacts == 0)
        return;

    if (fix_ms_)
    {
        int ibody, i_body_tag = fix_ms_->belongs_to(i);

        if(i_body_tag < 0)
            return;

        ibody = multisphere_->map(i_body_tag);

        if(ibody < 0)
            return;

        const int n_contact = (*n_ibody_contact_)(ibody);
        correction_factor = 1.-1./static_cast<double>(n_contact);
    }
    else
        correction_factor = 1.-1./static_cast<double>(n_wall_contacts);

    for(int iwallcontact = 0; iwallcontact < n_wall_contacts; iwallcontact++)
    {
        fix_contact_->contacthistory(i,iwallcontact,f_torque_contact);

        atom->f[i][0] -= correction_factor * f_torque_contact[0];
        atom->f[i][1] -= correction_factor * f_torque_contact[1];
        atom->f[i][2] -= correction_factor * f_torque_contact[2];
        atom->torque[i][0] -= correction_factor * f_torque_contact[3];
        atom->torque[i][1] -= correction_factor * f_torque_contact[4];
        atom->torque[i][2] -= correction_factor * f_torque_contact[5];
        
        FixMeshSurface* mesh = fix_contact_->getMesh();
        if (mesh) {
            
            // first simplified version
            
            double delta[3] = {0.0, 0.0, 0.0}, vrel[3] = {0.0, 0.0, 0.0};
            double f_corr[3] = {f_torque_contact[0], f_torque_contact[1], f_torque_contact[2]};
            vectorScalarMult3D(f_corr,-correction_factor,f_corr);
            mesh->add_particle_contribution(i,f_corr,delta,fix_contact_->partner(i,iwallcontact),vrel,fix_contact_->contacthistory(i,iwallcontact));
        }
    }
    
}

#endif
