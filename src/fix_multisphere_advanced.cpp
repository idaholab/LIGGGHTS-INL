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

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_multisphere_advanced.h"
#include "fix_contact_property_atom_wall.h"
#include "modify.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "pair_gran.h"
#include "fix_wall_gran.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMultisphereAdvanced::FixMultisphereAdvanced(LAMMPS *lmp, int narg, char **arg) :
  FixMultisphere(lmp, narg, arg)
{
    
    if(0 == strcmp(style,"concave/advanced"))
    {
        concave_ = true;
        int strln = strlen("multisphere/advanced") + 1;
        delete []style;
        style = new char[strln];
        strcpy(style,"multisphere/advanced");
    }

    do_modify_body_forces_torques_ = true;

    if(!force->pair_match("gran", 0))
          error->fix_error(FLERR,this,"Please use a granular pair style before using this fix");
    static_cast<PairGran*>(force->pair_match("gran", 0))->do_store_contact_forces();

    if(accepts_restart_data_from_style)
        delete []accepts_restart_data_from_style;
    accepts_restart_data_from_style = new char[12];
    sprintf(accepts_restart_data_from_style,"multisphere");
}

/* ---------------------------------------------------------------------- */

FixMultisphereAdvanced::~FixMultisphereAdvanced()
{
    if(accepts_restart_data_from_style)
    {
        delete []accepts_restart_data_from_style;
        accepts_restart_data_from_style = 0;
    }
}

/* ---------------------------------------------------------------------- */

int FixMultisphereAdvanced::setmask()
{
    int mask = FixMultisphere::setmask();
    mask |= MIN_POST_FORCE;
    mask |= POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

double FixMultisphereAdvanced::extend_cut_ghost()
{
    
    return 2.*max_r_bound();
}

/* ---------------------------------------------------------------------- */

void FixMultisphereAdvanced::init()
{
    FixMultisphere::init();

    if(0 == atom->map_style)
      error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");

    while(correctors_.size() > 0)
    {
        delete correctors_[0];
        correctors_.erase(correctors_.begin());
    }

    FixContactPropertyAtom *fix_contact_pair = static_cast<FixContactPropertyAtom*>(modify->find_fix_id("contactforces"));
    if(!fix_contact_pair)
        error->fix_error(FLERR,this,"Internal error: need fix contactforces");

    correctors_.push_back(new ContactForceCorrector(lmp,fix_contact_pair,this,&multisphere_,false,""));

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
                correctors_.push_back(new ContactForceCorrector(lmp,fix_contact_wall,this,&multisphere_,true,fixid));
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
            correctors_.push_back(new ContactForceCorrector(lmp,fix_contact_wall,this,&multisphere_,true,fixid));
        }
    }

    if(modify->my_index(this) <  modify->index_last_fix_of_style("wall/gran"))
        error->fix_error(FLERR,this,"needs to be after all fixes of type wall/gran");

    if(modify->n_fixes_style("setforce")>0 || modify->n_fixes_style("freeze")>0)
        error->fix_error(FLERR,this,"not compatible to setforce or freeze fixes");
}

/* ---------------------------------------------------------------------- */

void FixMultisphereAdvanced::setup(int dummy)
{
    
    post_force(dummy);

    FixMultisphere::setup(dummy);
}

/* ---------------------------------------------------------------------- */

void FixMultisphereAdvanced::post_force(int dummy)
{
    UNUSED(dummy);

    correctors_[0]->calculate_pair();

    int ncorr = correctors_.size();
    for(int icorr = 1; icorr < ncorr; icorr++)
        correctors_[icorr]->calculate_wall();
}

/* ---------------------------------------------------------------------- */

void FixMultisphereAdvanced::modify_body_forces_torques()
{
    int nall = atom->nlocal + atom->nghost;

    // correct forces and torques
    for (int i = 0; i < nall; i++)
    {
        
        if(body_[i] < 0) continue;

        correctors_[0]->modify_body_forces_torques_pair(i);

        int ncorr = correctors_.size();
        for(int icorr = 1; icorr < ncorr; icorr++)
            correctors_[icorr]->modify_body_forces_torques_wall(i);
    }
}
