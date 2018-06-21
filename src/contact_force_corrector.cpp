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
#include "contact_force_corrector.h"
#include "modify.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ContactForceCorrector::ContactForceCorrector(LAMMPS *lmp,class FixContactPropertyAtom *fpca,
                                            FixMultisphere const*fix_ms,class Multisphere *msp,
                                            bool wall, const char *wall_id) :
    Pointers(lmp),
    helper_(NULL),
    fix_contact_(fpca),
    fix_ms_(fix_ms),
    multisphere_(msp),
    max_n_ibody_contact_ 
    (
        wall ?
            NULL
        :
        (
            
            multisphere_->prop().addElementProperty< ScalarContainer<int> >
            ("max_n_ibody_contact","comm_exchange_borders","frame_invariant","restart_no",1,multisphere_->n_body())
        )
    ),
    n_ibody_contact_ 
    (
        fix_ms ?
        (
            wall ?
            (
                
                multisphere_->prop().addElementProperty< ScalarContainer<int> >
                (my_strcat("n_ibody_contact_",wall_id),"comm_exchange_borders","frame_invariant","restart_no",1,multisphere_->n_body())
            )
            :
            (
                
                multisphere_->prop().addElementProperty< ScalarContainer<int> >
                ("n_ibody_contact","comm_exchange_borders","frame_invariant","restart_no",1,multisphere_->n_body())
            )
        )
        :
            NULL
    ),
    is_wall_(wall),
    nbody_max_(0),
    jbodytags_(0),
    jbodytags_page_(0),
    ncontacts_(0),
    ncontacts_page_(0),
    pgsize_(0),
    oneatom_(0)
{ }

/* ---------------------------------------------------------------------- */

const char * ContactForceCorrector::my_strcat(const char* s1,const char* s2)
{
    helper_ = new char[200];
    sprintf(helper_,"%s%s",s1,s2);
    return helper_;
}

/* ---------------------------------------------------------------------- */

ContactForceCorrector::~ContactForceCorrector()
{
    if(!is_wall_ && fix_ms_)
        multisphere_->prop().removeElementProperty("max_n_ibody_contact");

    if (fix_ms_)
    {
        char container_id[200];
        n_ibody_contact_->id(container_id);
        multisphere_->prop().removeElementProperty(container_id);

        // delete locally stored arrays

        memory->sfree(jbodytags_);
        memory->sfree(ncontacts_);
        if(jbodytags_page_) delete [] jbodytags_page_;
        if(ncontacts_page_) delete [] ncontacts_page_;
    }
}

/* ----------------------------------------------------------------------
  create pages if first time or if neighbor pgsize/oneatom has changed
------------------------------------------------------------------------- */

void ContactForceCorrector::allocate_pages()
{
  
  if(is_wall_) return;

  int create = 0;
  if (jbodytags_page_ == NULL) create = 1;
  if (pgsize_ != neighbor->pgsize) create = 1;
  if (oneatom_ != neighbor->oneatom) create = 1;

  if (create) {
    delete [] jbodytags_page_;
    delete [] ncontacts_page_;

    pgsize_ = neighbor->pgsize;
    oneatom_ = neighbor->oneatom;
    int nmypage = comm->nthreads;
    jbodytags_page_ = new MyPage<int>[nmypage];
    ncontacts_page_ = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++) {
      jbodytags_page_[i].init(oneatom_,pgsize_);
      ncontacts_page_[i].init(oneatom_,pgsize_);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ContactForceCorrector::calculate_pair()
{
    if (!fix_ms_)
        return;

    allocate_pages();

    allocate_pair();

    // fill data structure
    calculate_contacts_pair();
}

/* ---------------------------------------------------------------------- */

void ContactForceCorrector::calculate_wall()
{
    if (!fix_ms_)
        return;
    
    // fill data structure
    calculate_contacts_wall();
}

/* ---------------------------------------------------------------------- */

void ContactForceCorrector::allocate_pair()
{
    int nbody = multisphere_->n_body();
    int nall = atom->nlocal + atom->nghost;
    int n_atom_contacts,i_body_tag,ibody,j_atom_tag,j,j_body_tag;
    std::vector<int>contact_body_ids;

    jbodytags_page_->reset();
    ncontacts_page_->reset();

    if(nbody_max_ < nbody)
    {
        nbody_max_ = nbody + 10000;
        jbodytags_ = (int **) memory->srealloc(jbodytags_,
                 nbody_max_*sizeof(int *),"ContactForceCorrector:jbodytags");
        ncontacts_ = (int **) memory->srealloc(ncontacts_,
                 nbody_max_*sizeof(int *),"ContactForceCorrector:ncontacts");
    }

    max_n_ibody_contact_->setAll(0);

    for(int i = 0; i < nall; i++)
    {
        i_body_tag = fix_ms_->belongs_to(i);

        if(i_body_tag < 0)
            continue;

        ibody = multisphere_->map(i_body_tag);

        if(ibody < 0)
            continue;

        contact_body_ids.clear();
        n_atom_contacts = fix_contact_->n_partner(i);

        for(int iatomcontact = 0; iatomcontact < n_atom_contacts; iatomcontact++)
        {
            
            j_atom_tag = fix_contact_->partner(i,iatomcontact);

            j = atom->map(j_atom_tag);

            if(j < 0)
            {
                
                continue;
            }

            j_body_tag = fix_ms_->belongs_to(j);

            if(j_body_tag > -1)
            {
                
                if(std::find(contact_body_ids.begin(), contact_body_ids.end(), j_body_tag) == contact_body_ids.end())
                {
                    contact_body_ids.push_back(j_body_tag);
                    
                }
            }
        }

        if(contact_body_ids.size() > 0)
            (*max_n_ibody_contact_)(ibody) += contact_body_ids.size();
    }

    for(int ibody = 0; ibody < nbody; ibody++)
    {
        
        const int max_contact = (*max_n_ibody_contact_)(ibody);
        jbodytags_[ibody] = jbodytags_page_->get(max_contact);
        ncontacts_[ibody] = ncontacts_page_->get(max_contact);

        vectorInitializeN(jbodytags_[ibody],max_contact,-1);
        vectorInitializeN(ncontacts_[ibody],max_contact,0);
    }
}

/* ---------------------------------------------------------------------- */

void ContactForceCorrector::calculate_contacts_pair()
{
    int nall = atom->nlocal + atom->nghost;
    int n_atom_contacts,ibody,i_body_tag,j_atom_tag,j,j_body_tag;

    n_ibody_contact_->setAll(0);

    for(int i = 0; i < nall; i++)
    {
        
        n_atom_contacts = fix_contact_->n_partner(i);

        i_body_tag = fix_ms_->belongs_to(i);

        if(i_body_tag < 0)
            continue;

        ibody = multisphere_->map(i_body_tag);

        if(ibody < 0)
            continue;

        for(int iatomcontact = 0; iatomcontact < n_atom_contacts; iatomcontact++)
        {
            
            j_atom_tag = fix_contact_->partner(i,iatomcontact);

            j = atom->map(j_atom_tag);

            j_body_tag = fix_ms_->belongs_to(j);

            if(j_body_tag > -1)
            {
                
                int ibody_contact = 0;
                bool new_in_list = true;

                while(ibody_contact < (*n_ibody_contact_)(ibody))
                {
                    if(j_body_tag == jbodytags_[ibody][ibody_contact])
                    {
                        new_in_list = false;
                        break;
                    }

                    ibody_contact++;
                }

                if(new_in_list)
                {
                    jbodytags_[ibody][ibody_contact] = j_body_tag;
                    
                    (*n_ibody_contact_)(ibody)++;
                }

                ncontacts_[ibody][ibody_contact]++;
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void ContactForceCorrector::calculate_contacts_wall()
{
    int nall = atom->nlocal + atom->nghost;
    int n_wall_contacts,ibody,i_body_tag;

    n_ibody_contact_->setAll(0);

    for(int i = 0; i < nall; i++)
    {
        
        n_wall_contacts = fix_contact_->n_partner(i);

        i_body_tag = fix_ms_->belongs_to(i);
        if(i_body_tag < 0)
            continue;

        ibody = multisphere_->map(i_body_tag);

        if(ibody > -1)
            (*n_ibody_contact_)(ibody) += n_wall_contacts;
    }

}
