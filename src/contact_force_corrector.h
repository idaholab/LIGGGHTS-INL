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

#ifndef LMP_CONTACT_FORCE_CORRECTOR_H
#define LMP_CONTACT_FORCE_CORRECTOR_H

#include "pointers.h"
#include "fix_multisphere.h"
#include "container.h"
#include "my_page.h"
#include "atom.h"
#include "fix_contact_property_atom.h"
#include "fix_mesh_surface.h"

namespace LAMMPS_NS {

class ContactForceCorrector : protected Pointers
{
    public:

      ContactForceCorrector(LAMMPS *lmp,class FixContactPropertyAtom *fpca,
                            FixMultisphere const*fix_ms,class Multisphere *msp,
                            bool wall,const char *wall_id);
      ~ContactForceCorrector();

      void allocate_pages();

      void calculate_pair();
      void calculate_wall();

      inline void modify_body_forces_torques_pair(int i);
      inline void modify_body_forces_torques_wall(int i);

    private:

      void allocate_pair();
      void calculate_contacts_pair();

      void calculate_contacts_wall();

      const char *my_strcat(const char* s1,const char* s2);

      char *helper_;

      class FixContactPropertyAtom *fix_contact_;

      FixMultisphere const *fix_ms_;
      class Multisphere *multisphere_;

      ScalarContainer<int> *max_n_ibody_contact_;

      ScalarContainer<int> *n_ibody_contact_;

      bool is_wall_;

      int nbody_max_;

      int **jbodytags_;
      MyPage<int> *jbodytags_page_;

      int **ncontacts_;
      MyPage<int> *ncontacts_page_;

      int pgsize_,oneatom_;

};

#include "contact_force_corrector_I.h"

}

#endif
