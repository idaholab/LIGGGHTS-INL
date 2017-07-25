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

#ifdef FIX_CLASS

FixStyle(insert/fragments,FixInsertFragments)

#else

#ifndef LMP_FIX_INSERT_FRAGMENTS_H
#define LMP_FIX_INSERT_FRAGMENTS_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertFragments : public FixInsert {
 public:

  FixInsertFragments(class LAMMPS *, int, char **);
  ~FixInsertFragments();

  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();

  void init_defaults();
  int calc_ninsert_this();

  BoundingBox getBoundingBox();

  void pre_neighbor();
  void setup_pre_neighbor();
  virtual void end_of_step();

  double insertion_fraction()
  { return 0.; }

  double extend_cut_ghost()
  { return 0.; }

 protected:

  // inherited functions
  virtual void calc_insertion_properties();
  virtual bool pre_insert();
  int distribute_ninsert_this(int);
  int is_nearby(int);
  void x_v_omega(int ninsert_this,int &ninserted_this,
                 int &ninserted_spheres_this,
                 double &mass_inserted_this);

  // functions declared in this class
  inline void generate_random(double *pos, double rad_broken,
                              double rad_insert);
  void print_stats_replace_during();

  // template holding data of the fragments
  class FixTemplateMultiplespheres *fix_fragments_;

  // flag if perform saftety check on bound sphere
  bool check_r_bound_;

  // number of fragments per particle
  int n_fragments_;

  // replacement flag
  class FixPropertyAtom *fix_replace_;

  // particle must be of this type to be replaced
  char *idregion_;
  class Region *ins_region_;

  // multisphere support
  class FixMultisphere *fix_ms_;
  class Multisphere *ms_;

  int type_replace_;

  // stats for replacements
  int n_replace_,n_replace_this_,n_replace_this_local_;
  double mass_replace_,mass_replace_this_,mass_replace_this_local_;

  // data for particles to be inserted
  double **replacedata_;
};

}

#endif
#endif
