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

#ifdef COMPUTE_CLASS

ComputeStyle(crosssection,ComputeCrosssection)

#else

#ifndef LMP_COMPUTE_CROSSSECTION_H
#define LMP_COMPUTE_CROSSSECTION_H

#include "compute_contact_atom.h"

namespace LAMMPS_NS {

class ComputeCrosssection : public ComputeContactAtom
{
 public:

  ComputeCrosssection(class LAMMPS *, int &iarg, int, char **);
  ~ComputeCrosssection();

  void init();

  void compute_vector();
  void compute_peratom();

  void write();

 private:

  void setup_cuts();
  void compute_convex_hull();
  inline int mod(double coo);
  double calc_ang();

  class ModifiedAndrew **mod_andrew_;
  int dim_, n_cuts_;
  double min_,max_;
  double cut_thickness_half_;
  double cut_dist_;

  bool angle_;

  FILE *file_;
  bool fileflag_;
  int iregion_;
  char *idregion_;
};

}

#endif
#endif
