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

#ifndef LMP_MULTISPHERE_PARALLEL_H
#define LMP_MULTISPHERE_PARALLEL_H

#define LMP_MULTISPHERE_PARALLEL_FLAG_

#include "multisphere.h"
#include "parallel_base.h"
#include "comm.h"

namespace LAMMPS_NS {

  class MultisphereParallel : public Multisphere {

    public:

      MultisphereParallel(LAMMPS *lmp);
      ~MultisphereParallel()
      { }

      void exchange()
      { mycomm_->exchange(); }

      void restart_serial(double *list);
      void write_restart_parallel(FILE *fp);
      void restart_parallel(double *list);

      // ParallelBase interface
      void clear_exchange()
      { mycomm_->copy_exchange_info_from_main_comm(); }
      int pack_exchange(const int dim);
      void unpack_exchange(const int nrecv, const int dim);
      void post_exchange();

      int pack_irregular(const int i, double * const buf);
      int unpack_irregular(double * const buf);

      int size_restart() const;
      int pack_restart(double *const buf) const;
      int unpack_restart(double *const buf);
      void finalize_restart();

    private:

      int pack_exchange_rigid(int i, double *buf);
      int unpack_exchange_rigid(double *buf);

      int nbody_all_restart;

      Comm *mycomm_;

  };

  // *************************************
  #include "multisphere_parallel_I.h"
  // *************************************
}

#endif
