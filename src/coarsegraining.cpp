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

#include "lmptype.h"
#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "modify.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "error.h"
#include "coarsegraining.h"

#define VERBOSE_COARSEGRAINING false
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Coarsegraining::Coarsegraining(LAMMPS *lmp) : Pointers(lmp)
{
}

/* ---------------------------------------------------------------------- */

Coarsegraining::~Coarsegraining()
{
}

/* ----------------------------------------------------------------------
   called as coarsegraining command in input script
------------------------------------------------------------------------- */

void Coarsegraining::command(int narg, char **arg)
{
  // error checks
  // ensure is first in input script, so all commands using length scales
  // can use it

  if(neighbor->skin != 0.3 || domain->box_exist)
    error->all(FLERR,"Illegal coarsegraining command, must use this command "
                     "first in input script");

  if(modify->nfix > 0)
    error->all(FLERR,"Illegal coarsegraining command, must use this command "
                     "before any 'fix' command");

  // parse arguments

  if (narg < 1) error->all(FLERR,"Illegal coarsegraining command");

  int iarg = 0;
  double cg = force->numeric(FLERR,arg[iarg++]);
  if(cg < 1.)
    error->all(FLERR,"Illegal coarsegraining command, cg > 1 expected");

  // set coarsegraining in force class
  bool typeSpecificCG(false);
  typeSpecificCG = force->setCG(cg);

#if VERBOSE_COARSEGRAINING
  printf("\n\nCoarsegraining:typeSpecificCG: %d, value: %g. \n", typeSpecificCG,cg);
#endif

  while(iarg < narg)
  {
      if(strcmp(arg[iarg],"model_check") == 0) {
          if (narg < iarg+2)
            error->all(FLERR,"Illegal coarsegraining command, not enough arguments");
          if(strcmp(arg[iarg+1],"error") == 0) {
            force->error_coarsegraining_ = true;
            force->warn_coarsegraining_ = false;
          } else if(strcmp(arg[iarg+1],"warn") == 0) {
            force->error_coarsegraining_ = false;
            force->warn_coarsegraining_ = true;
          } else if(strcmp(arg[iarg+1],"off") == 0) {
            force->error_coarsegraining_ = false;
            force->warn_coarsegraining_ = false;
          }
          else
            error->all(FLERR,"Illegal coarsegraining command, expecting 'error' or 'warn' after 'model_check'");
          iarg += 2;
      }
      else
      {
        cg = force->numeric(FLERR,arg[iarg++]);
        if(cg < 1.)
            error->all(FLERR,"Illegal coarsegraining command, cg > 1 expected, or illegal coarsegraining command, unknown keyword");

        typeSpecificCG = force->setCG(cg);
        #if VERBOSE_COARSEGRAINING
           printf("\n\nCoarsegraining:typeSpecificCG: %d, value: %g. \n", typeSpecificCG,cg);
        #endif
  }
  }
  #if VERBOSE_COARSEGRAINING
    force->reportCG();
  #endif
}
