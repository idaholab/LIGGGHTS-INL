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
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include <mpi.h>
#include "compute_sauter_diameter.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "group.h"
#include "error.h"

#include "math_extra_liggghts_nonspherical.h"

/* ---------------------------------------------------------------------- */

ComputeSauterDiameter::ComputeSauterDiameter(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg != iarg) error->compute_error(FLERR, this, "Not enough arguments for compute SauterDiameter");

  scalar_flag = 1;
  extscalar = 1;

  // error check

  if (!atom->superquadric_flag)
    error->all(FLERR,"Compute SauterDiameter requires atom style superquadric");

}

/* ---------------------------------------------------------------------- */

void ComputeSauterDiameter::init()
{
}

/* ---------------------------------------------------------------------- */

double ComputeSauterDiameter::compute_scalar()
{

  double *area = atom->area;
  double *vol = atom->volume;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double diameter = 0.0;
  int num = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      diameter += 6.0*vol[i]/area[i];
      num++;
    }
  double diameter_sum;
  int num_sum;
  MPI_Allreduce(&diameter,&diameter_sum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&num,&num_sum,1,MPI_INT,MPI_SUM,world);
  if(num_sum > 0)
    scalar = diameter_sum/static_cast<double>(num_sum);
  else
    scalar = 0.0;
  return scalar;
}
#endif
