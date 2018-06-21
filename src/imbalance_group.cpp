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

    Arno Mayrhofer (DCS Computing GmbH, Linz)

    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2018-     DCS Computing GmbH, Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include "imbalance_group.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* -------------------------------------------------------------------- */

ImbalanceGroup::ImbalanceGroup(LAMMPS *lmp) : Imbalance(lmp), id(0), factor(0)
{}

/* -------------------------------------------------------------------- */

ImbalanceGroup::~ImbalanceGroup()
{
  delete [] id;
  delete [] factor;
}

/* -------------------------------------------------------------------- */

int ImbalanceGroup::options(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal balance weight command");

  num = force->inumeric(FLERR,arg[0]);
  if (num < 1) error->all(FLERR,"Illegal balance weight command");
  if (2*num+1 > narg) error->all(FLERR,"Illegal balance weight command");

  id = new int[num];
  factor = new double[num];
  for (int i = 0; i < num; ++i) {
    id[i] = group->find(arg[2*i+1]);
    if (id[i] < 0)
      error->all(FLERR,"Unknown group in balance weight command");
    factor[i] = force->numeric(FLERR,arg[2*i+2]);
    if (factor[i] <= 0.0) error->all(FLERR,"Illegal balance weight command");
  }
  return 2*num+1;
}

/* -------------------------------------------------------------------- */

void ImbalanceGroup::compute(double *weight)
{
  const int * const mask = atom->mask;
  const int * const bitmask = group->bitmask;
  const int nlocal = atom->nlocal;

  if (num == 0) return;

  for (int i = 0; i < nlocal; ++i) {
    const int imask = mask[i];
    for (int j = 0; j < num; ++j) {
      if (imask & bitmask[id[j]])
        weight[i] *= factor[j];
    }
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceGroup::info(FILE *fp)
{
  if (num > 0) {
    const char * const * const names = group->names;

    fprintf(fp,"  group weights:");
    for (int i = 0; i < num; ++i)
      fprintf(fp," %s=%g",names[id[i]],factor[i]);
    fputs("\n",fp);
  }
}
