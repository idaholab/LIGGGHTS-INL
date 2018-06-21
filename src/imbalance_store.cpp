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

#include <string.h>
#include "imbalance_store.h"
#include "atom.h"
#include "modify.h"
#include "fix_property_atom.h"
#include "input.h"
#include "error.h"

using namespace LAMMPS_NS;

/* -------------------------------------------------------------------- */

ImbalanceStore::ImbalanceStore(LAMMPS *lmp) : Imbalance(lmp), name(0) {}

/* -------------------------------------------------------------------- */

ImbalanceStore::~ImbalanceStore()
{
  delete [] name;
}

/* -------------------------------------------------------------------- */

int ImbalanceStore::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal balance weight command");

  int len = strlen(arg[0]) + 1;
  name = new char[len];
  memcpy(name,arg[0],len);

  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::compute(double *weight)
{
  FixPropertyAtom *fix_weight = static_cast<FixPropertyAtom*>(modify->find_fix_property(name, "property/atom", "scalar", 0, 0, "Imbalance/store"));

  double * const prop = fix_weight->vector_atom;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i)
    prop[i] = weight[i];
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::info(FILE *fp)
{
  fprintf(fp,"  storing weight in atom property d_%s\n",name);
}
