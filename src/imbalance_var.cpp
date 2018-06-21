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
#include "imbalance_var.h"
#include "atom.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

// DEBUG
#include "update.h"

using namespace LAMMPS_NS;

/* -------------------------------------------------------------------- */

ImbalanceVar::ImbalanceVar(LAMMPS *lmp) : Imbalance(lmp), name(0) {}

/* -------------------------------------------------------------------- */

ImbalanceVar::~ImbalanceVar()
{
  delete [] name;
}

/* -------------------------------------------------------------------- */

int ImbalanceVar::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal balance weight command");

  int len = strlen(arg[0]) + 1;
  name = new char[len];
  memcpy(name,arg[0],len);
  init(0);

  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceVar::init(int flag)
{
  id = input->variable->find(name);
  if (id < 0) {
    error->all(FLERR,"Variable name for balance weight does not exist");
  } else {
    if (input->variable->atomstyle(id) == 0)
      error->all(FLERR,"Variable for balance weight has invalid style");
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceVar::compute(double *weight)
{
  const int all = group->find("all");
  if (all < 0) return;

  double *values;
  const int nlocal = atom->nlocal;
  memory->create(values,nlocal,"imbalance:values");

  input->variable->compute_atom(id,all,values,1,0);

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (values[i] <= 0.0) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->one(FLERR,"Balance weight <= 0.0");

  for (int i = 0; i < nlocal; i++) weight[i] *= values[i];

  memory->destroy(values);
}

/* -------------------------------------------------------------------- */

void ImbalanceVar::info(FILE *fp)
{
  fprintf(fp,"  weight variable: %s\n",name);
}
