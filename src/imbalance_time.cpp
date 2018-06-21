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

#include <mpi.h>
#include "imbalance_time.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* -------------------------------------------------------------------- */

ImbalanceTime::ImbalanceTime(LAMMPS *lmp) : Imbalance(lmp) {}

/* -------------------------------------------------------------------- */

int ImbalanceTime::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  factor = force->numeric(FLERR,arg[0]);
  if (factor <= 0.0) error->all(FLERR,"Illegal balance weight command");
  return 1;
}

/* ----------------------------------------------------------------------
   reset last and timers if necessary
------------------------------------------------------------------------- */

void ImbalanceTime::init(int flag)
{
  last = 0.0;

  // flag = 1 if called from FixBalance at start of run
  //   init Timer, so accumulated time not carried over from previous run
  // should NOT init Timer if called from Balance, it uses time from last run

  if (flag) timer->init(true);
}

/* -------------------------------------------------------------------- */

void ImbalanceTime::compute(double *weight)
{
  // cost = CPU time for relevant timers since last invocation
  // localwt = weight assigned to each owned atom
  // just return if no time yet tallied

  double cost = -last;
  cost += timer->cpu_elapsed(TIME_PAIR);
  cost += timer->cpu_elapsed(TIME_NEIGHBOR);
  cost += timer->cpu_elapsed(TIME_BOND);
  cost += timer->cpu_elapsed(TIME_KSPACE);

  double maxcost;
  MPI_Allreduce(&cost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);
  if (maxcost <= 0.0) return;

  int nlocal = atom->nlocal;
  double localwt = 0.0;
  if (nlocal) localwt = cost/nlocal;

  if (nlocal && localwt <= 0.0) error->one(FLERR,"Balance weight <= 0.0");

  // apply factor if specified != 1.0
  // wtlo,wthi = lo/hi values excluding 0.0 due to no atoms on this proc
  // lo value does not change
  // newhi = new hi value to give hi/lo ratio factor times larger/smaller
  // expand/contract all localwt values from lo->hi to lo->newhi

  if (factor != 1.0) {
    double wtlo,wthi;
    if (localwt == 0.0) localwt = BIG;
    MPI_Allreduce(&localwt,&wtlo,1,MPI_DOUBLE,MPI_MIN,world);
    if (localwt == BIG) localwt = 0.0;
    MPI_Allreduce(&localwt,&wthi,1,MPI_DOUBLE,MPI_MAX,world);
    if (wtlo == wthi) return;

    double newhi = wthi*factor;
    localwt = wtlo + ((localwt-wtlo)/(wthi-wtlo)) * (newhi-wtlo);
  }

  for (int i = 0; i < nlocal; i++) weight[i] *= localwt;

  // record time up to this point

  last += cost;
}

/* -------------------------------------------------------------------- */

void ImbalanceTime::info(FILE *fp)
{
  fprintf(fp,"  time weight factor: %g\n",factor);
}
