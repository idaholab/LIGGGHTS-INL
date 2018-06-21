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

    Christoph Kloss (JKU Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

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
#include <stdlib.h>
#include "fix_balance.h"
#include "balance.h"
#include "update.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "irregular.h"
#include "force.h"
#include "kspace.h"
#include "modify.h"
#include "rcb.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SHIFT,BISECTION};

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), balance(NULL), irregular(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix balance command");

  box_change_domain = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
  global_freq = 1;

  execute_flag = false;

  // parse required arguments

  int dimension = domain->dimension;

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix balance command");

  bool old_style = false;
  if (strcmp(arg[5],"shift") == 0) lbstyle = SHIFT;
  else if (strcmp(arg[5],"rcb") == 0) lbstyle = BISECTION;
  else
  {
    lbstyle = SHIFT;
    if (comm->me == 0)
      error->warning(FLERR, "No keyword 'shift' or 'rcb' found when initializing fix balance. Assuming old style fix balance input. Please consult the documentation and update your input script.");
    old_style = true;
  }

  if (!old_style)
    thresh = force->numeric(FLERR,arg[4]);

  int iarg = 5;
  if (lbstyle == SHIFT && !old_style) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal fix balance command");
    if (strlen(arg[iarg+1]) > 3)
      error->all(FLERR,"Illegal fix balance command");
    strcpy(bstr,arg[iarg+1]);
    nitermax = force->inumeric(FLERR,arg[iarg+2]);
    if (nitermax <= 0) error->all(FLERR,"Illegal fix balance command");
    stopthresh = force->numeric(FLERR,arg[iarg+3]);
    if (stopthresh < 1.0) error->all(FLERR,"Illegal fix balance command");
    iarg += 4;
  } else if (lbstyle == BISECTION) {
    iarg++;
  }

  if (old_style)
  {
    if (narg < 7)
      error->all(FLERR,"Illegal fix balance command");
    if (strlen(arg[4]) > 3) error->all(FLERR,"Illegal fix balance command");
    strcpy(bstr,arg[4]);
    nitermax = force->inumeric(FLERR,arg[5]);
    stopthresh = force->numeric(FLERR,arg[6]);

    thresh = stopthresh;
    iarg = 7;
  }

  // error checks

  if (lbstyle == SHIFT) {
    int blen = strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Fix balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Fix balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Fix balance shift string is invalid");
    }
    if (nitermax <= 0 || stopthresh < 1.0)
      error->all(FLERR,"Illegal fix balance command");
  }

  if (lbstyle == BISECTION && comm->style == 0)
    error->all(FLERR,"Fix balance rcb cannot be used with comm_style brick");

  // create instance of Balance class
  // if SHIFT, initialize it with params
  // process remaining optional args via Balance

  balance = new Balance(lmp);
  if (lbstyle == SHIFT) balance->shift_setup(bstr,nitermax,thresh);
  balance->options(iarg,narg,arg);
  wtflag = balance->wtflag;

  if (balance->varflag && nevery == 0)
    error->all(FLERR,"Fix balance nevery = 0 cannot be used with weight var");

  // create instance of Irregular class

  irregular = new Irregular(lmp);

  // only force reneighboring if nevery > 0

  if (nevery) force_reneighbor = 1;
  lastbalance = -1;

  // compute initial outputs

  itercount = 0;
  pending = 0;
  imbfinal = imbprev = maxloadperproc = 0.0;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  delete balance;
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::post_create()
{
  if (wtflag) balance->weight_storage(id);
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  const int my_ifix = modify->find_fix(id);
  const int nfix = modify->nfix;
  Fix **fix = modify->fix;
  for(int ifix = my_ifix+1; ifix < nfix; ifix++)
  {
      if (strcmp(fix[ifix]->style, "insert") == 0)
      {
          char errstr[200];
          sprintf(errstr,"Fix %s with id %s has to come before fix %s", fix[ifix]->style, fix[ifix]->id, style);
          error->fix_error(FLERR,this,errstr);
      }
  }

  execute_flag = true;

  balance->init_imbalance(1);
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup(int vflag)
{
  // compute final imbalance factor if setup_pre_exchange() invoked balancer
  // this is called at end of run setup, before output

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // insure atoms are in current box & update box via shrink-wrap
  // has to be be done before rebalance() invokes Irregular::migrate_all()
  //   since it requires atoms be inside simulation box
  //   even though pbc() will be done again in Verlet::run()
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshold exceeded

  balance->set_weights();
  imbnow = balance->imbalance_factor(maxloadperproc);
  if (imbnow > thresh && execute_flag) rebalance(true);

  execute_flag = true;

  // next_reneighbor = next time to force reneighboring

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  // return if not a rebalance timestep

  if (nevery && update->ntimestep < next_reneighbor) return;

  // do not allow rebalancing twice on same timestep
  // even if wanted to, can mess up elapsed time in ImbalanceTime

  if (update->ntimestep == lastbalance) return;
  lastbalance = update->ntimestep;

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshold exceeded
  // if weight variable is used, wrap weight setting in clear/add compute

  if (balance->varflag) modify->clearstep_compute();
  balance->set_weights();
  if (balance->varflag) modify->addstep_compute(update->ntimestep + nevery);

  imbnow = balance->imbalance_factor(maxloadperproc);

  if (imbnow > thresh) rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   compute final imbalance factor based on nlocal after comm->exchange()
   only do this if rebalancing just occured
------------------------------------------------------------------------- */

void FixBalance::pre_neighbor()
{
  if (!pending) return;
  imbfinal = balance->imbalance_factor(maxloadperproc);
  pending = 0;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::rebalance(const bool setupflag)
{
  imbprev = imbnow;

  // invoke balancer and reset comm->uniform flag

  int **sendproc = NULL;
  if (lbstyle == SHIFT) {
    itercount = balance->shift();
    comm->set_layout(LAYOUT_NONUNIFORM);
  } else if (lbstyle == BISECTION) {
    sendproc = balance->bisection();
    comm->set_layout(LAYOUT_TILED);
  }

  // output of new decomposition

  if (balance->outflag) balance->dumpout(update->ntimestep);

  // reset proc sub-domains
  // check and warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in exchange()

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();
  domain->subbox_too_small_check(neighbor->skin);

  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far
  // else allow caller's comm->exchange() to do it

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (lbstyle == BISECTION) irregular->migrate_all(setupflag,0,1,sendproc);
  else if (irregular->migrate_check())
    irregular->migrate_all(setupflag);
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // invoke KSpace setup_grid() to adjust to new proc sub-domains

  if (kspace_flag) force->kspace->setup_grid();

  // pending triggers pre_neighbor() to compute final imbalance factor
  // can only be done after atoms migrate in comm->exchange()

  pending = 1;
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return maxloadperproc;
  if (i == 1) return (double) itercount;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = irregular->memory_usage();
  if (balance->rcb) bytes += balance->rcb->memory_usage();
  return bytes;
}
