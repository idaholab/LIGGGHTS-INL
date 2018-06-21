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

#ifdef FIX_CLASS

FixStyle(balance,FixBalance)

#else

#ifndef LMP_FIX_BALANCE_H
#define LMP_FIX_BALANCE_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class LAMMPS *, int, char **);
  ~FixBalance();
  int setmask();
  void post_create();
  void init();
  void setup(int);
  void setup_pre_exchange();
  void pre_exchange();
  void pre_neighbor();
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  int nevery,lbstyle,nitermax;
  double thresh, stopthresh;

  bool execute_flag;

  char bstr[4];
  int wtflag;                   // 1 for weighted balancing

  double imbnow;                // current imbalance factor
  double imbprev;               // imbalance factor before last rebalancing
  double imbfinal;              // imbalance factor after last rebalancing
  double maxloadperproc;        // max load on any processor
  int itercount;                // iteration count of last call to Balance
  int kspace_flag;              // 1 if KSpace solver defined
  int pending;
  bigint lastbalance;           // last timestep balancing was attempted

  class Balance *balance;
  class Irregular *irregular;

  void rebalance(bool setupflag = false);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix balance string is invalid

The string can only contain the characters "x", "y", or "z".

E: Fix balance string is invalid for 2d simulation

The string cannot contain the letter "z".

E: Cannot open fix balance output file

Self-explanatory.

*/
