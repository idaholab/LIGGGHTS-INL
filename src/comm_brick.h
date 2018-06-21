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

    Arno Mayrhofer (DCS Computing GmbH)

    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2017-     DCS Computing GmbH, Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_COMM_BRICK_H
#define LMP_COMM_BRICK_H

#include "comm.h"
#include "atom.h"

#include <algorithm>

namespace LAMMPS_NS {

class CommBrick : public Comm {
 public:
  CommBrick(class LAMMPS *);
  CommBrick(class LAMMPS *, class Comm *);
  CommBrick(class Comm *, ParallelBase *);
  virtual ~CommBrick();

  virtual void init();
  virtual void setup();                        // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0);    // forward comm of atom coords
  virtual void reverse_comm();                 // reverse comm of forces
  virtual void exchange();                     // move atoms to new procs
  virtual void borders();                      // setup list of atoms to comm

  virtual void forward_comm_pair(class Pair *);    // forward comm from a Pair
  virtual void reverse_comm_pair(class Pair *);    // reverse comm from a Pair
  virtual void forward_comm_fix(class Fix *);
                                                   // forward comm from a Fix
  virtual void reverse_comm_fix(class Fix *);
                                                   // reverse comm from a Fix
  virtual void forward_comm_variable_fix(class Fix *);
                                     // variable size forward comm from a Fix
  virtual void forward_comm_compute(class Compute *);  // forward from a Compute
  virtual void reverse_comm_compute(class Compute *);  // reverse from a Compute
  virtual void forward_comm_dump(class Dump *);    // forward comm from a Dump
  virtual void reverse_comm_dump(class Dump *);    // reverse comm from a Dump

  void forward_comm_array(int, double **);         // forward comm of array
  virtual bigint memory_usage();

 protected:
  // NOTE: init_buffers is called from a constructor and must not be made virtual
  void init_buffers();

  int updown(int, int, int, double, int, double *);
                                            // compare cutoff to procs
  bool use_gran_opt();
  bool decide(int i,int dim,double lo,double hi,int ineed);
  bool decide_wedge(int i,int dim,double lo,double hi,int ineed);
  int get_type(const int i) const
  { return pb_ ? std::max(pb_->get_type(i),1) : atom->type[i]; }

  virtual void grow_list(int, int);         // reallocate one sendlist
  virtual void grow_swap(int);              // grow swap and multi arrays
  virtual void allocate_swap(int);          // allocate swap arrays
  virtual void allocate_multi(int);         // allocate multi arrays
  virtual void free_swap();                 // free swap arrays
  virtual void free_multi();                // free multi arrays
  virtual void exchangeEventsRecorder();    // Recorder for Exchange events
  virtual void exchangeEventsCorrector();   // Corrects receiving process ids

  int nswap;                        // # of swaps to perform = sum of maxneed
  int recvneed[3][2];               // # of procs away I recv atoms from
  int sendneed[3][2];               // # of procs away I send atoms to
  int maxneed[3];                   // max procs away any proc needs, per dim
  int maxswap;                      // max # of swaps memory is allocated for
  int *sendnum,*recvnum;            // # of atoms to send/recv in each swap
  int *sendproc,*recvproc;          // proc to send/recv to/from at each swap
  int *size_forward_recv;           // # of values to recv in each forward comm
  int *size_reverse_send;           // # to send in each reverse comm
  int *size_reverse_recv;           // # to recv in each reverse comm
  double *slablo,*slabhi;           // bounds of slab to send at each swap
  double **multilo,**multihi;       // bounds of slabs for multi-type swap
  double **cutghostmulti;           // cutghost on a per-type basis
  int *pbc_flag;                    // general flag for sending atoms thru PBC
  int **pbc;                        // dimension flags for PBC adjustments

  int *firstrecv;                   // where to put 1st recv atom in each swap
  int **sendlist;                   // list of atoms to send in each swap
  int *maxsendlist;                 // max size of send list for each swap

  int smax,rmax;             // max size in atoms of single borders send/recv

  // Domain wedge variables

  class DomainWedge *dw_;
  int ia,iphi;
  
  double nleft[2],nright[2],pleft[2],pright[2],c[2];
};

}

#endif

/* ERROR/WARNING messages:

E: Cannot change to comm_style brick from tiled layout

Self-explanatory.

*/
