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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
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

#ifndef LMP_COMM_H
#define LMP_COMM_H

#include "pointers.h"
#include "parallel_base.h"
#include <vector>
#include <list>

namespace LAMMPS_NS {

enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

class Comm : protected Pointers {
 friend class Info;

 public:
  int style;     // comm pattern: 0 = 6-way stencil, 1 = irregular tiling
  int mode;      // 0 = single cutoff, 1 = multi-type cutoff

  int me,nprocs;                    // proc info
  int ghost_velocity;               // 1 if ghost atoms have velocity, 0 if not
  double cutghost[3];               // cutoffs used for acquiring ghost atoms
  double cutghostuser;              // user-specified ghost cutoff (mode == 0)
  double *cutusermulti;            // per type user ghost cutoff (mode == 1)
  int recv_from_partition;          // recv proc layout from this partition
  int send_to_partition;            // send my proc layout to this partition
                                    // -1 if no recv or send
  int maxexchange_atom;             // max contribution to exchange from AtomVec
  int maxexchange_fix;              // max contribution to exchange from Fixes
  int nthreads;                     // OpenMP threads per MPI process

  // public settings specific to layout = UNIFORM, NONUNIFORM

  int procgrid[3];                  // procs assigned in each dim of 3d grid
  int user_procgrid[3];             // user request for procs in each dim
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs, 0/1 = left/right
  double *xsplit,*ysplit,*zsplit;   // fractional (0-1) sub-domain sizes
  int ***grid2proc;                 // which proc owns i,j,k loc in 3d grid

  //exchange events recorder
  bool               exchangeEvents;                 //main switch to record exchange events
  std::vector<int>   exchangeEventsLocalId;          //local ids that change processor
  std::vector<int>   exchangeEventsReceivingProcess; //process IDs that receive the atom
  std::vector<int>   exchangeEventsGlobalProblemIds; //list with global ids that change more than one dim

  // methods

  Comm(class LAMMPS *);
  virtual ~Comm();
  // NOTE: copy_arrays is called from a constructor and must not be made virtual
  void copy_pointers(class Comm *);
  void copy_arrays(class Comm *);
  virtual void init();
  void set(int, char **);         // set communication style; equivalent to LAMMPS comm_modify

  void set_processors(int, char **);      // set 3d processor grid attributes
  virtual void set_proc_grid(int outflag = 1); // setup 3d grid of procs

  virtual void setup() = 0;                       // setup 3d comm pattern
  virtual void forward_comm(int dummy = 0) = 0;   // forward comm of atom coords
  virtual void reverse_comm() = 0;                // reverse comm of forces
  virtual void exchange() = 0;                    // move atoms to new procs
  virtual void borders() = 0;                     // setup list of atoms to comm

  // forward/reverse comm from a Pair, Fix, Compute, Dump

  virtual void forward_comm_pair(class Pair *) = 0;    // forward comm from a Pair
  virtual void reverse_comm_pair(class Pair *) = 0;    // reverse comm from a Pair
  virtual void forward_comm_fix(class Fix *) = 0;      // forward comm from a Fix
  virtual void reverse_comm_fix(class Fix *) = 0;      // reverse comm from a Fix
  virtual void forward_comm_variable_fix(class Fix *) = 0; // variable-size variant
  virtual void forward_comm_compute(class Compute *) = 0;  // forward from a Compute
  virtual void reverse_comm_compute(class Compute *) = 0;  // reverse from a Compute
  virtual void forward_comm_dump(class Dump *) = 0;    // forward comm from a Dump
  virtual void reverse_comm_dump(class Dump *) = 0;    // reverse comm from a Dump

  // forward comm of an array
  // exchange of info on neigh stencil
  // set processor mapping options

  virtual void forward_comm_array(int, double **) = 0; // forward comm of array
  int binary(double, int, double *);

  // map a point to a processor, based on current decomposition

  virtual void coord2proc_setup() {}
  virtual int coord2proc(const double *const, int &, int &, int &);

  // memory usage

  virtual bigint memory_usage() = 0;

  // non-virtual functions common to all Comm styles

  void ring(int, int, void *, int, void (*)(int, char *),   // ring comm
            void *, int self = 1);
  int read_lines_from_file(FILE *, int, int, char *);  // read/bcast file lines

  double* get_buf_send()
  { return buf_send; }
  double* get_buf_recv()
  { return buf_recv; }

  double* check_grow_send(const int n, const int flag, const int requestextra = -1); // reallocate send buffer if required

  void register_subcomm(Comm *comm)
  { subcomms_.push_back(comm); }

  double get_mysplit(const int i, const int j) const
  { return mysplit[i][j]; }

  void set_mysplit(const double _mysplit[3][2]);
  void set_rcbnew(const int _rcbnew);
  void set_rcbcutfrac(const double _rcbcutfrac);
  void set_rcbcutdim(const int _rcbdim);
  void set_sizes(const int _size_border, const int _size_forward, const int _size_reverse, const int _maxforward, const int _maxreverse);
  void set_layout(const int _layout);

  int get_layout()
  { return layout; }

  // for tiled only, copies exchproc info from main comm to sub comms
  virtual void copy_exchange_info_from_main_comm()
  {}

 protected:
  void grow_send(int, int, const int requestextra = -1); // reallocate send buffer
  void grow_recv(int);                 // free/allocate recv buffer

  int layout;    // LAYOUT_UNIFORM = equal-sized bricks
                 // LAYOUT_NONUNIFORM = logical bricks, but diff sizes via LB
                 // LAYOUT_TILED = general tiling, due to RCB LB

  int bordergroup;                  // only communicate this group in borders

  int map_style;                    // non-0 if global->local mapping is done
  int comm_x_only,comm_f_only;      // 1 if only exchange x,f in for/rev comm

  int size_forward;                 // # of per-atom datums in forward comm
  int size_reverse;                 // # of datums in reverse comm
  int size_border;                  // # of datums in forward border comm

  int maxforward,maxreverse;        // max # of datums in forward/reverse comm
  int maxexchange;                  // max # of datums/atom in exchange comm

  int gridflag;                     // option for creating 3d grid
  int mapflag;                      // option for mapping procs to 3d grid
  char xyz[4];                      // xyz mapping of procs to 3d grid
  char *customfile;                 // file with custom proc map
  char *outfile;                    // proc grid/map output file

  int otherflag;                    // 1 if this partition dependent on another
  int other_style;                  // style of dependency
  int other_procgrid[3];            // proc layout of another partition
  int other_coregrid[3];            // core layout of another partition
  int ncores;                       // # of cores per node
  int coregrid[3];                  // 3d grid of cores within a node
  int user_coregrid[3];             // user request for cores in each dim

  double *buf_send;                 // send buffer for all comm
  double *buf_recv;                 // recv buffer for all comm
  int maxsend,maxrecv;              // current size of send/recv buffer
  int bufextra;                     // extra space beyond maxsend in send buffer

  // public settings specific to layout = TILED

  int rcbnew;                       // 1 if just reset by rebalance, else 0
  double mysplit[3][2];             // fractional (0-1) bounds of my sub-domain
  double rcbcutfrac;                // fractional RCB cut by this proc
  int rcbcutdim;                    // dimension of RCB cut

  // Extension of comm beyond particles
  ParallelBase *pb_;

  // Other comms
  bool is_subcomm;
  std::list<Comm*> subcomms_;
};

}

#endif

/* ERROR/WARNING messages:

W: OMP_NUM_THREADS environment is not set.

This environment variable must be set appropriately to use the
USER-OMP pacakge.

E: Bad grid of processors

The 3d grid of processors defined by the processors command does not
match the number of processors LAMMPS is being run on.

E: Processor count in z must be 1 for 2d simulation

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid group in communicate command

Self-explanatory.

E: Communicate group != atom_modify first group

Self-explanatory.

E: Invalid cutoff in communicate command

Specified cutoff must be >= 0.0.

E: Specified processors != physical processors

The 3d grid of processors defined by the processors command does not
match the number of processors LAMMPS is being run on.

E: Cannot use processors part command without using partitions

See the command-line -partition switch.

E: Invalid partitions in processors part command

Valid partitions are numbered 1 to N and the sender and receiver
cannot be the same partition.

E: Sending partition in processors part command is already a sender

Cannot specify a partition to be a sender twice.

E: Receiving partition in processors part command is already a receiver

Cannot specify a partition to be a receiver twice.

E: Processors grid numa and map style are incompatible

Using numa for gstyle in the processors command requires using
cart for the map option.

E: Processors part option and grid style are incompatible

Cannot use gstyle numa or custom with the part option.

*/
