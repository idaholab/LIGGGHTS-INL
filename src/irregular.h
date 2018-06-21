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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_IRREGULAR_H
#define LMP_IRREGULAR_H

#include "pointers.h"
#include "parallel_base.h"
#include <list>

namespace LAMMPS_NS {

class Irregular : protected Pointers {
 public:

#if defined(LMP_QSORT)
  // static variable across all Irregular objects, for qsort callback

  static int *proc_recv_copy;
#endif

  Irregular(class LAMMPS *);
  ~Irregular();
  void migrate_all(const bool setupflag = false, int sortflag = 0, int preassign = 0,
                     int **procassign = NULL);
  int migrate_check();
  int create_data(int, int *, int sortflag = 0);
  void exchange_data(char *, int, char *);
  void destroy_data();
  bigint memory_usage();

  void migrate_pb(ParallelBase* pb, int sortflag = 0, int preassign = 0,
                  int *procassign = NULL);

 private:
  void migrate_atoms(int sortflag = 0, int preassign = 0,
                     int *procassign = NULL);

  int migrate_check_atom();
  int migrate_check_pb(ParallelBase* pb);

  int me,nprocs;
  int triclinic;
  int map_style;

  int maxsend,maxrecv;              // size of buf send/recv in # of doubles
  double *buf_send,*buf_recv;       // bufs used in migrate_atoms
  int maxdbuf;                      // size of double buf in bytes
  double *dbuf;                     // double buf for largest single atom send
  int maxbuf;                       // size of char buf in bytes
  char *buf;                        // char buf for largest single data send

  int *mproclist,*msizes;           // persistent vectors in migrate_atoms
  int maxlocal;                     // allocated size of mproclist and msizes

  int *work1,*work2;                // work vectors

  // plan params for irregular communication of atoms or datums
  // no params refer to atoms/data copied to self

  int nsend_proc;            // # of messages to send
  int nrecv_proc;            // # of messages to recv
  int sendmax_proc;          // # of doubles/datums in largest send message
  int *proc_send;            // list of procs to send to
  int *num_send;             // # of atoms/datums to send to each proc
  int *index_send;           // list of which atoms/datums to send to each proc
  int *proc_recv;            // list of procs to recv from
  MPI_Request *request;      // MPI requests for posted recvs
  MPI_Status *status;        // MPI statuses for WaitAll

  // extra plan params plan for irregular communication of atoms
  // no params refer to atoms copied to self

  int *length_send;          // # of doubles to send to each proc
  int *length_recv;          // # of doubles to recv from each proc
  int *offset_send;          // where each atom starts in send buffer

  // extra plan params plan for irregular communication of datums
  // 2 self params refer to data copied to self

  int *num_recv;             // # of datums to recv from each proc
  int num_self;              // # of datums to copy to self
  int *index_self;           // list of which datums to copy to self

  // parallel base pointers
  std::list<ParallelBase*> pb_;

  // private methods

  int create_atom(int, int *, int *, int);
  void exchange_atom(double *, int *, double *);
  void destroy_atom();

  int binary(double, int, double *);

  void grow_send(int,int);          // reallocate send buffer
  void grow_recv(int);              // free/allocate recv buffer
};

}

#endif

/* ERROR/WARNING messages:

*/
