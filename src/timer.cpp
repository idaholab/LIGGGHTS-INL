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

#include <mpi.h>
#include "timer.h"
#include "memory.h"
#include "modify.h"

#ifdef _WIN32
#include <windows.h>
#include <stdint.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

using namespace LAMMPS_NS;

// Return the CPU time for the current process in seconds very
// much in the same way as MPI_Wtime() returns the wall time.

static double CPU_Time()
{
  double rv = 0.0;

#ifdef _WIN32

  // from MSD docs.
  FILETIME ct,et,kt,ut;
  union { FILETIME ft; uint64_t ui; } cpu;
  if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }

#else /* ! _WIN32 */

  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }

#endif /* ! _WIN32 */

  return rv;
}

/* ---------------------------------------------------------------------- */

Timer::Timer(LAMMPS *lmp) :
  Pointers(lmp),
  cpu_array(NULL)
{
  memory->create(array,TIME_N,"array");
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  memory->destroy(array);
  if (cpu_array)
    memory->destroy(cpu_array);
}

/* ---------------------------------------------------------------------- */

void Timer::init(bool save_cpu)
{
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;

  if(modify->timing) {
    for (int i = 0; i < modify->nfix; i++) modify->fix[i]->reset_time_recording();
  }

  if (save_cpu)
  {
      memory->create(cpu_array,TIME_N,"cpu_array");
      for (int i = 0; i < TIME_N; i++) cpu_array[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_time = MPI_Wtime();
  if (cpu_array)
      previous_cpu_time = CPU_Time();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(int which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  const double current_time = MPI_Wtime();
  array[which] += current_time - previous_time;
  previous_time = current_time;
  if (cpu_array)
  {
    const double current_cpu_time = CPU_Time();
    cpu_array[which] += current_cpu_time - previous_cpu_time;
    previous_cpu_time = current_cpu_time;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start(int which)
{
  MPI_Barrier(world);
  array[which] = MPI_Wtime();
  if (cpu_array)
    cpu_array[which] = CPU_Time();
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(int which)
{
  MPI_Barrier(world);
  const double current_time = MPI_Wtime();
  array[which] = current_time - array[which];
  if (cpu_array)
    cpu_array[which] = CPU_Time() - cpu_array[which];
}

/* ---------------------------------------------------------------------- */

double Timer::cpu_elapsed(int which)
{
  if (!cpu_array)
    return -1.0;
  const double current_time = CPU_Time();
  return (current_time - cpu_array[which]);
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(int which)
{
  const double current_time = MPI_Wtime();
  return (current_time - array[which]);
}
