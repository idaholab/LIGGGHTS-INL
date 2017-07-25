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
    (if not contributing author is listed, this file has been contributed
    by the core developer)
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2017-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <mpi.h>
#include "compute_bond_counter.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBondCounter::ComputeBondCounter(LAMMPS *lmp, int &iarg, int narg, char **arg) :
    Compute(lmp, iarg, narg, arg),
    pair_created(0),
    pair_broken(0),
    pair_total(0),
    wall_created(0),
    wall_broken(0),
    wall_total(0)
{
    if (narg != iarg) error->all(FLERR,"Illegal compute bond/counter command");

    vector_flag = 1;
    size_vector = 6;
    vector = new double[size_vector];

    if (modify->find_compute_style_strict("bond/counter", 0))
        error->all(FLERR, "Only one compute of type bond/counter allowed");
}

/* ---------------------------------------------------------------------- */

ComputeBondCounter::~ComputeBondCounter()
{
    delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeBondCounter::init()
{
    invoked_vector = update->ntimestep;

    pair_created = 0;
    pair_broken = 0;
    pair_total = 0;
    wall_created = 0;
    wall_broken = 0;
    wall_total = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeBondCounter::compute_vector()
{
    invoked_vector = update->ntimestep;

    unsigned int tmp = pair_created;
    MPI_Reduce(&tmp, &pair_created, 1, MPI_UNSIGNED, MPI_SUM, 0, world);
    tmp = pair_broken;
    MPI_Reduce(&tmp, &pair_broken, 1, MPI_UNSIGNED, MPI_SUM, 0, world);
    tmp = pair_total;
    MPI_Reduce(&tmp, &pair_total, 1, MPI_UNSIGNED, MPI_SUM, 0, world);
    tmp = wall_created;
    MPI_Reduce(&tmp, &wall_created, 1, MPI_UNSIGNED, MPI_SUM, 0, world);
    tmp = wall_broken;
    MPI_Reduce(&tmp, &wall_broken, 1, MPI_UNSIGNED, MPI_SUM, 0, world);
    tmp = wall_total;
    MPI_Reduce(&tmp, &wall_total, 1, MPI_UNSIGNED, MPI_SUM, 0, world);

    // since count was computed on the timestep after the previous output
    // it does not necessarily reflect the current number and needs to be
    // corrected by the bonds created and destroyed since then
    pair_total += pair_created - pair_broken;
    wall_total += wall_created - wall_broken;

    vector[0] = static_cast<double>(pair_created);
    vector[1] = static_cast<double>(pair_broken);
    vector[2] = static_cast<double>(pair_total);
    vector[3] = static_cast<double>(wall_created);
    vector[4] = static_cast<double>(wall_broken);
    vector[5] = static_cast<double>(wall_total);

    pair_created = 0;
    pair_broken = 0;
    pair_total = 0;
    wall_created = 0;
    wall_broken = 0;
    wall_total = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeBondCounter::bond_created(const bool is_wall, const int i, const int j)
{
    if (is_wall)
        wall_created++;
    else if (force->newton_bond || j < atom->nlocal || atom->tag[i] < atom->tag[j])
        pair_created++;
}

/* ---------------------------------------------------------------------- */

void ComputeBondCounter::bond_broken(const bool is_wall, const int i, const int j)
{
    if (is_wall)
        wall_broken++;
    else if (force->newton_bond || j < atom->nlocal || atom->tag[i] < atom->tag[j])
        pair_broken++;
}

/* ---------------------------------------------------------------------- */

void ComputeBondCounter::bond_count(const bool is_wall, const int i, const int j)
{
    // count only if recently reset
    if (invoked_vector == update->ntimestep-1)
    {
        if (is_wall)
            wall_total++;
        else if (force->newton_bond || j < atom->nlocal || atom->tag[i] < atom->tag[j])
            pair_total++;
    }
}
