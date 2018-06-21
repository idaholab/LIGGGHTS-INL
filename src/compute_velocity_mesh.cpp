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
#include "compute_velocity_mesh.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVelocityMesh::ComputeVelocityMesh(LAMMPS *lmp, int &iarg, int narg, char **arg) :
    Compute(lmp, iarg, narg, arg),
    fms(NULL)
{
    if (narg != iarg+1) error->all(FLERR,"Illegal compute velocity/mesh command");

    vector_flag = 1;
    size_vector = 6;
    extvector = 1;
    vector = new double[size_vector];

    // the only argument to this compute is a mesh id
    int ifix = modify->find_fix(arg[iarg++]);
    fms = dynamic_cast<FixMeshSurface*>(modify->fix[ifix]);
    if (!fms)
        error->all(FLERR, "compute velocity/mesh could not find mesh id");
    fms->triMesh()->set_store_vel();
    fms->triMesh()->set_store_omega();
}

/* ---------------------------------------------------------------------- */

ComputeVelocityMesh::~ComputeVelocityMesh()
{
    if (fms && fms->triMesh())
    {
        fms->triMesh()->unset_store_vel();
        fms->triMesh()->unset_store_omega();
    }
    delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeVelocityMesh::init()
{ }

/* ---------------------------------------------------------------------- */

void ComputeVelocityMesh::compute_vector()
{
    if (!fms)
        error->all(FLERR, "compute velocity/mesh referrs to a mesh that was removed");
    invoked_vector = update->ntimestep;

    fms->triMesh()->get_global_vel(&(vector[0]));
    fms->triMesh()->get_global_omega(&(vector[3]));
}
