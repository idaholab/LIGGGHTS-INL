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

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK

#include <string.h>
#include "dump_mesh_vtm.h"
#include "tri_mesh.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "comm.h"
#include <stdint.h>
#include <vtkXMLMultiBlockDataWriter.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpMeshVTM::DumpMeshVTM(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    filecurrent(NULL),
    multiname_ex(NULL),
    dumpMesh(NULL)
{
    if (narg < 7)
        error->all(FLERR,"Illegal dump mesh/vtm command");

    //INFO: CURRENTLY ONLY PROC 0 writes

    format_default = NULL;

    if (!vtkMultiProcessController::GetGlobalController())
    {
        vtkMPIController *vtkController = vtkMPIController::New();
        vtkController->Initialize();
        vtkMultiProcessController::SetGlobalController(vtkController);
    }
    vtkMPIController * controller = getLocalController();

    std::list<int> allowed_extensions;
    allowed_extensions.push_back(VTK_FILE_FORMATS::PVTM);
    DumpVTK::identify_file_type(filename, allowed_extensions, style, multiproc, nclusterprocs, filewriter, fileproc, world, clustercomm);

    int ioptional = 5;
    dumpMesh = new DumpMesh(lmp, nclusterprocs, multiproc, filewriter, fileproc, controller);
    ioptional += dumpMesh->parse_parameters(narg-ioptional, &(arg[ioptional]));

    if (ioptional < narg)
        error->all(FLERR,"Invalid attribute in dump mesh/vtm command");

    char *ptr = strchr(filename,'%');
    if (ptr) {
      multiname_ex = new char[strlen(filename) + 16];
      *ptr = '\0';
      sprintf(multiname_ex,"%s_%d%s",filename,me,ptr+1);
      *ptr = '%';
    }
}

/* ---------------------------------------------------------------------- */

DumpMeshVTM::~DumpMeshVTM()
{
    if (filecurrent)
        delete [] filecurrent;
    delete dumpMesh;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::init_style()
{
    size_one = dumpMesh->init_style();
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::write_header(bigint ndump)
{
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTM::count()
{
    return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::pack(int *ids)
{
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::setFileCurrent()
{
    DumpVTK::setFileCurrent(filecurrent, filename, multifile, padflag);
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::write()
{
    write_data(0, NULL);
}

/* ---------------------------------------------------------------------- */

void DumpMeshVTM::write_data(int n, double *mybuf)
{
    setFileCurrent();

    vtkSmartPointer<vtkMultiBlockDataSet> mbSet = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    dumpMesh->prepare_mbSet(mbSet);

    if (!filewriter)
        return;

    DumpVTK::write_vtm(mbSet, filecurrent);
}

/* ---------------------------------------------------------------------- */

int DumpMeshVTM::modify_param(int narg, char **arg)
{
    const int mvtk = DumpVTK::modify_param(narg, arg);
    if (mvtk > 0)
        return mvtk;

    return 0;
}

#endif // LAMMPS_VTK
