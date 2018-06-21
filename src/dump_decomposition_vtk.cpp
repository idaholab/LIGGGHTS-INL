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

    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK)

#include <string.h>
#include "dump_decomposition_vtk.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "comm.h"
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkVoxel.h>
#include <vtkPoints.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>

// include last to ensure correct macros
#include "domain_definitions.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpDecompositionVTK::DumpDecompositionVTK(LAMMPS *lmp, int narg, char **arg) :
    Dump(lmp, narg, arg),
    DumpVTK(lmp),
    xdata(NULL),
    ydata(NULL),
    zdata(NULL),
    lasttimestep(-1),
    filecurrent(NULL)
{
    if (narg != 5)
        error->all(FLERR,"Illegal dump decomposition command");

    //INFO: CURRENTLY ONLY PROC 0 writes

    format_default = NULL;

    //number of properties written out in one line with buff
    size_one=1;  //dont use buff

    std::list<int> allowed_extensions;
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTK);
    allowed_extensions.push_back(VTK_FILE_FORMATS::VTU);
    DumpVTK::identify_file_type(filename, allowed_extensions, style, multiproc, nclusterprocs, filewriter, fileproc, world, clustercomm);

    xdata = new double[2*comm->nprocs];
    ydata = new double[2*comm->nprocs];
    zdata = new double[2*comm->nprocs];
    all_nlocal_ = new int[comm->nprocs];
}

/* ---------------------------------------------------------------------- */

DumpDecompositionVTK::~DumpDecompositionVTK()
{
    delete []xdata;
    delete []ydata;
    delete []zdata;
    delete []all_nlocal_;
    delete []filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::init_style()
{
    if (domain->triclinic == 1)
        error->all(FLERR,"Can not perform dump decomposition for triclinic box");
    if (binary)
        error->all(FLERR,"Can not perform dump decomposition in binary mode");

    // default format not needed

    delete [] format;
    format = new char[150];

    // open single file, one time only

    if (multifile == 0)
        openfile();

    delete []xdata;
    delete []ydata;
    delete []zdata;

    xdata = new double[2*comm->nprocs];
    ydata = new double[2*comm->nprocs];
    zdata = new double[2*comm->nprocs];
}

/* ---------------------------------------------------------------------- */

int DumpDecompositionVTK::modify_param(int narg, char **arg)
{
    const int mvtk = DumpVTK::modify_param(narg, arg);
    if (mvtk > 0)
        return mvtk;
    return 0;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_header(bigint ndump)
{
    if (comm->me!=0)
        return;
}

/* ---------------------------------------------------------------------- */

int DumpDecompositionVTK::count()
{
    if (comm->me!=0)
        return 0;
    return 1;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::pack(int *ids)
{
    
    double to_send[2];
    to_send[0] = domain->sublo[0];
    to_send[1] = domain->subhi[0];
    MPI_Gather(to_send, 2, MPI_DOUBLE, xdata, 2, MPI_DOUBLE, 0 ,world);
    to_send[0] = domain->sublo[1];
    to_send[1] = domain->subhi[1];
    MPI_Gather(to_send, 2, MPI_DOUBLE, ydata, 2, MPI_DOUBLE, 0 ,world);
    to_send[0] = domain->sublo[2];
    to_send[1] = domain->subhi[2];
    MPI_Gather(to_send, 2, MPI_DOUBLE, zdata, 2, MPI_DOUBLE, 0 ,world);
    MPI_Gather(&(atom->nlocal), 1, MPI_INT, all_nlocal_, 1, MPI_INT, 0, world);

    return;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_data(int n, double *mybuf)
{
    
    if (comm->me!=0)
        return;

    //ensure it is only written once in multi-proc (work-around)
    if(lasttimestep==update->ntimestep)
        return;
    lasttimestep=update->ntimestep;
    DumpVTK::setFileCurrent(filecurrent, filename, multifile, padflag);

    vtkSmartPointer<vtkUnstructuredGrid> procGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    procGrid->Allocate(comm->nprocs, 1);

    vtkSmartPointer<vtkPoints> procPoints = vtkSmartPointer<vtkPoints>::New();
    procPoints->SetNumberOfPoints(8*comm->nprocs);

    vtkSmartPointer<vtkIntArray> nlocal = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> proc = vtkSmartPointer<vtkIntArray>::New();
    nlocal->SetName("nlocal");
    proc->SetName("proc");

    for (int i = 0; i < comm->nprocs; i++)
    {
        procPoints->InsertPoint(8*i+0, xdata[2*i+0], ydata[2*i+0], zdata[2*i+0]);
        procPoints->InsertPoint(8*i+1, xdata[2*i+0], ydata[2*i+0], zdata[2*i+1]);
        procPoints->InsertPoint(8*i+2, xdata[2*i+0], ydata[2*i+1], zdata[2*i+0]);
        procPoints->InsertPoint(8*i+3, xdata[2*i+0], ydata[2*i+1], zdata[2*i+1]);
        procPoints->InsertPoint(8*i+4, xdata[2*i+1], ydata[2*i+0], zdata[2*i+0]);
        procPoints->InsertPoint(8*i+5, xdata[2*i+1], ydata[2*i+0], zdata[2*i+1]);
        procPoints->InsertPoint(8*i+6, xdata[2*i+1], ydata[2*i+1], zdata[2*i+0]);
        procPoints->InsertPoint(8*i+7, xdata[2*i+1], ydata[2*i+1], zdata[2*i+1]);

        vtkSmartPointer<vtkVoxel> procCell = vtkSmartPointer<vtkVoxel>::New();
        procCell->GetPointIds()->SetId(0, 8*i+0);
        procCell->GetPointIds()->SetId(1, 8*i+1);
        procCell->GetPointIds()->SetId(2, 8*i+2);
        procCell->GetPointIds()->SetId(3, 8*i+3);
        procCell->GetPointIds()->SetId(4, 8*i+4);
        procCell->GetPointIds()->SetId(5, 8*i+5);
        procCell->GetPointIds()->SetId(6, 8*i+6);
        procCell->GetPointIds()->SetId(7, 8*i+7);

        procGrid->InsertNextCell(procCell->GetCellType(), procCell->GetPointIds());

        nlocal->InsertNextValue(all_nlocal_[i]);
        proc->InsertNextValue(i);
    }
    procGrid->SetPoints(procPoints);
    procGrid->GetCellData()->AddArray(nlocal);
    procGrid->GetCellData()->AddArray(proc);

    //vtkSmartPointer<vtkDataObject> procGrid = mbSet->GetBlock(1);
    if (vtk_file_format_ == VTK_FILE_FORMATS::VTU)
        DumpVTK::write_vtu(procGrid, filecurrent);
    else
        DumpVTK::write_vtk_unstructured_grid(procGrid, filecurrent);
}

#endif
