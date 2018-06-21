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

#ifdef DUMP_CLASS

DumpStyle(decomposition/vtk,DumpDecompositionVTK)

#else

#ifndef LMP_DUMP_DECOMPOSITION_VTK_H
#define LMP_DUMP_DECOMPOSITION_VTK_H

#include "dump.h"
#include "dump_vtk.h"

namespace LAMMPS_NS
{

class DumpDecompositionVTK : public Dump, public DumpVTK
{
public:
    DumpDecompositionVTK(LAMMPS *, int, char**);
    ~DumpDecompositionVTK();
    void init_style();

private:
    int modify_param(int, char **);
    void write_header(bigint);
    int count();
    void pack(int *);
    void write_data(int, double *);

    double *xdata;
    double *ydata;
    double *zdata;
    int *all_nlocal_;

    int lasttimestep;
    char *filecurrent;
};

}

#endif
#endif

#endif // LAMMPS_VTK
