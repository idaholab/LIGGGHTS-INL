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

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "compute_pair_gran_local_bond.h"
#include "error.h"
#include "fix_wall_gran.h"
#include "pair_gran_proxy.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePairGranLocalBond::ComputePairGranLocalBond(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  ComputePairGranLocal(lmp, iarg, narg, arg),
  bond_history_offset_(0)
{

  // store everything by default , including history, fn, ft
  posflag =  idflag = fflag = fnflag = ftflag = torqueflag = torquenflag = torquetflag = histflag = 1;

  // do not store vel, area, heat flux, delta
  velflag = areaflag = deltaflag = heatflag = 0;

  wall = 0;
  if(strcmp(style,"wall/gran/local/bond") == 0) wall = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocalBond::init()
{
    ComputePairGranLocal::init();

    if (wall)
        bond_history_offset_ = static_cast<FixWallGran*>(fixwall)->get_history_offset("bond_contactflag");
    else
        bond_history_offset_ = static_cast<PairGranProxy*>(pairgran)->get_history_offset("bond_contactflag");

    //error->all(FLERR,"TODO here");
    //offset_bond_history_ = pairgran->offset_bond_history();

    if(bond_history_offset_ == -1)
        error->compute_error(FLERR,this,"requires a cohesion bond model to work with");
}

/* ---------------------------------------------------------------------- */

bool ComputePairGranLocalBond::decide_add(double *hist, double * &contact_pos)
{
    
    // ensure the right ones are picked, i.e. only those who actually have cohesive bond force, not just hertz etc

    if(hist[bond_history_offset_] > 0.5)
    {
        contact_pos = &hist[bond_history_offset_+2];
        return true;
    }
    return false;
}
