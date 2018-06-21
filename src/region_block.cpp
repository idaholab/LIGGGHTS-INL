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

#include <stdlib.h>
#include <string.h>
#include "region_block.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "variable.h"
#include "update.h"

// include last to ensure correct macros
#include "domain_definitions.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegBlock::RegBlock(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg),
  xloisvar(false),
  xlovar(0),
  xlostr(NULL),
  xhiisvar(false),
  xhivar(0),
  xhistr(NULL),
  yloisvar(false),
  ylovar(0),
  ylostr(NULL),
  yhiisvar(false),
  yhivar(0),
  yhistr(NULL),
  zloisvar(false),
  zlovar(0),
  zlostr(NULL),
  zhiisvar(false),
  zhivar(0),
  zhistr(NULL),
  xlo(0.0),
  xhi(0.0),
  ylo(0.0),
  yhi(0.0),
  zlo(0.0),
  zhi(0.0)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
    else if (domain->triclinic == 0) xlo = domain->boxlo[0];
    else xlo = domain->boxlo_bound[0];
  }
  else
  {
    if (strncmp(arg[2], "v_", 2) == 0)
    {
      xlostr = new char[strlen(arg[2])-2+1];
      strcpy(xlostr, &arg[2][2]);
      xloisvar = true;
      varshape = 1;
    }
    else
      xlo = xscale*force->numeric(FLERR,arg[2]);
  }

  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else if (domain->triclinic == 0) xhi = domain->boxhi[0];
    else xhi = domain->boxhi_bound[0];
  }
  else
  {
    if (strncmp(arg[3], "v_", 2) == 0)
    {
      xhistr = new char[strlen(arg[3])-2+1];
      strcpy(xhistr, &arg[3][2]);
      xhiisvar = true;
      varshape = 1;
    }
    else
      xhi = xscale*force->numeric(FLERR,arg[3]);
  }

  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
    else if (domain->triclinic == 0) ylo = domain->boxlo[1];
    else ylo = domain->boxlo_bound[1];
  }
  else
  {
    if (strncmp(arg[4], "v_", 2) == 0)
    {
      ylostr = new char[strlen(arg[4])-2+1];
      strcpy(ylostr, &arg[4][2]);
      yloisvar = true;
      varshape = 1;
    }
    else
      ylo = yscale*force->numeric(FLERR,arg[4]);
  }

  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[5],"INF") == 0) yhi = BIG;
    else if (domain->triclinic == 0) yhi = domain->boxhi[1];
    else yhi = domain->boxhi_bound[1];
  }
  else
  {
    if (strncmp(arg[5], "v_", 2) == 0)
    {
      yhistr = new char[strlen(arg[5])-2+1];
      strcpy(yhistr, &arg[5][2]);
      yhiisvar = true;
      varshape = 1;
    }
    else
      yhi = yscale*force->numeric(FLERR,arg[5]);
  }

  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else if (domain->triclinic == 0) zlo = domain->boxlo[2];
    else zlo = domain->boxlo_bound[2];
  }
  else
  {
    if (strncmp(arg[6], "v_", 2) == 0)
    {
      zlostr = new char[strlen(arg[6])-2+1];
      strcpy(zlostr, &arg[6][2]);
      zloisvar = true;
      varshape = 1;
    }
    else
      zlo = zscale*force->numeric(FLERR,arg[6]);
  }

  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else if (domain->triclinic == 0) zhi = domain->boxhi[2];
    else zhi = domain->boxhi_bound[2];
  }
  else
  {
    if (strncmp(arg[7], "v_", 2) == 0)
    {
      zhistr = new char[strlen(arg[7])-2+1];
      strcpy(zhistr, &arg[7][2]);
      zhiisvar = true;
      varshape = 1;
    }
    else
      zhi = zscale*force->numeric(FLERR,arg[7]);
  }

  // error check
  if (varshape)
  {
    variable_check();
    shape_update();
  }

  if (xlo > xhi || ylo > yhi || zlo > zhi)
    error->all(FLERR,"Illegal region block command");

  // extent of block

  if (interior) {
    bboxflag = 1;
    extent_xlo = xlo;
    extent_xhi = xhi;
    extent_ylo = ylo;
    extent_yhi = yhi;
    extent_zlo = zlo;
    extent_zhi = zhi;
  } else bboxflag = 0;

  // particle could be close to all 6 planes

  cmax = 6;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegBlock::~RegBlock()
{
  delete [] contact;
  if (xlostr)
    delete[] xlostr;
  if (xhistr)
    delete[] xhistr;
  if (ylostr)
    delete[] ylostr;
  if (yhistr)
    delete[] yhistr;
  if (zlostr)
    delete[] zlostr;
  if (zhistr)
    delete[] zhistr;
}

/* ---------------------------------------------------------------------- */

void RegBlock::init()
{
  Region::init();
  if (varshape)
    variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegBlock::inside(double x, double y, double z)
{
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of block
   can be one contact for each of 6 faces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

int RegBlock::surface_interior(double *x, double cutoff)
{
  double delta;

  // x is exterior to block

  if (x[0] < xlo || x[0] > xhi || x[1] < ylo || x[1] > yhi ||
      x[2] < zlo || x[2] > zhi) return 0;

  // x is interior to block or on its surface

  int n = 0;

  delta = x[0] - xlo;
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delx = delta;
    contact[n].dely = contact[n].delz = 0.0;
    n++;
  }
  delta = xhi - x[0];
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delx = -delta;
    contact[n].dely = contact[n].delz = 0.0;
    n++;
  }

  delta = x[1] - ylo;
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].dely = delta;
    contact[n].delx = contact[n].delz = 0.0;
    n++;
  }
  delta = yhi - x[1];
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].dely = -delta;
    contact[n].delx = contact[n].delz = 0.0;
    n++;
  }

  delta = x[2] - zlo;
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delz = delta;
    contact[n].delx = contact[n].dely = 0.0;
    n++;
  }
  delta = zhi - x[2];
  if (delta < cutoff) {
    contact[n].r = delta;
    contact[n].delz = -delta;
    contact[n].delx = contact[n].dely = 0.0;
    n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of block
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

int RegBlock::surface_exterior(double *x, double cutoff)
{
  double xp,yp,zp;

  // x is far enough from block that there is no contact
  // x is interior to block

  if (x[0] <= xlo-cutoff || x[0] >= xhi+cutoff ||
      x[1] <= ylo-cutoff || x[1] >= yhi+cutoff ||
      x[2] <= zlo-cutoff || x[2] >= zhi+cutoff) return 0;
  if (x[0] > xlo && x[0] < xhi && x[1] > ylo && x[1] < yhi &&
      x[2] > zlo && x[2] < zhi) return 0;

  // x is exterior to block or on its surface
  // xp,yp,zp = point on surface of block that x is closest to
  //            could be edge or corner pt of block
  // do not add contact point if r >= cutoff

  if (x[0] < xlo) xp = xlo;
  else if (x[0] > xhi) xp = xhi;
  else xp = x[0];
  if (x[1] < ylo) yp = ylo;
  else if (x[1] > yhi) yp = yhi;
  else yp = x[1];
  if (x[2] < zlo) zp = zlo;
  else if (x[2] > zhi) zp = zhi;
  else zp = x[2];

  add_contact(0,x,xp,yp,zp);
  if (contact[0].r < cutoff) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void RegBlock::variable_check()
{
  if (xloisvar)
  {
      xlovar = input->variable->find(xlostr);
      if (xlovar < 0)
        error->all(FLERR,"Variable name for region block variable xlo does not exist");
      if (!input->variable->equalstyle(xlovar))
        error->all(FLERR,"Variable for region block variable xlo is invalid style");
  }
  if (xhiisvar)
  {
    xhivar = input->variable->find(xhistr);
    if (xhivar < 0)
      error->all(FLERR,"Variable name for region block variable xhi does not exist");
    if (!input->variable->equalstyle(xhivar))
      error->all(FLERR,"Variable for region block variable xhi is invalid style");
  }

  if (yloisvar)
  {
    ylovar = input->variable->find(ylostr);
    if (ylovar < 0)
      error->all(FLERR,"Variable name for region block variable ylo does not exist");
    if (!input->variable->equalstyle(ylovar))
      error->all(FLERR,"Variable for region block variable ylo is invalid style");
  }
  if (yhiisvar)
  {
    yhivar = input->variable->find(yhistr);
    if (yhivar < 0)
      error->all(FLERR,"Variable name for region block variable yhi does not exist");
    if (!input->variable->equalstyle(yhivar))
      error->all(FLERR,"Variable for region block variable yhi is invalid style");
  }

  if (zloisvar)
  {
    zlovar = input->variable->find(zlostr);
    if (zlovar < 0)
      error->all(FLERR,"Variable name for region block variable zlo does not exist");
    if (!input->variable->equalstyle(zlovar))
      error->all(FLERR,"Variable for region block variable zlo is invalid style");
  }
  if (zhiisvar)
  {
    zhivar = input->variable->find(zhistr);
    if (zhivar < 0)
      error->all(FLERR,"Variable name for region block variable zhi does not exist");
    if (!input->variable->equalstyle(zhivar))
      error->all(FLERR,"Variable for region block variable zhi is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void RegBlock::shape_update()
{
  if (xloisvar)
    xlo = xscale * input->variable->compute_equal(xlovar);
  if (xhiisvar)
    xhi = xscale * input->variable->compute_equal(xhivar);
  if (xlo > xhi)
    error->one(FLERR,"Variable evaluation in region block x-axis gave bad value");
  if (yloisvar)
    ylo = yscale * input->variable->compute_equal(ylovar);
  if (yhiisvar)
    yhi = yscale * input->variable->compute_equal(yhivar);
  if (ylo > yhi)
    error->one(FLERR,"Variable evaluation in region block y-axis gave bad value");
  if (zloisvar)
    zlo = zscale * input->variable->compute_equal(zlovar);
  if (zhiisvar)
    zhi = zscale * input->variable->compute_equal(zhivar);
  if (zlo > zhi)
    error->one(FLERR,"Variable evaluation in region block z-axis gave bad value");
}
