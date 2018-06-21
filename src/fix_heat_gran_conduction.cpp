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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "fix_heat_gran_conduction.h"

#include "atom.h"
#include "compute_pair_gran_local.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_scalar_transport_equation.h"
#include "force.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "properties.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "pair_gran_proxy.h"
#include "global_properties.h"
#include <cmath>
#include <algorithm>

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "superquadric.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHeatGranCond::FixHeatGranCond(class LAMMPS *lmp, int narg, char **arg) :
  FixHeatGran(lmp, narg, arg),
  store_contact_data_(false),
  fix_conduction_contact_area_(0),
  fix_n_conduction_contacts_(0),
  fix_wall_heattransfer_coeff_(0),
  fix_wall_temperature_(0),
  conduction_contact_area_(0),
  n_conduction_contacts_(0),
  wall_heattransfer_coeff_(0),
  wall_temp_(0),
  temp_max_(-1.),
  heatPhy(lmp)
{
  iarg_ = 5;

  bool hasargs = true;
  while(iarg_ < narg && hasargs)
  {
    hasargs = false;

    if(strcmp(arg[iarg_],"store_contact_data") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'store_contact_data'");
      if(strcmp(arg[iarg_+1],"yes") == 0)
        store_contact_data_ = true;
      else if(strcmp(arg[iarg_+1],"no") == 0)
        store_contact_data_ = false;
      else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'store_contact_data'");
      iarg_ += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg_],"temp_max") == 0)
    {
          if (iarg_+2 > narg)
              error->fix_error(FLERR,this,"not enough arguments for keyword 'temp_max'");
          temp_max_ = force->numeric(FLERR,arg[iarg_+1]);
          if (temp_max_ < 0)
              error->fix_error(FLERR,this,"temp_max needs to be >= 0 K");
          iarg_ += 2;
          hasargs = true;
    }
  }

  heatPhy.parseArgs(iarg_,narg,arg); // create after own argument parsing

  if( iarg_ < narg && strcmp(style,"heat/gran/conduction") == 0 )
         error->fix_error(FLERR,this,"unknown keyword");
}

/* ---------------------------------------------------------------------- */

FixHeatGranCond::~FixHeatGranCond()
{

}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::post_create()
{
  FixHeatGran::post_create();

  // register contact storage
  fix_conduction_contact_area_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("contactAreaConduction","property/atom","scalar",0,0,this->style,false));
  if(!fix_conduction_contact_area_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="contactAreaConduction";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="contactAreaConduction";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_conduction_contact_area_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_n_conduction_contacts_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("nContactsConduction","property/atom","scalar",0,0,this->style,false));
  if(!fix_n_conduction_contacts_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="nContactsConduction";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="nContactsConduction";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_n_conduction_contacts_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_wall_heattransfer_coeff_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("wallHeattransferCoeff","property/atom","scalar",0,0,this->style,false));
  if(!fix_wall_heattransfer_coeff_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="wallHeattransferCoeff";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="wallHeattransferCoeff";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_wall_heattransfer_coeff_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  fix_wall_temperature_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("wallTemp","property/atom","scalar",0,0,this->style,false));
  if(!fix_wall_temperature_ && store_contact_data_)
  {
    const char* fixarg[10];
    fixarg[0]="wallTemp";
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]="wallTemp";
    fixarg[4]="scalar";
    fixarg[5]="no";
    fixarg[6]="yes";
    fixarg[7]="no";
    fixarg[8]="0.";
    fix_wall_temperature_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  if(store_contact_data_ && (!fix_conduction_contact_area_ || !fix_n_conduction_contacts_ || !fix_wall_heattransfer_coeff_ || !fix_wall_temperature_))
    error->one(FLERR,"internal error");
}

/* ---------------------------------------------------------------------- */

int FixHeatGranCond::setmask()
{
  int mask = FixHeatGran::setmask();
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::updatePtrs()
{
  FixHeatGran::updatePtrs();

  if(store_contact_data_)
  {
    conduction_contact_area_ = fix_conduction_contact_area_->vector_atom;
    n_conduction_contacts_ = fix_n_conduction_contacts_->vector_atom;
    wall_heattransfer_coeff_ = fix_wall_heattransfer_coeff_->vector_atom;
    wall_temp_ = fix_wall_temperature_->vector_atom;
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::init()
{
  // initialize base class
  FixHeatGran::init();

  heatPhy.init();

  if(temp_max_ > 0.)
  {
      const char * newarg[2];
      char * newarg1 = new char[30];
      sprintf(newarg1,"%f",temp_max_);
      newarg[0] = "quantity_max";
      newarg[1] = newarg1;
      fix_ste->modify_param(2, const_cast<char**>(newarg) );
      delete [] newarg1;
  }

  updatePtrs();

  // error checks on coarsegraining
  
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::pre_force(int vflag)
{
    
    if(store_contact_data_)
    {
        fix_wall_heattransfer_coeff_->set_all(0.);
        fix_wall_temperature_->set_all(0.);
    }
}

/* ---------------------------------------------------------------------- */

void FixHeatGranCond::post_force(int vflag)
{

    if (history_flag == 0)
        post_force_eval<0>(vflag);
    else
        post_force_eval<1>(vflag);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG>
void FixHeatGranCond::post_force_eval(int vflag)
{
  double delx,dely,delz;
  double radj,radsum,rsq,r;

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0)
    error->warning(FLERR,"Fix heat/gran/conduction implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0)
    error->warning(FLERR,"Fix heat/gran/conduction implementation may not be valid for pair style hybrid/overlay");

  const int inum = pair_gran->list->inum;
  int *ilist = pair_gran->list->ilist;
  int *numneigh = pair_gran->list->numneigh;
  int **firstneigh = pair_gran->list->firstneigh;

  NeighList *listgranhistory = NULL;
  int *contact_flag = NULL,**first_contact_flag = NULL;
  double *all_contact_hist = NULL,**first_contact_hist = NULL;
  int dnum = 0;
  if(HISTFLAG) {
      listgranhistory = pair_gran->listgranhistory;
      first_contact_flag = listgranhistory->firstneigh;
      first_contact_hist = listgranhistory->firstdouble;
      dnum = listgranhistory->dnum;
  }

  double *radius = atom->radius;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs();

  if(store_contact_data_)
  {
    fix_conduction_contact_area_->set_all(0.);
    fix_n_conduction_contacts_->set_all(0.);
  }

  heatPhy.beginPass();
  updateRequestedCallbacks();
  callBeginPass();

  // loop over neighbors of my atoms
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const double radi = radius[i];
    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    if(HISTFLAG) {
        contact_flag = first_contact_flag[i];
        all_contact_hist = first_contact_hist[i];
    }

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      if(!HISTFLAG)
      {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj;
      }

      if ((HISTFLAG && contact_flag[jj]) || (!HISTFLAG && (rsq < radsum*radsum))) {  //contact
        
        if(HISTFLAG)
        {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
          if(rsq >= radsum*radsum) continue;
        }

        r = sqrt(rsq);

        double contactArea;
        const double hc = heatPhy.computeHtPropertiesPP(i,j,r,radsum,all_contact_hist,dnum*jj,contactArea);

        const double flux = (Temp[j]-Temp[i])*hc;
        const double dirFlux[3] = {flux*delx, flux*dely, flux*delz};

        //Add half of the flux (located at the contact) to each particle in contact
        heatFlux[i] += flux;
        directionalHeatFlux[i][0] += 0.50 * dirFlux[0];
        directionalHeatFlux[i][1] += 0.50 * dirFlux[1];
        directionalHeatFlux[i][2] += 0.50 * dirFlux[2];

        if(store_contact_data_)
        {
          conduction_contact_area_[i] += contactArea;
          n_conduction_contacts_[i] += 1.;
        }
        if (newton_pair || j < nlocal)
        {
          heatFlux[j] -= flux;
          directionalHeatFlux[j][0] += 0.50 * dirFlux[0];
          directionalHeatFlux[j][1] += 0.50 * dirFlux[1];
          directionalHeatFlux[j][2] += 0.50 * dirFlux[2];

          if(store_contact_data_)
          {
            conduction_contact_area_[j] += contactArea;
            n_conduction_contacts_[j] += 1.;
          }
        }

        call_add_heat_pp_callback(i,j,flux);
      }
    }
  }

  callEndPass();

  if(newton_pair)
  {
    fix_heatFlux->do_reverse_comm();
    fix_directionalHeatFlux->do_reverse_comm();
    fix_conduction_contact_area_->do_reverse_comm();
    fix_n_conduction_contacts_->do_reverse_comm();
  }

  if(store_contact_data_)
  for(int i = 0; i < nlocal; i++)
  {
     if(n_conduction_contacts_[i] > 0.5)
        conduction_contact_area_[i] /= n_conduction_contacts_[i];
  }
}
