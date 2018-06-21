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

#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "compute_crosssection.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "region.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "modified_andrew.h"
#include "math_extra_liggghts.h"

#if defined(_WIN32) || defined(_WIN64)
double inline round(double d);
#endif

using namespace LAMMPS_NS;
using MODIFIED_ANDREW_AUX::Circle;

/* ---------------------------------------------------------------------- */

ComputeCrosssection::ComputeCrosssection(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  ComputeContactAtom(lmp, iarg, narg, arg),
  angle_(0),
  file_(0),
  fileflag_(false),
  iregion_(-1),
  idregion_(0)
{
  if (narg < iarg+10)
    error->compute_error(FLERR,this,"need at least 10 args");

  if(strcmp(arg[iarg++],"dim"))
    error->compute_error(FLERR,this,"expecting keyword 'dim'");

  if(strcmp(arg[iarg],"x") == 0)
    dim_ = 0;
  else if(strcmp(arg[iarg],"y") == 0)
    dim_ = 1;
  else if(strcmp(arg[iarg],"z") == 0)
    dim_ = 2;
  else
    error->compute_error(FLERR,this,"expecting 'x', 'y' or 'z' after 'dim'");
  iarg++;

  if(strcmp(arg[iarg++],"min"))
    error->compute_error(FLERR,this,"expecting keyword 'min'");
  min_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"max"))
    error->compute_error(FLERR,this,"expecting keyword 'max'");
  max_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"n_cuts"))
    error->compute_error(FLERR,this,"expecting keyword 'n_cuts'");
  n_cuts_ = atoi(arg[iarg++]);
  if(n_cuts_ < 2)
    error->compute_error(FLERR,this,"'n_cuts' >= 2 required");

  if(strcmp(arg[iarg++],"cut_thickness"))
    error->compute_error(FLERR,this,"expecting keyword 'cut_thickness'");
  cut_thickness_half_ = atof(arg[iarg++]) / 2.;

  // parse args
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
      hasargs = false;
      if(strcmp(arg[iarg],"file") == 0) {
          if(narg < iarg+2)
              error->compute_error(FLERR,this,"not enough arguments for 'file'");
          if(0 == comm->me)
          {
            file_ = fopen(arg[iarg+1],"w");
            if(!file_)
                error->one(FLERR,"Illegal compute crosssection command, cannot open file");
          }
          fileflag_ = true;
          iarg += 2;
          hasargs = true;
      } else if(strcmp(arg[iarg],"region") == 0) {
            if(narg < iarg+2)
                error->compute_error(FLERR,this,"not enough arguments for 'region'");
            int iregion = domain->find_region(arg[iarg+1]);
            if (iregion == -1)
                error->compute_error(FLERR,this,"region ID does not exist");
            int n = strlen(arg[iarg+1]) + 1;
            idregion_ = new char[n];
            strcpy(idregion_,arg[iarg+1]);
            iarg += 2;
            hasargs = true;
      }
      else error->compute_error(FLERR,this,"unknown keyword");
  }

  // calculate slices
  setup_cuts();

  vector_flag = 1;

  if(fileflag_)
    size_vector = n_cuts_+1;
  else
    size_vector = n_cuts_;

  vector = new double[n_cuts_+1];
}

/* ---------------------------------------------------------------------- */

ComputeCrosssection::~ComputeCrosssection()
{
    for(int i = 0; i < n_cuts_; i++)
        delete mod_andrew_[i];
    delete [] mod_andrew_;

    delete [] vector;

    if(file_) fclose(file_);

    if(idregion_) delete []idregion_;
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::init()
{
  ComputeContactAtom::init();

  // set index and check validity of region

  if (idregion_) {
    iregion_ = domain->find_region(idregion_);
    if (iregion_ == -1)
      error->compute_error(FLERR,this,"region ID does not exist");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::setup_cuts()
{
    cut_dist_ = (max_ - min_) / (n_cuts_-1);

    mod_andrew_ = new ModifiedAndrew*[n_cuts_];
    for(int i = 0; i < n_cuts_; i++)
        mod_andrew_[i] = new ModifiedAndrew(lmp);
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::compute_vector()
{
  if(invoked_vector == update->ntimestep)
    return;

  invoked_vector = update->ntimestep;

  compute_peratom();

  for(int i = 0; i < n_cuts_; i++)
  {
      
      vector[i] = mod_andrew_[i]->area();
  }

  if(fileflag_)
  {
      vector[n_cuts_] = calc_ang();
      //fprintf(file_,"%f\n",ang);

      // only proc 0 writes
      if(file_)
        write();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  ComputeContactAtom::compute_peratom();

  compute_convex_hull();
}

/* ---------------------------------------------------------------------- */

inline int ComputeCrosssection::mod(double coo)
{
    double coo_shift = coo - min_;
    int result = round( coo_shift / cut_dist_);
    double remainder = coo_shift - static_cast<double>( result ) * cut_dist_;

    if(remainder > -cut_thickness_half_ && remainder < cut_thickness_half_)
        return result;
    return -1;
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::compute_convex_hull()
{
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double *radius = atom->radius;
    double coo;
    int m;
    Circle circle;

    double mi = min_ - cut_thickness_half_;
    double ma = max_ + cut_thickness_half_;

    for(int i = 0; i < nlocal; i++)
    {
        coo = x[i][dim_];
        if(contact[i] >= 2 && coo > mi && coo < ma)
        {
           if (iregion_ >= 0 && !domain->regions[iregion_]->match(x[i][0],x[i][1],x[i][2]))
             continue;

            m = mod(coo);
            if(m >= 0)
            {
                circle.x = x[i][(1+dim_)%3];
                circle.y = x[i][(2+dim_)%3];
                circle.r = radius[i];
                mod_andrew_[m]->add_contact(circle);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::write()
{
    
    for(int i = 0; i < n_cuts_; i++)
    {
        double coo = min_ + static_cast<double>(i)*cut_dist_;
        fprintf(file_,"%f %f %f\n",coo,vector[i],sqrt(vector[i]/M_PI));
    }
    fflush(file_);
}

/* ---------------------------------------------------------------------- */

double ComputeCrosssection::calc_ang()
{
    
    int ilo = 0;
    int ihi = n_cuts_-1;
    int imid, ilomid, ihimid, idelta;
    double rhimid, rlomid, del, ang;

    const double threshold = cut_thickness_half_*cut_thickness_half_*1e-12;
    for(int i = 0; i < n_cuts_; i++)
    {
        if(MathExtraLiggghts::compDouble(vector[i],0.,threshold))
            ilo += 1;
        else
            break;
    }

    for(int i = n_cuts_-1; i >= 0; i--)
    {
        if(MathExtraLiggghts::compDouble(vector[i],0.,threshold))
            ihi -= 1;
        else
            break;
    }

    if(ilo == ihi)
        error->one(FLERR,"Compute crossection could not calculate angle (1)");

    for(int i = ilo; i <= ihi; i++)
    {
        if(MathExtraLiggghts::compDouble(vector[i],0.,threshold))
            error->one(FLERR,"Compute crossection could not calculate angle - internal error");
    }

    imid = static_cast<int>(0.5*static_cast<double>(ilo+ihi));
    idelta = static_cast<int>(0.25*static_cast<double>(ihi-ilo));

    ilomid = imid-idelta;
    ihimid = imid+idelta;

    if(ilomid == ihimid || ilomid == ilo || ihimid == ihi)
    {
        //ilomid = ilo;
        //ihimid = ihi;
        
        error->one(FLERR,"Compute crossection could not calculate angle (2)");
    }

    rlomid = sqrt(vector[ilomid]/M_PI);
    rhimid = sqrt(vector[ihimid]/M_PI);
    del = (ihimid-ilomid)*cut_dist_;

    if(rhimid > rlomid)
        ang = 90. + 180./M_PI*atan((rhimid-rlomid)/del);
    else
        ang = 180./M_PI*atan(del/(rlomid-rhimid));

    return ang;
}
