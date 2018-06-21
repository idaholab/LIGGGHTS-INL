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
    (if no contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <set>
#include "fix_insert_fragments.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "comm.h"
#include "region.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_multisphere.h"
#include "multisphere_parallel.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_multiplespheres.h"
#include "particleToInsert.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixInsertFragments::FixInsertFragments(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg),
  fix_fragments_(0),
  check_r_bound_(true),
  n_fragments_(0),
  fix_replace_(0),
  idregion_(0),
  ins_region_(0),
  fix_ms_(0),
  ms_(0),
  type_replace_(-1),
  n_replace_(0),n_replace_this_(0),n_replace_this_local_(0),
  mass_replace_(0.),mass_replace_this_(0.),mass_replace_this_local_(0.),
  replacedata_(0)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments for 'region'");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->fix_error(FLERR,this,"Region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion_ = new char[n];
      strcpy(idregion_,arg[iarg+1]);
      ins_region_ = domain->regions[iregion];
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"type_replace") == 0) {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'type_replace'");
      type_replace_ = atoi(arg[iarg+1]);
      if(type_replace_ < 1 || type_replace_ > atom->ntypes)
        error->fix_error(FLERR,this,"1 <= type_replace <= #atom types required");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"check_bound_sphere") == 0) {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'check_bound_sphere'");
      if(0 == strcmp(arg[iarg+1],"yes"))
        check_r_bound_ = true;
      else if(0 == strcmp(arg[iarg+1],"no"))
        check_r_bound_ = false;
      else
        error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'check_bound_sphere'");
      iarg += 2;
      hasargs = true;
    } else error->fix_error(FLERR,this,"unknown keyword");
  }

  // no fixed insertion target (such as # particles) since insertion triggered by break-ups
  ninsert_exists = 0;

  // since fragments are normalized
  warn_boxentent = false;

  // turn off overlap check and turn off start stats since both not required
  check_ol_flag = 0;
  print_stats_start_flag = 0;

  nevery = 1;

  if(maxrad >= 1.)
    error->fix_error(FLERR,this,"Particle distribution must be relative, max radius mus be < 1");

  // check for fix holding particle fragments
  if(fix_distribution->n_particletemplates() != 1)
      error->fix_error(FLERR,this,"fix of type particledistribution/discrete must hold exactly one template");
  if(strcmp(fix_distribution->particletemplates()[0]->style,"particletemplate/multiplespheres"))
      error->fix_error(FLERR,this,"fix of type particledistribution/discrete must hold exactly "
                                    "one template of type fix particletemplate/multiplespheres");
  fix_fragments_ = static_cast<FixTemplateMultiplespheres*>(fix_distribution->particletemplates()[0]);

  if(!fix_fragments_->is_relative())
      error->fix_error(FLERR,this,"fix particletemplate/multiplespheres used by this fix must use "
                                  "option 'relative = yes'");

  // get number of fragments from fix
  n_fragments_ = fix_fragments_->number_spheres();

  if(check_r_bound_ && fix_fragments_->max_r_bound() > 1.)
  {
      char errstr[200];
      sprintf(errstr,"Bounding sphere of particle template %s is %f, but should be 1 at maximum",
                     fix_fragments_->id,fix_fragments_->max_r_bound());
      error->fix_error(FLERR,this,errstr);
  }

  //if(!fix_fragments_->all_overlap_atleast_one_slightly())
  //    error->fix_error(FLERR,this,"fix particletemplate/multiplespheres - all spheres must overlap slightly at least one sphere");
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::post_create()
{
  if(!fix_replace_)
  {
        char *replace_name = new char[9+strlen(id)];
        sprintf(replace_name,"replace_%s",id);

        const char * fixarg[9];
        
        fixarg[0]=replace_name;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=replace_name;
        fixarg[4]="scalar"; 
        fixarg[5]="yes";    
        fixarg[6]="yes";    
        fixarg[7]="yes";    
        fixarg[8]="0.";
        fix_replace_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete [] replace_name;
  }
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::pre_delete(bool unfixflag)
{
    if(unfixflag)
        modify->delete_fix(fix_replace_->id);
}

/* ---------------------------------------------------------------------- */

FixInsertFragments::~FixInsertFragments()
{
    if(replacedata_)
        memory->destroy(replacedata_);
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::init_defaults()
{
    check_r_bound_ = true;
    n_fragments_ = 0;
    fix_fragments_ = 0;
    fix_replace_ = 0;
    idregion_ = 0;
    ins_region_ = 0;
    type_replace_ = -1;
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::init()
{
    FixInsert::init();

    // check validity of region
    if (idregion_) {
      int iregion = domain->find_region(idregion_);
      if (iregion == -1)
        error->fix_error(FLERR,this,"Region ID does not exist");
      ins_region_ = domain->regions[iregion];
    }

    // multisphere support
    fix_ms_ = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
    if(modify->n_fixes_style("multisphere") > 1)
      error->fix_error(FLERR,this,"does not support more than one fix multisphere.");
    if(fix_ms_)
      ms_ = &(fix_ms_->data());

}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
   called via setup
------------------------------------------------------------------------- */

void FixInsertFragments::calc_insertion_properties()
{
    // error checks
    if(nflowrate > 0. || massflowrate > 0.)
        error->fix_error(FLERR,this,"specifying 'nflowrate' or 'massflowrate' is not allowed");
    if(ninsert > 0 || massinsert > 0.)
        error->fix_error(FLERR,this,"specifying 'nparticles' or 'mass' is not allowed");
    if(insert_every <= 0)
        error->fix_error(FLERR,this,"specifying 'every' must be > 0");

    // do not need ninsert_per
}

/* ---------------------------------------------------------------------- */

int FixInsertFragments::setmask()
{
    int mask = FixInsert::setmask();
    mask |= PRE_NEIGHBOR;
    mask |= END_OF_STEP;
    return mask;
}

/* ---------------------------------------------------------------------- */

BoundingBox FixInsertFragments::getBoundingBox() {

    if(ins_region_)
    {
        BoundingBox bb(ins_region_->extent_xlo, ins_region_->extent_xhi,
                 ins_region_->extent_ylo, ins_region_->extent_yhi,
                 ins_region_->extent_zlo, ins_region_->extent_zhi);
        const double extend = 2*neighbor->skin;
        bb.extendByDelta(extend);
        return bb;
    }
    else
    {
        BoundingBox bb(domain->sublo[0], domain->boxhi[0],
                     domain->boxlo[1], domain->boxhi[1],
                     domain->boxlo[2], domain->boxhi[2]);
        const double extend = 2*neighbor->skin;
        bb.extendByDelta(extend);
        return bb;
    }
}

/* ---------------------------------------------------------------------- */

inline int FixInsertFragments::is_nearby(int i)
{
    // need not check overlap with existing particles since we
    // know space originally taken by deleted particles is free
    return 0;
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::pre_neighbor()
{
    
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    int nall = nlocal+nghost;
    double *flag = fix_replace_->vector_atom;

    for(int i = nlocal; i < nall; i++)
    {
        flag[i] = 0.;
    }
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::setup_pre_neighbor()
{
    pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::end_of_step()
{
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    int nall = nlocal+nghost;
    double *flag = fix_replace_->vector_atom;
    int *mask = atom->mask;
    int *type = atom->type;
    double **x = atom->x;
    double xcm[3];

    // check breakage criterion for all local particles

    for(int i = 0; i < nall; i++)
    {
        if (mask[i] & groupbit)
        {
            // skip if type not applicable
            if(type_replace_ > -1 && type_replace_ != type[i])
                continue;

            // case single sphere
            // do not handle multi-sphere
            if (
                 (!fix_ms_ || fix_ms_->belongs_to(i) < 0) &&
                 (!ins_region_ || ins_region_->match(x[i][0],x[i][1],x[i][2]))
               )
            {
                // do not need to care about ghosts for non-MS case
                if(i >= nlocal)
                    continue;
                flag[i] = 1.;
                
            }

            // case multi-sphere
            // owned atoms of not-owned bodies will not be flagged at this point!

            if(fix_ms_ && fix_ms_->belongs_to(i) >= 0)
            {
                int ibody = ms_->map(fix_ms_->belongs_to(i));
                if(ibody >= 0)
                {
                    ms_->xcm(xcm,ibody);

                    if(!ins_region_ || ins_region_->match(xcm[0],xcm[1],xcm[2]))
                    {
                        flag[i] = 1.;
                        
                    }
                }
            }
        }
    }

    // reverse comm of flag so owned atoms of not-owned bodies will be flagged here
    // in the unlikely case anything is double-counted, reset it here
    fix_replace_->do_reverse_comm();

    for(int i = 0; i < nlocal; i++)
        if (flag[i] > (1. - 1e-6))
            flag[i] = 1.;
}

/* ---------------------------------------------------------------------- */

bool FixInsertFragments::pre_insert()
{
    int i,ireplace;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    double **x = atom->x;
    double **v = atom->v;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *flag = fix_replace_->vector_atom;
    AtomVec *avec = atom->avec;
    double xcm[3],vcm[3],quat[4],rbound;
    double unitquat[4];
    quatIdentity4D(unitquat);
    std::set<int> bodytag_list;

    // count # of particles to remove

    bodytag_list.clear();
    n_replace_this_local_ = 0;
    mass_replace_this_local_ = 0.;
    for(i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit && MathExtraLiggghts::compDouble(flag[i],1.,1e-6))
        {
            // case single sphere
            if(!fix_ms_ || fix_ms_->belongs_to(i) < 0)
            {
                n_replace_this_local_++;
                mass_replace_this_local_ += rmass[i];
            }

            // case multisphere
            // count each MS body only if owned on this proc, and only once
            else if(fix_ms_)
            {
                int body_tag = fix_ms_->belongs_to(i) ;
                const bool is_counted = bodytag_list.find(body_tag) != bodytag_list.end();
                int ibody = ms_->map(body_tag);

                // body owned not counted already
                if((ibody >= 0) && !is_counted)
                {
                    bodytag_list.insert(body_tag);
                    n_replace_this_local_++;
                    mass_replace_this_local_ += rmass[i];
                }
            }
        }
    }

   // tally stats

   MPI_Sum_Scalar(n_replace_this_local_,n_replace_this_,world);
   n_replace_ += n_replace_this_;
   MPI_Sum_Scalar(mass_replace_this_local_,mass_replace_this_,world);
   mass_replace_ += mass_replace_this_;

   // allocate breakage data
   
   if(replacedata_)
     memory->destroy(replacedata_);
   memory->create(replacedata_,n_replace_this_local_,ms_?11:7,"FixInsertFragments::replacedata");

   // fill breakage data and remove particles
   // clear bodytags again here to exclude potential double

   bodytag_list.clear();
   i = ireplace = 0;
   while (i < nlocal)
   {
      if (mask[i] & groupbit && MathExtraLiggghts::compDouble(flag[i],1.,1e-6))
      {
          
          // case single sphere
          if(!fix_ms_ || fix_ms_->belongs_to(i) < 0)
          {
              //copy data needed for insertion
              vectorCopy3D(x[i],&replacedata_[ireplace][0]);
              vectorCopy3D(v[i],&replacedata_[ireplace][3]);
              replacedata_[ireplace][6] = radius[i];
              if(ms_) vectorCopy4D(unitquat,&replacedata_[ireplace][7]);
              ireplace++;

              // delete particle
              avec->copy(nlocal-1,i,1);
              nlocal--;
          }
          // case multisphere
          else if(fix_ms_)
          {
                // only add to replacedata if owned body
                int body_tag = fix_ms_->belongs_to(i);
                int ibody = ms_->map(body_tag);
                if(ibody >= 0)
                {
                    const bool is_in = bodytag_list.find(body_tag) != bodytag_list.end();

                    // bodytag_list contains tag, i.e. already processed
                    if(!is_in)
                    {
                        // add tag to set
                        bodytag_list.insert(body_tag);

                        // copy data needed for insertion
                        // IMPORTANT: USING XCM HERE, NOT XBOUND
                        // THIS IS JUST A CONVENTION!!!
                        // but is ok since template is internally transferred
                        // into coordinate system with xcm = 0/0/0

                        ms_->xcm(xcm,ibody);
                        ms_->vcm(vcm,ibody);
                        ms_->quat(quat,ibody);
                        rbound = ms_->r_bound(ibody);
                        
                        vectorCopy3D(xcm,&replacedata_[ireplace][0]);
                        vectorCopy3D(vcm,&replacedata_[ireplace][3]);
                        replacedata_[ireplace][6] = rbound;
                        vectorCopy4D(quat,&replacedata_[ireplace][7]);
                        ireplace++;
                    }
                }

                // delete particle in any case, also if it belongs to a body which is not owned
                avec->copy(nlocal-1,i,1);
                nlocal--;
          }
      }
      else i++;
   }

   // update local and global # particles
   
   atom->nlocal = nlocal;
   double rlocal = static_cast<double>(atom->nlocal);
   MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);

   // print stats
   print_stats_replace_during();

   return true;
}

/* ---------------------------------------------------------------------- */

int FixInsertFragments::calc_ninsert_this()
{
    // number of ptis to insert this timestep
    // will effectively insert n_break_this * n_fragments spheres
    return n_replace_this_;
}

/* ---------------------------------------------------------------------- */

int FixInsertFragments::distribute_ninsert_this(int)
{
    
    return n_replace_this_local_;
}

/* ---------------------------------------------------------------------- */

void FixInsertFragments::print_stats_replace_during()
{
  bigint step = update->ntimestep;

  if (me == 0 && n_replace_this_ > 0)
  {
    if (screen)
      fprintf(screen ,"Particle replacement: replaced %d particle templates (mass %f) at step " BIGINT_FORMAT "\n "
                      " - a total of %d particles (mass %f) replaced so far.\n",
              n_replace_this_,mass_replace_this_,step,n_replace_,mass_replace_);

    if (logfile)
      fprintf(logfile,"Particle replacement: replaced %d particle templates (mass %f) at step " BIGINT_FORMAT "\n "
                      " - a total of %d particles (mass %f) replaced so far.\n",
              n_replace_this_,mass_replace_this_,step,n_replace_,mass_replace_);
  }
}

/* ----------------------------------------------------------------------
   generate new particles at positions where old particles were deleted
   function is executed locally on each process

   overlap check is not needed since space around broken particles is empty

   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertFragments::x_v_omega(int ninsert_this_local,int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
    double pos_ins[3],v_ins[3],omega_ins[3],quat_ins[4],rad_replaced;
    int iparticle, nins;
    ParticleToInsert *pti;

    vectorZeroize3D(omega_ins);
    quatIdentity4D(quat_ins);

    // local insertion
    ninserted_this_local = ninserted_spheres_this_local = 0;
    mass_inserted_this_local = 0.;

    // n_break_this ptis with n_fragments spheres each
    // n_break_this * n_fragments spheres to be generated

    iparticle = 0;

    while(iparticle < n_replace_this_local_)
    {
        
        // get position, velocity and radius of broken particle
        vectorCopy3D(&replacedata_[iparticle][0],pos_ins);
        vectorCopy3D(&replacedata_[iparticle][3],v_ins);
        // omega neglected
        rad_replaced = replacedata_[iparticle][6];

        if(ms_) vectorCopy4D(&replacedata_[iparticle][7],quat_ins);

        // get pti and scale it down with radius of broken particle
        pti = fix_distribution->pti_list[iparticle];
        pti->scale_pti(rad_replaced);

        nins = pti->set_x_v_omega(pos_ins,v_ins,omega_ins,quat_ins);

        // tally stats
        ninserted_spheres_this_local += nins;
        mass_inserted_this_local += pti->mass_ins;
        ninserted_this_local++;

        iparticle++;
    }

    // tally stats, have to do this since operation is locally on each process
    // as opposed to e.g. FixInsertPack::x_v_omega()

}
