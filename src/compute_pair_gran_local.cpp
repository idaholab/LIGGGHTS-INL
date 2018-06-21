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

#include "compute_pair_gran_local.h"
#include "atom.h"
#include "error.h"
#include "fix_heat_gran_conduction.h"
#include "fix_multisphere.h"
#include "fix_wall_gran.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "pair_gran.h"
#include "pair_gran_proxy.h"
#include "update.h"
#include "vector_liggghts.h"

#include <cmath>

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::ComputePairGranLocal(LAMMPS *lmp, int &iarg, int narg, char **arg) :
    Compute(lmp, iarg, narg, arg),
    wall(0), // default: pair data
    reference_exists(0),
    pairgran(NULL),
    fixheat(NULL),
    fixwall(NULL),
    fix_ms(NULL),
    verbose(false),
    dnum(0),
    nmax(0),
    vector(NULL),
    array(NULL)
{
    if (narg < iarg) error->all(FLERR,"Illegal compute pair/gran/local or wall/gran/local command");

    local_flag = 1;
    timeflag = 1;

    // store everything by default expect heat flux
    posflag = velflag = idflag = fflag = torqueflag = histflag = areaflag = 1;

    // do not store fn, ft, heat flux, delta by default
    fnflag = ftflag  = torquenflag = torquetflag = deltaflag = heatflag = cpflag = msidflag = 0;

    //no extra distance for building the list of pairs

    // if further args, store only the properties that are listed
    if(narg > iarg)
        posflag = velflag = idflag = fflag = fnflag = ftflag = torqueflag = torquenflag = torquetflag = histflag = areaflag = deltaflag = heatflag = cpflag = msidflag = 0;

    for (; iarg < narg; iarg++)
    {
        //int i = iarg-3;
        if (strcmp(arg[iarg],"pos") == 0) posflag = 1;
        else if (strcmp(arg[iarg],"vel") == 0) velflag = 1;
        else if (strcmp(arg[iarg],"id") == 0) idflag = 1;
        else if (strcmp(arg[iarg],"force") == 0) fflag = 1;
        else if (strcmp(arg[iarg],"force_normal") == 0) fnflag = 1;
        else if (strcmp(arg[iarg],"force_tangential") == 0) ftflag = 1;
        else if (strcmp(arg[iarg],"torque") == 0) torqueflag = 1;
        else if (strcmp(arg[iarg],"torque_normal") == 0) torquenflag = 1;
        else if (strcmp(arg[iarg],"torque_tangential") == 0) torquetflag = 1;
        else if (strcmp(arg[iarg],"history") == 0) histflag = 1;
        else if (strcmp(arg[iarg],"contactArea") == 0) areaflag = 1;
        else if (strcmp(arg[iarg],"delta") == 0) deltaflag = 1;
        else if (strcmp(arg[iarg],"heatFlux") == 0) heatflag = 1;
        else if (strcmp(arg[iarg],"contactPoint") == 0) cpflag = 1;
        else if (strcmp(arg[iarg],"ms_id") == 0) msidflag = 1;
        else if (strcmp(arg[iarg],"verbose") == 0) verbose = true;
        else if (strcmp(arg[iarg],"extraSurfDistance") == 0) error->all(FLERR,"this keyword is deprecated; neighbor->contactDistanceFactor is now used directly");
        else if(0 == strcmp(style,"wall/gran/local") || 0 == strcmp(style,"pair/gran/local"))
            error->compute_error(FLERR,this,"illegal/unrecognized keyword");
    }

    if(update->ntimestep > 0 && !modify->fix_restart_in_progress())
        error->compute_error(FLERR,this,"Need to define this compute before first run");

}

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::~ComputePairGranLocal()
{
    memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::post_create()
{
    // check if wall data requested
    if(strcmp(style,"wall/gran/local") == 0) wall = 1;

    // initialize once as dump parses in constructor for length of per-atom data
    
    init_cpgl(false);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::pre_delete(bool uncomputeflag)
{
    if(uncomputeflag)
    {
        if(wall == 0) pairgran->unregister_post_force_callback(this);
        else fixwall->unregister_post_force_callback(this);
        if(fixheat) fixheat->unregister_post_force_callback(this);
    }
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init()
{
    init_cpgl(false); 
    
    fix_ms = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init_cpgl(bool requestflag)
{
    newton_pair = force->newton_pair;

    if (idflag && atom->tag_enable == 0)
        error->all(FLERR,"Compute pair/gran/local (or derived) requested to compute IDs, this requires atoms have IDs.");

    // if available from previous run, unregister
    if(pairgran && reference_exists)
        pairgran->unregister_post_force_callback(this);

    if(fixwall && reference_exists)
        fixwall->unregister_post_force_callback(this);
    if(fixheat && reference_exists)
          fixheat->unregister_post_force_callback(this);

    fixwall = NULL;
    fixheat = NULL;
    pairgran = NULL;

    // register heat transfer fix for pair data
    if(wall == 0)
    {
        if (force->pair == NULL)
            error->all(FLERR,"No pair style is defined for compute pair/gran/local");

        pairgran = (PairGran*)force->pair_match("gran",0);

        if (pairgran == NULL)
            error->all(FLERR,"No valid granular pair style found for use with compute pair/gran/local");

        if (pairgran->cplenable() == 0)
            error->all(FLERR,"Pair style does not support compute pair/gran/local");

        pairgran->register_post_force_callback(this);
        dnum = pairgran->dnum_pair();

        if(requestflag)
        {
            
            int irequest = neighbor->request((void *) this);
            neighbor->requests[irequest]->pair = 0;
            neighbor->requests[irequest]->compute = 1;
            neighbor->requests[irequest]->half = 0;
            neighbor->requests[irequest]->gran = 1;
            //neighbor->requests[irequest]->granhistory = 1;
            neighbor->requests[irequest]->occasional = 1;
        }

        // register heat transfer fix if applicable
        if(heatflag)
        {
            for(int ifix = 0; ifix < modify->nfix; ifix++)
            {
                if(
                        (0 == strcmp(modify->fix[ifix]->style,"heat/gran")) ||
                        (0 == strcmp(modify->fix[ifix]->style,"heat/gran/conduction"))
                        )
                {
                    fixheat = static_cast<FixHeatGranCond*>(modify->fix[ifix]);
                }
            }
            if(!fixheat) error->all(FLERR,"Compute pair/gran/local can not calculate heat flux values since no fix heat/gran/conduction not compute them");

            // group of this compute and heat transfer fix must be same so same number of pairs is computed
            if(groupbit != fixheat->groupbit) error->all(FLERR,"Compute pair/gran/local group and fix heat/gran/conduction group cannot be different");
            fixheat->register_post_force_callback(this);
        }
    }
    // register with granular wall, only accept mesh walls
    else
    {
        const int n_wall_fixes = modify->n_fixes_style("wall/gran");

        bool printed_primitive_warning = false;
        for(int ifix = 0; ifix < n_wall_fixes; ifix++)
        {
            FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
            if(fwg->is_mesh_wall())
                fixwall = fwg;
            else if (!printed_primitive_warning)
            {
                error->warning(FLERR, "Compute wall/gran/local detected primitive wall, output will only happen for mesh wall - particle pairs");
                printed_primitive_warning = true;
            }
        }

        if(!fixwall) error->all(FLERR,"Compute wall/gran/local requires a fix of type wall/gran using one or more mesh walls. This fix has to come before the compute in the script");
        fixwall->register_post_force_callback(this);
        dnum = fixwall->dnum();
    }

    // at this point we know that ptr is valid
    reference_exists = 1;

    if(histflag && dnum == 0) error->all(FLERR,"Compute pair/gran/local or wall/gran/local can not calculate history values since pair or wall style does not compute them");
    // standard values: pos1,pos2,id1,id2,extra id for mesh wall,force,torque,contact area

    nvalues = posflag*6 + velflag*6 + idflag*3 + fflag*3 + fnflag*3 + ftflag*3 + torqueflag*3 + torquenflag*3 + torquetflag*3 + histflag*dnum + areaflag + deltaflag + heatflag + cpflag*3 + msidflag*2;
    size_local_cols = nvalues;

}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init_list(int id, NeighList *ptr)
{
    list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::deleteReference()
{
    reference_exists = 0;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::compute_local()
{
    //  invoked_local = update->ntimestep; // moved this to actual calculation -> beginPass

    if(!reference_exists) error->one(FLERR,"Compute pair/gran/local or wall/gran/local reference does no longer exist (pair or fix deleted)");

    // check if compute was calculated in this timestep
    // otherwise it's now too late
    if (invoked_local != update->ntimestep)
        error->all(FLERR,"Compute is not current!");
}

/* ----------------------------------------------------------------------
   count pairs on this proc
------------------------------------------------------------------------- */

int ComputePairGranLocal::count_pairs(int & nCountSurfacesIntersect)
{
    double **x = atom->x;
    double *radius = atom->radius;
    int *mask = atom->mask;
    const int nlocal = atom->nlocal;
    const int nall = nlocal + atom->nghost;

    //neighbor->build_one(list->index);

    list = pairgran->list;
    const int inum = list->inum;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    const double contactDistanceFactorSqr = (neighbor->contactDistanceFactor+1e-16) * (neighbor->contactDistanceFactor+1e-16);

    // loop over neighbors of my atoms
    // skip if I or J are not in group

    int m, n;
    m = n = nCountSurfacesIntersect = 0;
    for (int ii = 0; ii < inum; ii++)
    {
        const int i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        int *jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; jj++)
        {
            int j = jlist[jj];

            if (j >= nall) j %= nall;

            if (!(mask[j] & groupbit)) continue;
            if (newton_pair == 0 && j >= nlocal && atom->tag[i] <= atom->tag[j]) continue;

            const double delx = xtmp - x[j][0];
            const double dely = ytmp - x[j][1];
            const double delz = ztmp - x[j][2];
            const double rsq = delx*delx + dely*dely + delz*delz;
            
            if (rsq <  (radius[i]+radius[j])*(radius[i]+radius[j]))
                nCountSurfacesIntersect++;
            if (rsq > (radius[i]+radius[j])*(radius[i]+radius[j])*contactDistanceFactorSqr)
                continue;
            m++;
        }
    }
    if(verbose)
        printf("ComputePairGranLocal::count_pairs: detected %d pairs , and %d pairs with surface intersection. \n",
               m, nCountSurfacesIntersect);

    return m;
}

/* ----------------------------------------------------------------------
   count wall contacts on this proc
------------------------------------------------------------------------- */

int ComputePairGranLocal::count_wallcontacts(int & nCountWithOverlap)
{
    
    /* during setup the number of partners (npartner_) may not set correctly,
         * since this is done during the first force loop
         * Perform a first estimation on the number of neighbors
         */
    nCountWithOverlap = 0;
    return fixwall->n_neighs_local();
}

/* ----------------------------------------------------------------------
   add data from particle-particle contact on this proc
------------------------------------------------------------------------- */

void ComputePairGranLocal::add_heat_pp(const int i, const int j,const double hf)
{
    if (newton_pair == 0 && j >= atom->nlocal && atom->tag[i] <= atom->tag[j])
        return;

    if (!(atom->mask[i] & groupbit)) return;
    if (!(atom->mask[j] & groupbit)) return;

    if(heatflag)
    {
        // heat flux is always last value
        array[ipair][nvalues-1] = hf;
        
    }
    else
        error->one(FLERR,"Illegal situation in ComputePairGranLocal::add_heat");

    // inc counter again, since reset in compute_local() after getting pair data
    ipair++;
}

/* ----------------------------------------------------------------------
   add data from particle-wall contact on this proc
------------------------------------------------------------------------- */

void ComputePairGranLocal::add_heat_pw(const LCM::SurfacesIntersectData &sidata,double hf)
{
    const int ip = sidata.i;

    if (!(atom->mask[ip] & groupbit))
        return;

    if(heatflag)
    {
        // heat flux is always last value
        // use ipair -1 , add_heat_wall is not always called
        array[ipair-1][nvalues-1] = hf;
    }
    //else error->one(FLERR,"Illegal situation in ComputePairGranLocal::add_heat_wall");
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::pair_finalize()
{
    // ensure have the right number of entries - i.e. the # of sidata.hasforceupdate occurrancies!

    ncount_added_via_pair = ipair;
    size_local_rows = ncount_added_via_pair;

    //    fprintf(screen,"%s pair_finalize() ipair %d\n",id,ipair);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::reallocate(int n)
{
    // grow vector or array and indices array

    while (nmax < n)
        nmax += DELTA;

    memory->destroy(array);
    memory->create(array,nmax,nvalues,"pair/local:array");
    array_local = array;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePairGranLocal::memory_usage()
{
    double bytes = nmax*nvalues * sizeof(double);
    return bytes;
}

int ComputePairGranLocal::get_history_offset(const char * const name)
{
    if (pairgran)
        return static_cast<PairGranProxy*>(pairgran)->get_history_offset(name);
    else if (fixwall)
        return fixwall->get_history_offset(name);
    else
    {
        error->all(FLERR, "Internal error");
        return -1;
    }
}

void ComputePairGranLocal::beginPass()
{
    
    if(!reference_exists)
        error->one(FLERR,"Compute pair/gran/local or wall/gran/local reference does no longer exist (pair or fix deleted)");

    // do nothing if this function was already once called this timestep
    if (invoked_local == update->ntimestep)
        return;

    // count local entries and compute pair info

    int nCountSurfacesIntersect(0);
    if(wall == 0)
        ncount = count_pairs(nCountSurfacesIntersect);            // # pairs is ensured to be the same for pair and heat
    else
        ncount = count_wallcontacts(nCountSurfacesIntersect);     // # wall contacts ensured to be same for wall/gran and heat

    // only consider rows with overlap (but allocate memory for all)
    if (ncount > nmax)
        reallocate(ncount);
    size_local_rows = nCountSurfacesIntersect; 

    ipair = 0; //reset
}

void ComputePairGranLocal::compute_post_force(const LCM::SurfacesIntersectData &sidata, const LCM::ForceData &i_forces, const LCM::ForceData &j_forces)
{
    
    if (sidata.is_wall) 
        post_force_pw(sidata,i_forces,j_forces);
    else
        post_force_pp(sidata,i_forces,j_forces);
}

void ComputePairGranLocal::endPass()
{
    pair_finalize();
    // update after calculation
    invoked_local = update->ntimestep;
}

void ComputePairGranLocal::post_force_pp(const LCM::SurfacesIntersectData &sidata, const LCM::ForceData &i_forces, const LCM::ForceData &j_forces)
{

    const int i = sidata.i;
    const int j = sidata.j;

    if (!(atom->mask[i] & groupbit))
        return;
    if (!(atom->mask[j] & groupbit))
        return;

    const int nlocal = atom->nlocal;
    if (newton_pair == 0 && j >= nlocal && atom->tag[i] <= atom->tag[j])
        return;

    double *contact_pos; // unused for pairs
    if(!decide_add(sidata.contact_history, contact_pos))
        return;

    const double *xi = atom->x[i];
    const double *xj = atom->x[j];
    const double *vi = atom->v[i];
    const double *vj = atom->v[j];

    const double fx = i_forces.delta_F[0];
    const double fy = i_forces.delta_F[1];
    const double fz = i_forces.delta_F[2];
    const double tor1 = i_forces.delta_torque[0];
    const double tor2 = i_forces.delta_torque[1];
    const double tor3 = i_forces.delta_torque[2];

    if(ipair>=nmax)
    {
        if(screen) fprintf(screen,"ipair: %d, ncount: %d. \n", ipair, ncount);
        error->one(FLERR,"Attempt to add_pair, but number of pairs (nmax) is too small. Try to increase contact_distance_factor via an appropriate neigh_modify command!");
    }

    int n = 0;
    if(posflag)
    {
        double xi_w[3],xj_w[3];
        vectorCopy3D(xi,xi_w);
        vectorCopy3D(xj,xj_w);
        domain->remap(xi_w);
        domain->remap(xj_w);
        vectorToBuf3D(xi_w,array[ipair],n);
        vectorToBuf3D(xj_w,array[ipair],n);
    }
    if(velflag)
    {
        vectorToBuf3D(vi,array[ipair],n);
        vectorToBuf3D(vj,array[ipair],n);
    }
    if(idflag)
    {
        array[ipair][n++] = static_cast<double>(atom->tag[i]);
        array[ipair][n++] = static_cast<double>(atom->tag[j]);
        if(i < nlocal && j < nlocal)
            array[ipair][n++] = 0.;
        else
        {
            if(domain->is_periodic_ghost(i) || domain->is_periodic_ghost(j))
                array[ipair][n++] = 1.;
            else
                array[ipair][n++] = 0.;
        }
    }
    if(fflag)
    {
        array[ipair][n++] = fx;
        array[ipair][n++] = fy;
        array[ipair][n++] = fz;
        
    }
    if(fnflag)
    {
        double fc[3],fn[3],normal[3];

        vectorConstruct3D(fc,fx,fy,fz);
        vectorSubtract3D(atom->x[j],atom->x[i],normal);
        vectorNormalize3D(normal);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);

        array[ipair][n++] = fn[0];
        array[ipair][n++] = fn[1];
        array[ipair][n++] = fn[2];
    }
    if(ftflag)
    {
        double fc[3],fn[3],ft[3],normal[3];

        vectorConstruct3D(fc,fx,fy,fz);
        vectorSubtract3D(atom->x[j],atom->x[i],normal);
        vectorNormalize3D(normal);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);
        vectorSubtract3D(fc,fn,ft);

        array[ipair][n++] = ft[0];
        array[ipair][n++] = ft[1];
        array[ipair][n++] = ft[2];

    }
    if(torqueflag)
    {
        array[ipair][n++] = tor1;
        array[ipair][n++] = tor2;
        array[ipair][n++] = tor3;
    }
    if(torquenflag)
    {
        double torc[3],torn[3],normal[3];

        vectorConstruct3D(torc,tor1,tor2,tor3);
        vectorSubtract3D(atom->x[j],atom->x[i],normal);
        vectorNormalize3D(normal);
        const double tornmag = vectorDot3D(torc,normal);
        vectorScalarMult3D(normal,tornmag,torn);

        array[ipair][n++] = torn[0];
        array[ipair][n++] = torn[1];
        array[ipair][n++] = torn[2];
    }
    if(torquetflag)
    {
        double torc[3],torn[3],tort[3],normal[3];

        vectorConstruct3D(torc,tor1,tor2,tor3);
        vectorSubtract3D(atom->x[j],atom->x[i],normal);
        vectorNormalize3D(normal);
        const double tornmag = vectorDot3D(torc,normal);
        vectorScalarMult3D(normal,tornmag,torn);
        vectorSubtract3D(torc,torn,tort);

        array[ipair][n++] = tort[0];
        array[ipair][n++] = tort[1];
        array[ipair][n++] = tort[2];
    }
    if(histflag)
    {
        for(int d = 0; d < dnum; d++)
            array[ipair][n++] = sidata.contact_history[d];
    }
    if(areaflag)
    {
        const double radi = atom->radius[i];
        const double radj = atom->radius[j];
        double del[3];
        vectorSubtract3D(atom->x[i],atom->x[j],del);
        const double rsq = vectorMag3DSquared(del);
        const double r = sqrt(rsq);
        const double contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/rsq;
        array[ipair][n++] = contactArea;
        
    }
    if(deltaflag)
    {
        if(atom->superquadric_flag)
        {
            int alpha1_offset = get_history_offset("alpha1_offset");
            int alpha2_offset = get_history_offset("alpha2_offset");
            array[ipair][n++] = sidata.contact_history[alpha1_offset]+sidata.contact_history[alpha2_offset];
        }
        else
        {
            const double radi = atom->radius[i];
            const double radj = atom->radius[j];
            double del[3];
            vectorSubtract3D(atom->x[i],atom->x[j],del);
            array[ipair][n++] = radi+radj-vectorMag3D(del);
            
        }
    }

    if(cpflag)
    {
        double cp[3];
        if(atom->superquadric_flag)
        {
            int contact_point_offset = get_history_offset("contact_point_offset");
            cp[0] = sidata.contact_history[contact_point_offset];
            cp[1] = sidata.contact_history[contact_point_offset + 1];
            cp[2] = sidata.contact_history[contact_point_offset + 2];
        }
#ifdef NONSPHERICAL_ACTIVE_FLAG
        else if (atom->shapetype_flag)
        {
            vectorCopy3D(sidata.contact_point, cp);
        }
#endif
        else
        {
            const double radi = atom->radius[i];
            const double radj = atom->radius[j];
            for(int dim = 0; dim < 3; dim++)
                cp[dim] = (atom->x[i][dim]*radj + atom->x[j][dim]*radi)/(radi+radj);
        }

        array[ipair][n++] = cp[0];
        array[ipair][n++] = cp[1];
        array[ipair][n++] = cp[2];
    }
    if(msidflag)
    {
        array[ipair][n++] = static_cast<double>(fix_ms->belongs_to(i));
        array[ipair][n++] = static_cast<double>(fix_ms->belongs_to(j));
    }

    ipair++;
}

void ComputePairGranLocal::post_force_pw(const LCM::SurfacesIntersectData &sidata, const LCM::ForceData &i_forces, const LCM::ForceData &j_forces)
{
    const int i = sidata.i;

    const double fx = i_forces.delta_F[0];
    const double fy = i_forces.delta_F[1];
    const double fz = i_forces.delta_F[2];
    const double tor1 = i_forces.delta_torque[0]*sidata.area_ratio;
    const double tor2 = i_forces.delta_torque[1]*sidata.area_ratio;
    const double tor3 = i_forces.delta_torque[2]*sidata.area_ratio;
    double normal[3];
    vectorCopy3D(sidata.en, normal);
    vectorNegate3D(normal);

    if (!(atom->mask[i] & groupbit))
        return;

    int n = 0;

    double contact_pos[3];
#ifdef NONSPHERICAL_ACTIVE_FLAG
    if(atom->superquadric_flag)
        vectorCopy3D(sidata.contact_point,contact_pos); // default value
    else
#endif
        vectorSubtract3D(atom->x[i],sidata.delta,contact_pos);

    //    if(!decide_add(sidata.contact_history, contact_pos))
    //        return;

    if(posflag)
    {
        array[ipair][n++] = contact_pos[0];
        array[ipair][n++] = contact_pos[1];
        array[ipair][n++] = contact_pos[2];
        array[ipair][n++] = atom->x[i][0];
        array[ipair][n++] = atom->x[i][1];
        array[ipair][n++] = atom->x[i][2];
    }
    if(velflag)
    {
        vectorCopy3D(sidata.v_j,&array[ipair][n]);
        n += 3;
        vectorCopy3D(sidata.v_i,&array[ipair][n]);
        n += 3;
    }
    if(idflag)
    {
        array[ipair][n++] = static_cast<double>(sidata.iMesh);
        array[ipair][n++] = static_cast<double>(sidata.j); // iTri
        array[ipair][n++] = static_cast<double>(atom->tag[i]);
        
    }
    if(fflag)
    {
        array[ipair][n++] = fx;
        array[ipair][n++] = fy;
        array[ipair][n++] = fz;
    }
    if(fnflag)
    {
        double fc[3],fn[3];

        vectorConstruct3D(fc,fx,fy,fz);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);

        array[ipair][n++] = fn[0];
        array[ipair][n++] = fn[1];
        array[ipair][n++] = fn[2];
    }
    if(ftflag)
    {
        double fc[3],fn[3],ft[3];

        vectorConstruct3D(fc,fx,fy,fz);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);
        vectorSubtract3D(fc,fn,ft);

        array[ipair][n++] = ft[0];
        array[ipair][n++] = ft[1];
        array[ipair][n++] = ft[2];
    }
    if(torqueflag)
    {
        array[ipair][n++] = tor1;
        array[ipair][n++] = tor2;
        array[ipair][n++] = tor3;
    }
    if(torquenflag)
    {
        double torc[3],torn[3];

        vectorConstruct3D(torc,tor1,tor2,tor3);
        const double tornmag = vectorDot3D(torc,normal);
        vectorScalarMult3D(normal,tornmag,torn);

        array[ipair][n++] = torn[0];
        array[ipair][n++] = torn[1];
        array[ipair][n++] = torn[2];
    }
    if(torquetflag)
    {
        double torc[3],torn[3],tort[3];

        vectorConstruct3D(torc,tor1,tor2,tor3);
        const double tornmag = vectorDot3D(torc,normal);
        vectorScalarMult3D(normal,tornmag,torn);
        vectorSubtract3D(torc,torn,tort);

        array[ipair][n++] = tort[0];
        array[ipair][n++] = tort[1];
        array[ipair][n++] = tort[2];
    }
    if(histflag)
    {
        for(int d = 0; d < dnum; d++)
            array[ipair][n++] = sidata.contact_history[d];
    }
    if(areaflag)
    {
        const double contactArea = (atom->radius[i]*atom->radius[i]-sidata.rsq)*M_PI;
        array[ipair][n++] = contactArea;
        
    }
    if(deltaflag)
    {
        array[ipair][n++] = atom->radius[i]-sqrt(sidata.rsq);

    }
    if(msidflag)
    {
        n++;
        array[ipair][n++] = static_cast<double>(fix_ms->belongs_to(i));
    }

    ipair++;
}
